"""关键词驱动的酶搜索 + DLKcat 流水线 CLI。"""

from __future__ import annotations

import argparse
import logging
import re
import shutil
import sys
from pathlib import Path
from typing import Any

from .alphafold_db import collect_structure_accessions, download_structures_for_fasta
from .blast_runner import run_blast
from .dedupe import merge_fastas, run_cdhit
from .dlkcat_runner import run_dlkcat_prediction
from .foldseek_runner import run_foldseek_all_vs_all, run_foldseek_seq_vs_afdb
from .io_utils import (
    ROOT,
    accession_from_header,
    atomic_write_text,
    ensure_dir,
    load_config,
    parse_fasta,
    resolve_path,
    resolve_threads,
    result_fail,
    result_ok,
    setup_logging,
    write_fasta,
)
from .substrates import resolve_substrates
from .uniprot_client import download_seed

logger = logging.getLogger(__name__)

STEPS = (
    "uniprot",
    "seed_dedupe",
    "blast1",
    "afdb",
    "foldseek_afdb",
    "foldseek1",
    "dedupe",
    "blast2",
    "foldseek2",
    "finalize",
    "dlkcat",
)


def _safe_keyword(keyword: str) -> str:
    s = re.sub(r"[^A-Za-z0-9._-]+", "_", (keyword or "default").strip())
    return s.strip("_") or "default"


def query_from_keyword(keyword: str) -> str:
    kw = keyword.strip()
    return f"(gene:{kw}) OR (protein_name:{kw}) OR {kw}"


def apply_keyword_to_cfg(cfg: dict[str, Any], keyword: str | None) -> dict[str, Any]:
    """注入 keyword、自动 query，并设置分目录路径。"""
    kw = _safe_keyword(keyword or str(cfg.get("keyword") or "default"))
    cfg["keyword"] = kw
    cfg["_project_root"] = str(ROOT)
    if not str(cfg.get("query") or "").strip():
        cfg["query"] = query_from_keyword(kw)

    paths = dict(cfg.get("paths") or {})
    data_root = resolve_path(paths.get("data_root", "data"))
    results_root = resolve_path(paths.get("results_root", "results"))
    paths["data_raw"] = str(data_root / kw / "raw")
    paths["structures"] = str(data_root / kw / "structures")
    paths["results"] = str(results_root / kw)
    paths["substrates"] = str(resolve_path(paths.get("substrates", "data/substrates")))
    cfg["paths"] = paths
    return cfg


def _cfg_paths(cfg: dict[str, Any]) -> dict[str, Path]:
    paths = cfg.get("paths") or {}
    return {
        "raw": resolve_path(paths.get("data_raw", "data/raw")),
        "structures": resolve_path(paths.get("structures", "data/structures")),
        "results": resolve_path(paths.get("results", "results")),
        "substrates": resolve_path(paths.get("substrates", "data/substrates")),
    }


def _seed_prefix(cfg: dict[str, Any]) -> str:
    return _safe_keyword(str(cfg.get("keyword") or "seed"))


def _seed_fasta(cfg: dict[str, Any]) -> Path:
    return _cfg_paths(cfg)["raw"] / f"{_seed_prefix(cfg)}_seed.fasta"


def _seed_unique_fasta(cfg: dict[str, Any]) -> Path:
    return _cfg_paths(cfg)["raw"] / f"{_seed_prefix(cfg)}_seed_unique.fasta"


def _seed_tsv(cfg: dict[str, Any]) -> Path:
    return _cfg_paths(cfg)["raw"] / f"{_seed_prefix(cfg)}_uniprot.tsv"


def _blast_query_fasta(cfg: dict[str, Any]) -> Path:
    """优先用种子 CD-HIT 代表集；不存在则退回原始 seed。"""
    unique = _seed_unique_fasta(cfg)
    if unique.exists() and parse_fasta(unique):
        return unique
    return _seed_fasta(cfg)


def step_uniprot(cfg: dict[str, Any], force: bool) -> dict[str, Any]:
    paths = _cfg_paths(cfg)
    up = cfg.get("uniprot") or {}
    query = cfg.get("query") or query_from_keyword(str(cfg.get("keyword") or "gmas"))
    r = download_seed(
        query,
        paths["raw"],
        force=force,
        timeout=float(up.get("timeout_sec", 120)),
        max_retries=int(up.get("max_retries", 5)),
        page_size=int(up.get("page_size", 500)),
        seed_prefix=_seed_prefix(cfg),
    )
    if not r.get("ok"):
        return r
    if int(r.get("count") or 0) == 0:
        return result_fail("UniProt seed 为空，流水线中止", **r)
    return r


def step_seed_dedupe(cfg: dict[str, Any], force: bool) -> dict[str, Any]:
    """对 UniProt 种子做 CD-HIT，缩小后续 DIAMOND/BLAST 查询量。"""
    seed = _seed_fasta(cfg)
    if not seed.exists() or not parse_fasta(seed):
        return result_fail("种子 FASTA 不存在或为空，无法 seed_dedupe", fasta=str(seed))
    unique = _seed_unique_fasta(cfg)
    cluster_map = _cfg_paths(cfg)["raw"] / f"{_seed_prefix(cfg)}_seed_cluster_map.tsv"
    dcfg = cfg.get("dedupe") or {}
    r = run_cdhit(
        seed,
        unique,
        cluster_map,
        identity=float(dcfg.get("identity", 0.9)),
        word_size=int(dcfg.get("word_size", 5)),
        threads=resolve_threads(dcfg.get("threads", 0)),
        force=force,
    )
    if r.get("ok"):
        n_in = len(parse_fasta(seed))
        n_out = int(r.get("count") or len(parse_fasta(unique)))
        logger.info("种子去重: %s → %s (identity=%s)", n_in, n_out, dcfg.get("identity", 0.9))
        r = {**r, "input_count": n_in, "count": n_out, "unique_fasta": str(unique)}
    return r


def step_blast(
    cfg: dict[str, Any],
    query_fasta: Path,
    out_dir: Path,
    *,
    force: bool,
    source_tag: str,
) -> dict[str, Any]:
    b = cfg.get("blast") or {}
    db = str(b.get("database", "swissprot"))
    db_path = resolve_path(db) if not Path(db).is_absolute() else Path(db)
    n_threads = resolve_threads(b.get("num_threads", 0))
    try:
        mt_mode = int(b.get("mt_mode", 1))
    except (TypeError, ValueError):
        mt_mode = 1
    engine = str(b.get("engine") or b.get("program") or "blastp")
    return run_blast(
        query_fasta,
        out_dir,
        program=str(b.get("program", "blastp")),
        database=str(db_path),
        evalue=float(b.get("evalue", 1e-5)),
        max_target_seqs=int(b.get("max_target_seqs", 100)),
        hitlist_size=int(b.get("hitlist_size", 50)),
        mode=str(b.get("mode", "remote")),
        sleep_sec=float(b.get("sleep_sec", 3)),
        force=force,
        source_tag=source_tag,
        num_threads=n_threads,
        mt_mode=mt_mode,
        engine=engine,
        diamond_sensitivity=str(b.get("diamond_sensitivity") or ""),
    )


def step_afdb(cfg: dict[str, Any], fasta: Path, out_dir: Path, force: bool) -> dict[str, Any]:
    af = cfg.get("alphafold") or {}
    return download_structures_for_fasta(
        fasta,
        out_dir,
        versions=list(af.get("versions") or ["v4", "v6", "v3", "v2"]),
        timeout=float(af.get("timeout_sec", 60)),
        force=force,
        workers=int(af.get("workers", 16)),
        uniprot_only=bool(af.get("uniprot_only", True)),
    )


def step_foldseek(
    cfg: dict[str, Any],
    structure_dir: Path,
    out_dir: Path,
    *,
    force: bool,
    source_tag: str,
) -> dict[str, Any]:
    fs = cfg.get("foldseek") or {}
    return run_foldseek_all_vs_all(
        structure_dir,
        out_dir,
        sensitivity=float(fs.get("sensitivity", 7.5)),
        evalue=float(fs.get("evalue", 0.001)),
        threads=resolve_threads(fs.get("threads", 0)),
        force=force,
        source_tag=source_tag,
    )


def step_foldseek_afdb(
    cfg: dict[str, Any],
    query_fasta: Path,
    out_dir: Path,
    structure_dir: Path,
    *,
    force: bool,
) -> dict[str, Any]:
    """用序列搜索 Foldseek-AFDB，补全非 UniProt（如 GenBank）命中的结构。"""
    fs = cfg.get("foldseek") or {}
    af = cfg.get("alphafold") or {}
    db = str(fs.get("afdb_database") or "/db/ncbi/foldseek/afdb_swissprot")
    db_path = resolve_path(db) if not Path(db).is_absolute() else Path(db)
    diamond_db = str(fs.get("afdb_diamond_db") or "/db/ncbi/foldseek/swissprot")
    diamond_path = resolve_path(diamond_db) if not Path(diamond_db).is_absolute() else Path(diamond_db)
    return run_foldseek_seq_vs_afdb(
        query_fasta,
        db_path,
        out_dir,
        structure_dir,
        sensitivity=float(fs.get("sensitivity", 7.5)),
        evalue=float(fs.get("evalue", 0.001)),
        threads=resolve_threads(fs.get("threads", 0)),
        max_seqs=int(fs.get("afdb_max_seqs", 5)),
        force=force,
        source_tag="foldseek_afdb",
        download_structures=bool(fs.get("afdb_download_structures", True)),
        af_versions=list(af.get("versions") or ["v4", "v6", "v3", "v2"]),
        af_timeout=float(af.get("timeout_sec", 45)),
        af_workers=int(af.get("workers", 16)),
        backend=str(fs.get("afdb_backend") or "auto"),
        diamond_db=diamond_path,
    )


def step_dedupe(cfg: dict[str, Any], force: bool) -> dict[str, Any]:
    paths = _cfg_paths(cfg)
    results = paths["results"]
    dedup_dir = ensure_dir(results / "dedup")
    merged = dedup_dir / "merged_round1.fasta"
    unique = dedup_dir / "unique.fasta"
    cluster_map = dedup_dir / "cluster_map.tsv"

    inputs = [
        _seed_fasta(cfg),
        results / "round1" / "blast_hits.fasta",
        results / "round1" / "foldseek_afdb_hits.fasta",
        results / "round1" / "foldseek_hits.fasta",
    ]
    mr = merge_fastas(inputs, merged)
    if not mr.get("ok") or int(mr.get("count") or 0) == 0:
        return result_fail("一轮合并后序列为空", **mr)

    dcfg = cfg.get("dedupe") or {}
    return run_cdhit(
        merged,
        unique,
        cluster_map,
        identity=float(dcfg.get("identity", 0.9)),
        word_size=int(dcfg.get("word_size", 5)),
        threads=resolve_threads(dcfg.get("threads", 0)),
        force=force,
    )


def step_finalize(cfg: dict[str, Any], force: bool) -> dict[str, Any]:
    """最终产物 = NCBI 序列 ∩ AlphaFold DB 结构，再 CD-HIT 去重。"""
    paths = _cfg_paths(cfg)
    results = paths["results"]
    final_dir = ensure_dir(results / "final")
    candidates_fasta = final_dir / "ncbi_candidates.fasta"
    out_fasta = final_dir / "complete_enzymes.fasta"
    out_tsv = final_dir / "complete_summary.tsv"
    dropped_tsv = final_dir / "dropped_no_afdb.tsv"
    struct_link = ensure_dir(final_dir / "structures")
    struct_all = ensure_dir(paths["structures"] / "all")
    seed_fa = _seed_fasta(cfg)
    seed_tsv = _seed_tsv(cfg)

    inputs = [
        seed_fa,
        results / "round1" / "blast_hits.fasta",
        results / "round1" / "foldseek_afdb_hits.fasta",
        results / "round1" / "foldseek_hits.fasta",
        results / "round2" / "blast_hits.fasta",
        results / "round2" / "foldseek_hits.fasta",
        results / "dedup" / "unique.fasta",
    ]
    mr = merge_fastas([p for p in inputs if p.exists()], candidates_fasta)
    if not mr.get("ok") or int(mr.get("count") or 0) == 0:
        return result_fail("NCBI 侧候选序列为空", **mr)
    logger.info("NCBI 侧候选（去重前合并）: %s 条", mr.get("count"))

    logger.info("=== 最终: 为全部 NCBI 候选下载 AlphaFold DB ===")
    af = step_afdb(cfg, candidates_fasta, struct_all, force=False)
    if not af.get("ok"):
        logger.warning("最终 AFDB 批量下载返回异常: %s", af.get("reason"))

    struct_ok = collect_structure_accessions(
        paths["structures"] / "seed",
        paths["structures"] / "round2",
        paths["structures"] / "merged_r2",
        struct_all,
        struct_link,
    )

    source_map: dict[str, set[str]] = {}

    def _tag_fasta(path: Path, tag: str) -> None:
        if not path.exists():
            return
        for header, _ in parse_fasta(path):
            acc = accession_from_header(header).upper()
            if acc:
                source_map.setdefault(acc, set()).add(tag)

    _tag_fasta(seed_fa, "uniprot_seed")
    _tag_fasta(results / "round1" / "blast_hits.fasta", "blast_r1")
    _tag_fasta(results / "round1" / "foldseek_afdb_hits.fasta", "foldseek_afdb")
    _tag_fasta(results / "round1" / "foldseek_hits.fasta", "foldseek_r1")
    _tag_fasta(results / "round2" / "blast_hits.fasta", "blast_r2")
    _tag_fasta(results / "round2" / "foldseek_hits.fasta", "foldseek_r2")

    meta: dict[str, dict[str, str]] = {}
    if seed_tsv.exists():
        lines = seed_tsv.read_text(encoding="utf-8", errors="ignore").splitlines()
        if lines:
            cols = lines[0].split("\t")
            for line in lines[1:]:
                parts = line.split("\t")
                row = {cols[i]: parts[i] if i < len(parts) else "" for i in range(len(cols))}
                acc = (row.get("accession") or "").upper()
                if acc:
                    meta[acc] = row

    both_records: list[tuple[str, str]] = []
    dropped_lines = ["accession\tsources\treason"]
    for header_fa, seq in parse_fasta(candidates_fasta):
        acc = accession_from_header(header_fa).upper()
        sources = ";".join(sorted(source_map.get(acc, {"unknown"})))
        if acc in struct_ok:
            both_records.append((f"{acc} {header_fa}", seq))
        else:
            dropped_lines.append(f"{acc}\t{sources}\tno_alphafold_db")

    atomic_write_text(dropped_tsv, "\n".join(dropped_lines) + "\n")
    both_fasta = final_dir / "ncbi_afdb_intersection.fasta"
    write_fasta(both_fasta, both_records)
    logger.info(
        "NCBI∩AFDB 交集: %s 条（丢弃无结构 %s 条）",
        len(both_records),
        len(dropped_lines) - 1,
    )
    if not both_records:
        return result_fail(
            "NCBI 与 AlphaFold DB 无共同结果",
            candidates=int(mr.get("count") or 0),
            dropped=str(dropped_tsv),
        )

    dcfg = cfg.get("dedupe") or {}
    cluster_map = final_dir / "cluster_map.tsv"
    dr = run_cdhit(
        both_fasta,
        out_fasta,
        cluster_map,
        identity=float(dcfg.get("identity", 0.9)),
        word_size=int(dcfg.get("word_size", 5)),
        threads=resolve_threads(dcfg.get("threads", 0)),
        force=True,
    )
    if not dr.get("ok"):
        return result_fail(dr.get("reason", "最终去重失败"), **dr)

    final_accs = {
        accession_from_header(h).upper() for h, _ in parse_fasta(out_fasta)
    }
    header = (
        "accession\tsources\tgene_names\tprotein_name\torganism_name\t"
        "length\tec\treviewed\thas_afdb\thas_ncbi"
    )
    rows = [header]
    for header_fa, seq in parse_fasta(out_fasta):
        acc = accession_from_header(header_fa).upper()
        m = meta.get(acc, {})
        sources = ";".join(sorted(source_map.get(acc, {"unknown"})))
        rows.append(
            "\t".join(
                [
                    acc,
                    sources,
                    m.get("gene_names", m.get("gene", "")),
                    m.get("protein_name", ""),
                    m.get("organism_name", m.get("organism", "")),
                    m.get("length", str(len(seq))),
                    m.get("ec", ""),
                    m.get("reviewed", ""),
                    "yes",
                    "yes",
                ]
            )
        )
    atomic_write_text(out_tsv, "\n".join(rows) + "\n")

    for old in list(struct_link.glob("AF-*-F1.pdb")) + list(struct_link.glob("AF-*-F1.cif")):
        try:
            old.unlink()
        except OSError:
            pass
    src_dirs = [
        struct_all,
        paths["structures"] / "seed",
        paths["structures"] / "round2",
        paths["structures"] / "merged_r2",
    ]
    copied = 0
    for acc in sorted(final_accs):
        dest = struct_link / f"AF-{acc}-F1.pdb"
        for sdir in src_dirs:
            src = sdir / f"AF-{acc}-F1.pdb"
            if src.exists() and src.stat().st_size > 1000:
                try:
                    shutil.copy2(src, dest)
                    copied += 1
                    break
                except OSError as e:
                    logger.warning("复制结构失败 %s: %s", src, e)

    index_lines = ["accession\tpdb_path\tstatus"]
    for acc in sorted(final_accs):
        pdb = struct_link / f"AF-{acc}-F1.pdb"
        if pdb.exists():
            index_lines.append(f"{acc}\t{pdb}\tok")
        else:
            index_lines.append(f"{acc}\t\tmissing")
    atomic_write_text(final_dir / "structure_index.tsv", "\n".join(index_lines) + "\n")

    n = len(parse_fasta(out_fasta))
    if n == 0:
        return result_fail("最终结果为空", fasta=str(out_fasta), tsv=str(out_tsv))
    logger.info(
        "最终完整结果（NCBI∩AFDB 去重）: %s 条，PDB=%s → %s",
        n,
        copied,
        out_tsv,
    )
    return result_ok(
        fasta=str(out_fasta),
        tsv=str(out_tsv),
        structures=str(struct_link),
        count=n,
        pdb_count=copied,
        ncbi_candidates=int(mr.get("count") or 0),
        intersection=len(both_records),
        dropped_no_afdb=len(dropped_lines) - 1,
        dropped_tsv=str(dropped_tsv),
    )


def step_dlkcat(cfg: dict[str, Any], force: bool) -> dict[str, Any]:
    paths = _cfg_paths(cfg)
    final_dir = paths["results"] / "final"
    enzyme_fasta = final_dir / "complete_enzymes.fasta"
    if not enzyme_fasta.exists() or not parse_fasta(enzyme_fasta):
        return result_fail(
            "缺少最终酶集 complete_enzymes.fasta，请先跑 finalize",
            fasta=str(enzyme_fasta),
        )
    # 先确保底物可用
    sr = resolve_substrates(cfg.get("substrates"), paths["substrates"], force=False)
    if not sr.get("ok"):
        return sr
    return run_dlkcat_prediction(enzyme_fasta, final_dir, cfg=cfg, force=force)


def run_pipeline(cfg: dict[str, Any], steps: list[str], force: bool) -> dict[str, Any]:
    paths = _cfg_paths(cfg)
    ensure_dir(paths["raw"])
    ensure_dir(paths["structures"])
    ensure_dir(paths["substrates"])
    ensure_dir(paths["results"] / "round1")
    ensure_dir(paths["results"] / "round2")
    ensure_dir(paths["results"] / "dedup")
    ensure_dir(paths["results"] / "final")

    summary: dict[str, Any] = {}
    seed_fasta = _seed_fasta(cfg)
    blast_query = _blast_query_fasta(cfg)

    if "uniprot" in steps:
        logger.info("=== 步骤: UniProtKB 搜索下载 (keyword=%s) ===", cfg.get("keyword"))
        r = step_uniprot(cfg, force)
        summary["uniprot"] = r
        if not r.get("ok"):
            return result_fail(r.get("reason", "UniProt 失败"), steps=summary)

    if "seed_dedupe" in steps:
        logger.info("=== 步骤: 种子 CD-HIT 去重 ===")
        r = step_seed_dedupe(cfg, force)
        summary["seed_dedupe"] = r
        if not r.get("ok"):
            return result_fail(r.get("reason", "种子去重失败"), steps=summary)
        blast_query = _blast_query_fasta(cfg)

    if "blast1" in steps:
        logger.info("=== 步骤: NCBI BLAST 一轮 (query=%s) ===", blast_query.name)
        r = step_blast(
            cfg,
            blast_query,
            paths["results"] / "round1",
            force=force,
            source_tag="blast_r1",
        )
        summary["blast1"] = r
        if not r.get("ok"):
            logger.warning("BLAST 一轮失败: %s（继续后续步骤）", r.get("reason"))

    if "afdb" in steps:
        logger.info("=== 步骤: AlphaFold DB 结构下载 ===")
        af_fasta = paths["results"] / "round1" / "afdb_input.fasta"
        merge_fastas(
            [
                blast_query,
                paths["results"] / "round1" / "blast_hits.fasta",
            ],
            af_fasta,
        )
        r = step_afdb(cfg, af_fasta, paths["structures"] / "seed", force)
        summary["afdb"] = r
        if not r.get("ok"):
            logger.warning("AFDB 下载失败: %s", r.get("reason"))

    if "foldseek_afdb" in steps:
        logger.info("=== 步骤: Foldseek 序列搜 AFDB ===")
        af_fasta = paths["results"] / "round1" / "afdb_input.fasta"
        if not af_fasta.exists() or not parse_fasta(af_fasta):
            merge_fastas(
                [
                    blast_query,
                    paths["results"] / "round1" / "blast_hits.fasta",
                ],
                af_fasta,
            )
        r = step_foldseek_afdb(
            cfg,
            af_fasta,
            paths["results"] / "round1",
            paths["structures"] / "seed",
            force=force,
        )
        summary["foldseek_afdb"] = r
        if not r.get("ok"):
            logger.warning("Foldseek-AFDB 失败: %s（继续后续步骤）", r.get("reason"))

    if "foldseek1" in steps:
        logger.info("=== 步骤: Foldseek 一轮 ===")
        r = step_foldseek(
            cfg,
            paths["structures"] / "seed",
            paths["results"] / "round1",
            force=force,
            source_tag="foldseek_r1",
        )
        summary["foldseek1"] = r
        if not r.get("ok"):
            logger.warning("Foldseek 一轮失败: %s", r.get("reason"))

    if "dedupe" in steps:
        logger.info("=== 步骤: 去重 ===")
        r = step_dedupe(cfg, force)
        summary["dedupe"] = r
        if not r.get("ok"):
            return result_fail(r.get("reason", "去重失败"), steps=summary)

    unique_fasta = paths["results"] / "dedup" / "unique.fasta"

    if "blast2" in steps:
        logger.info("=== 步骤: NCBI BLAST 二轮 ===")
        if not unique_fasta.exists() or not parse_fasta(unique_fasta):
            logger.warning("去重代表集为空，跳过 BLAST 二轮")
            summary["blast2"] = result_fail("unique.fasta 为空")
        else:
            r = step_blast(
                cfg,
                unique_fasta,
                paths["results"] / "round2",
                force=force,
                source_tag="blast_r2",
            )
            summary["blast2"] = r

    if "foldseek2" in steps:
        logger.info("=== 步骤: Foldseek 二轮（含新结构下载）===")
        if unique_fasta.exists():
            step_afdb(cfg, unique_fasta, paths["structures"] / "round2", force)
            merged_struct = paths["structures"] / "merged_r2"
            ensure_dir(merged_struct)
            for src in (paths["structures"] / "seed", paths["structures"] / "round2"):
                if not src.exists():
                    continue
                for pdb in list(src.glob("*.pdb")) + list(src.glob("*.cif")):
                    dest = merged_struct / pdb.name
                    if force or not dest.exists():
                        try:
                            shutil.copy2(pdb, dest)
                        except OSError:
                            pass
            r = step_foldseek(
                cfg,
                merged_struct,
                paths["results"] / "round2",
                force=force,
                source_tag="foldseek_r2",
            )
            summary["foldseek2"] = r
        else:
            summary["foldseek2"] = result_fail("unique.fasta 不存在")

    if "finalize" in steps:
        logger.info("=== 步骤: 汇总最终结果（NCBI∩AFDB 去重）===")
        r = step_finalize(cfg, force)
        summary["finalize"] = r
        if not r.get("ok"):
            return result_fail(r.get("reason", "汇总失败"), steps=summary)

    if "dlkcat" in steps:
        logger.info("=== 步骤: DLKcat kcat 预测 ===")
        r = step_dlkcat(cfg, force)
        summary["dlkcat"] = r
        if not r.get("ok"):
            return result_fail(r.get("reason", "DLKcat 失败"), steps=summary)

    return result_ok(steps=summary)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="关键词酶搜索流水线：UniProtKB → NCBI → AFDB → 去重 → DLKcat",
    )
    p.add_argument(
        "--config",
        default=str(ROOT / "config.yaml"),
        help="配置文件路径",
    )
    p.add_argument(
        "--keyword",
        default=None,
        help="搜索关键词（优先于 config.keyword；自动生成 UniProt 查询）",
    )
    p.add_argument(
        "--step",
        default="all",
        help="执行步骤: all 或逗号分隔 " + "|".join(STEPS),
    )
    p.add_argument("--force", action="store_true", help="强制重跑，覆盖已有结果")
    p.add_argument("-v", "--verbose", action="store_true", help="调试日志")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    setup_logging(logging.DEBUG if args.verbose else logging.INFO)

    try:
        cfg = load_config(args.config)
    except Exception as e:
        logger.error("加载配置失败: %s", e)
        return 1

    cfg = apply_keyword_to_cfg(cfg, args.keyword)
    force = bool(args.force or (cfg.get("pipeline") or {}).get("force"))
    if args.step == "all":
        steps = list(STEPS)
    else:
        steps = [s.strip() for s in args.step.split(",") if s.strip()]
        unknown = [s for s in steps if s not in STEPS]
        if unknown:
            logger.error("未知步骤: %s；可选: %s", unknown, STEPS)
            return 1

    logger.info(
        "开始流水线 keyword=%s query=%s steps=%s force=%s threads(blast/foldseek/cdhit)=%s/%s/%s",
        cfg.get("keyword"),
        cfg.get("query"),
        steps,
        force,
        resolve_threads((cfg.get("blast") or {}).get("num_threads", 0)),
        resolve_threads((cfg.get("foldseek") or {}).get("threads", 0)),
        resolve_threads((cfg.get("dedupe") or {}).get("threads", 0)),
    )
    try:
        result = run_pipeline(cfg, steps, force)
    except Exception as e:
        logger.exception("流水线异常退出")
        print({"ok": False, "reason": str(e)}, file=sys.stderr)
        return 1

    if result.get("ok"):
        logger.info("流水线完成")
        fin = (result.get("steps") or {}).get("finalize") or {}
        if fin.get("tsv"):
            logger.info("最终汇总: %s", fin["tsv"])
        dl = (result.get("steps") or {}).get("dlkcat") or {}
        if dl.get("predictions"):
            logger.info("DLKcat 预测: %s", dl["predictions"])
        return 0

    logger.error("流水线失败: %s", result.get("reason"))
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
