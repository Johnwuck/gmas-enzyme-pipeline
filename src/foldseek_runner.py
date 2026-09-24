"""Foldseek 结构比对。"""

from __future__ import annotations

import logging
import shutil
import subprocess
from pathlib import Path
from typing import Any

from .io_utils import (
    accession_from_header,
    atomic_write_text,
    ensure_dir,
    result_fail,
    result_ok,
    should_skip,
)
from .uniprot_client import fetch_sequences_by_accessions

logger = logging.getLogger(__name__)


def _which(cmd: str) -> str | None:
    return shutil.which(cmd)


def _list_pdbs(structure_dir: Path) -> list[Path]:
    return sorted(structure_dir.glob("*.pdb")) + sorted(structure_dir.glob("*.cif"))


def run_foldseek_all_vs_all(
    structure_dir: Path,
    out_dir: Path,
    *,
    sensitivity: float = 7.5,
    evalue: float = 0.001,
    threads: int = 4,
    force: bool = False,
    source_tag: str = "foldseek_r1",
) -> dict[str, Any]:
    """对目录内 AFDB 结构做 Foldseek all-vs-all easy-search。"""
    ensure_dir(out_dir)
    out_tsv = out_dir / "foldseek.tsv"
    hits_fasta = out_dir / "foldseek_hits.fasta"

    if should_skip(out_tsv, force) and should_skip(hits_fasta, force):
        logger.info("跳过 Foldseek（结果已存在）: %s", out_dir)
        return result_ok(
            foldseek_tsv=str(out_tsv),
            hits_fasta=str(hits_fasta),
            skipped=True,
            source=source_tag,
        )

    foldseek = _which("foldseek")
    pdbs = _list_pdbs(structure_dir)
    if not pdbs:
        # 无结构时写空结果，不阻断流水线
        header = (
            "query\ttarget\tidentity\talnlen\tmismatch\tgapopen\t"
            "qstart\tqend\ttstart\ttend\tevalue\tbits\tsource"
        )
        atomic_write_text(out_tsv, header + "\n")
        atomic_write_text(hits_fasta, "")
        return result_ok(
            foldseek_tsv=str(out_tsv),
            hits_fasta=str(hits_fasta),
            hit_count=0,
            reason="无可用 PDB 结构，跳过 Foldseek",
            source=source_tag,
            ok=True,
        )

    if not foldseek:
        logger.warning("未找到 foldseek，写出空表（可 conda install foldseek）")
        header = (
            "query\ttarget\tidentity\talnlen\tmismatch\tgapopen\t"
            "qstart\tqend\ttstart\ttend\tevalue\tbits\tsource"
        )
        atomic_write_text(out_tsv, header + "\n")
        # 仍可根据已有结构 accession 回拉序列作为结构侧候选
        accs = []
        for p in pdbs:
            # AF-P12345-F1.pdb
            name = p.stem
            if name.startswith("AF-"):
                accs.append(name.split("-")[1])
            else:
                accs.append(accession_from_header(name))
        fr = fetch_sequences_by_accessions(accs, hits_fasta)
        return result_ok(
            foldseek_tsv=str(out_tsv),
            hits_fasta=str(hits_fasta),
            hit_count=0,
            seq_count=fr.get("count", 0),
            reason="foldseek 未安装，仅保留已下载结构对应序列",
            source=source_tag,
        )

    db_dir = out_dir / "fs_db"
    ensure_dir(db_dir)
    db_path = db_dir / "structures"
    tmp_dir = out_dir / "fs_tmp"
    ensure_dir(tmp_dir)

    # 用目录创建 DB
    try:
        subprocess.run(
            [foldseek, "createdb", str(structure_dir), str(db_path)],
            check=True,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError as e:
        logger.exception("foldseek createdb 失败")
        return result_fail(f"foldseek createdb 失败: {e.stderr or e}")

    logger.info("Foldseek easy-search threads=%s pdbs=%s", threads, len(pdbs))
    raw_out = out_dir / "foldseek_raw.tsv"
    try:
        subprocess.run(
            [
                foldseek,
                "easy-search",
                str(db_path),
                str(db_path),
                str(raw_out),
                str(tmp_dir),
                "-s",
                str(sensitivity),
                "-e",
                str(evalue),
                "--threads",
                str(threads),
                "--format-output",
                "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits",
            ],
            check=True,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError as e:
        logger.exception("foldseek easy-search 失败")
        return result_fail(f"foldseek easy-search 失败: {e.stderr or e}")

    # 加 source 列并解析 accession
    header = (
        "query\ttarget\tidentity\talnlen\tmismatch\tgapopen\t"
        "qstart\tqend\ttstart\ttend\tevalue\tbits\tsource\ttarget_accession"
    )
    lines = [header]
    accessions: list[str] = []
    hit_count = 0
    if raw_out.exists():
        for line in raw_out.read_text(encoding="utf-8", errors="ignore").splitlines():
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 12:
                continue
            query, target = parts[0], parts[1]
            if query == target:
                continue  # 跳过自比对
            t_acc = _acc_from_structure_id(target)
            if t_acc:
                accessions.append(t_acc)
            lines.append("\t".join(parts[:12] + [source_tag, t_acc]))
            hit_count += 1
    atomic_write_text(out_tsv, "\n".join(lines) + "\n")

    fr = fetch_sequences_by_accessions(accessions, hits_fasta)
    logger.info("Foldseek 完成: pairs=%s seqs=%s", hit_count, fr.get("count"))
    return result_ok(
        foldseek_tsv=str(out_tsv),
        hits_fasta=str(hits_fasta),
        hit_count=hit_count,
        seq_count=fr.get("count", 0),
        source=source_tag,
    )


def _acc_from_structure_id(sid: str) -> str:
    # AF-P12345-F1 or file stem
    s = Path(sid).stem if "/" in sid or sid.endswith(".pdb") else sid
    if s.startswith("AF-"):
        parts = s.split("-")
        if len(parts) >= 2:
            return parts[1].upper()
    return accession_from_header(s).upper()


def run_foldseek_seq_vs_afdb(
    query_fasta: Path,
    afdb_db: Path,
    out_dir: Path,
    structure_dir: Path,
    *,
    sensitivity: float = 7.5,
    evalue: float = 0.001,
    threads: int = 8,
    max_seqs: int = 5,
    force: bool = False,
    source_tag: str = "foldseek_afdb",
    download_structures: bool = True,
    af_versions: list[str] | None = None,
    af_timeout: float = 45,
    af_workers: int = 16,
    backend: str = "auto",
    diamond_db: Path | None = None,
) -> dict[str, Any]:
    """用氨基酸序列映射到 AFDB 结构。

    backend:
      - diamond: DIAMOND 搜 Swiss-Prot（本地），再下 AFDB PDB（推荐，query 为 FASTA）
      - foldseek: 仅当 query 为结构时可用；FASTA 会自动回退 diamond
        （Foldseek easy-search 输入须为 PDB/mmCIF，不能解析氨基酸 FASTA）
      - auto: 默认 diamond（序列→AFDB）
    """
    from .alphafold_db import download_structures_for_fasta
    from .io_utils import parse_fasta, write_fasta

    ensure_dir(out_dir)
    ensure_dir(structure_dir)
    out_tsv = out_dir / "foldseek_afdb.tsv"
    hits_fasta = out_dir / "foldseek_afdb_hits.fasta"
    raw_out = out_dir / "foldseek_afdb_raw.tsv"

    if should_skip(out_tsv, force) and should_skip(hits_fasta, force):
        logger.info("跳过 Foldseek-AFDB 序列搜索（结果已存在）")
        return result_ok(
            foldseek_tsv=str(out_tsv),
            hits_fasta=str(hits_fasta),
            skipped=True,
            source=source_tag,
        )

    records = parse_fasta(query_fasta)
    if not records:
        return result_fail("Foldseek-AFDB query FASTA 为空", fasta=str(query_fasta))

    filtered: list[tuple[str, str]] = []
    for header, seq in records:
        acc = accession_from_header(header).upper()
        pdb = structure_dir / f"AF-{acc}-F1.pdb"
        if pdb.exists() and pdb.stat().st_size > 1000:
            continue
        if seq:
            filtered.append((header, seq))
    query_use = out_dir / "foldseek_afdb_query.fasta"
    if filtered:
        write_fasta(query_use, filtered)
    else:
        write_fasta(query_use, records)
        filtered = records

    db = Path(str(afdb_db))
    foldseek_ok = db.exists() or Path(str(db) + ".dbtype").exists()
    mode = (backend or "auto").strip().lower()
    # Foldseek easy-search 不能把氨基酸 FASTA 当 query（会报 No structures found）
    if mode in {"auto", "foldseek"}:
        if mode == "foldseek":
            logger.warning(
                "afdb_backend=foldseek 但 query 为 FASTA；"
                "Foldseek easy-search 仅接受 PDB/mmCIF，回退 DIAMOND→Swiss-Prot→AFDB"
            )
        mode = "diamond"
    if mode == "foldseek" and not foldseek_ok:
        logger.warning("Foldseek-AFDB 库不存在，回退 DIAMOND→Swiss-Prot→AFDB")
        mode = "diamond"

    logger.info(
        "序列搜 AFDB: backend=%s queries=%s/%s threads=%s",
        mode,
        len(filtered),
        len(records),
        threads,
    )

    accessions: list[str] = []
    hit_count = 0
    lines = [
        "query\ttarget\tidentity\talnlen\tmismatch\tgapopen\t"
        "qstart\tqend\ttstart\ttend\tevalue\tbits\tsource\ttarget_accession"
    ]

    if mode == "foldseek":
        foldseek = _which("foldseek")
        if not foldseek:
            return result_fail("未找到 foldseek")
        tmp_dir = out_dir / "fs_afdb_tmp"
        ensure_dir(tmp_dir)
        try:
            subprocess.run(
                [
                    foldseek,
                    "easy-search",
                    str(query_use),
                    str(db),
                    str(raw_out),
                    str(tmp_dir),
                    "-s",
                    str(sensitivity),
                    "-e",
                    str(evalue),
                    "--max-seqs",
                    str(max(1, int(max_seqs))),
                    "--threads",
                    str(max(1, int(threads))),
                    "--format-output",
                    "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits",
                ],
                check=True,
                capture_output=True,
                text=True,
            )
        except subprocess.CalledProcessError as e:
            logger.exception("foldseek easy-search AFDB 失败")
            return result_fail(f"foldseek AFDB 搜索失败: {e.stderr or e}")
        if raw_out.exists():
            for line in raw_out.read_text(encoding="utf-8", errors="ignore").splitlines():
                if not line.strip():
                    continue
                parts = line.split("\t")
                if len(parts) < 12:
                    continue
                t_acc = _acc_from_structure_id(parts[1])
                if t_acc:
                    accessions.append(t_acc)
                lines.append("\t".join(parts[:12] + [source_tag, t_acc]))
                hit_count += 1
    else:
        # DIAMOND → Swiss-Prot UniProt accession → AFDB PDB
        diamond = _which("diamond")
        if not diamond:
            return result_fail("未找到 diamond（序列搜 AFDB 回退需要）")
        ddb = Path(str(diamond_db or "/db/ncbi/foldseek/swissprot"))
        if not ddb.exists() and not Path(str(ddb) + ".dmnd").exists():
            # 尝试常见路径
            for cand in (
                ddb,
                Path("/db/ncbi/foldseek/swissprot"),
                Path(str(ddb) + ".dmnd"),
            ):
                if cand.exists() or Path(str(cand) + ".dmnd").exists():
                    ddb = cand.with_suffix("") if cand.suffix == ".dmnd" else cand
                    break
            else:
                return result_fail(
                    f"DIAMOND Swiss-Prot 库不存在: {ddb} "
                    "（可从本地 BLAST swissprot 导出后 diamond makedb）"
                )
        dmnd_path = ddb if str(ddb).endswith(".dmnd") else Path(str(ddb) + ".dmnd")
        if not dmnd_path.exists() and ddb.exists():
            dmnd_path = ddb
        try:
            subprocess.run(
                [
                    diamond,
                    "blastp",
                    "-q",
                    str(query_use),
                    "-d",
                    str(ddb),
                    "-o",
                    str(raw_out),
                    "--outfmt",
                    "6",
                    "qseqid",
                    "sseqid",
                    "pident",
                    "length",
                    "mismatch",
                    "gapopen",
                    "qstart",
                    "qend",
                    "sstart",
                    "send",
                    "evalue",
                    "bitscore",
                    "-e",
                    str(evalue),
                    "-k",
                    str(max(1, int(max_seqs))),
                    "--threads",
                    str(max(1, int(threads))),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
        except subprocess.CalledProcessError as e:
            logger.exception("DIAMOND Swiss-Prot 搜索失败")
            return result_fail(f"DIAMOND→Swiss-Prot 失败: {e.stderr or e}")
        if raw_out.exists():
            for line in raw_out.read_text(encoding="utf-8", errors="ignore").splitlines():
                if not line.strip():
                    continue
                parts = line.split("\t")
                if len(parts) < 12:
                    continue
                t_acc = accession_from_header(parts[1]).upper()
                if t_acc:
                    accessions.append(t_acc)
                lines.append("\t".join(parts[:12] + [f"{source_tag}_diamond", t_acc]))
                hit_count += 1

    atomic_write_text(out_tsv, "\n".join(lines) + "\n")

    uniq_accs = sorted({a.upper() for a in accessions if a})
    stub = out_dir / "foldseek_afdb_targets.fasta"
    write_fasta(stub, [(a, "M") for a in uniq_accs])

    downloaded = 0
    if download_structures and uniq_accs:
        fr = fetch_sequences_by_accessions(uniq_accs, hits_fasta)
        dl_fasta = hits_fasta if fr.get("ok") and int(fr.get("count") or 0) > 0 else stub
        dr = download_structures_for_fasta(
            dl_fasta,
            structure_dir,
            versions=af_versions or ["v4", "v6", "v3", "v2"],
            timeout=af_timeout,
            force=force,
            workers=af_workers,
            uniprot_only=True,
        )
        downloaded = int(dr.get("ok_count") or 0)
    else:
        write_fasta(hits_fasta, [])

    logger.info(
        "序列搜 AFDB 完成: backend=%s pairs=%s unique_targets=%s structures=%s",
        mode,
        hit_count,
        len(uniq_accs),
        downloaded,
    )
    return result_ok(
        foldseek_tsv=str(out_tsv),
        hits_fasta=str(hits_fasta),
        hit_count=hit_count,
        target_count=len(uniq_accs),
        structure_count=downloaded,
        source=source_tag,
        backend=mode,
    )
