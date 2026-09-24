"""NCBI BLAST 比对（remote NCBIWWW 或本地 blastp）。"""

from __future__ import annotations

import logging
import shutil
import subprocess
import time
from io import StringIO
from pathlib import Path
from typing import Any

from Bio.Blast import NCBIXML

from .io_utils import (
    accession_from_header,
    atomic_write_text,
    ensure_dir,
    looks_like_uniprot_acc,
    parse_fasta,
    result_fail,
    result_ok,
    should_skip,
    write_fasta,
)
from .uniprot_client import fetch_sequences_by_accessions

logger = logging.getLogger(__name__)


def _which(cmd: str) -> str | None:
    return shutil.which(cmd)


def _parse_blast_xml(xml_text: str) -> list[dict[str, Any]]:
    hits: list[dict[str, Any]] = []
    if not xml_text.strip():
        return hits
    for record in NCBIXML.parse(StringIO(xml_text)):
        qid = record.query or ""
        for aln in record.alignments:
            hit_id = aln.hit_id or aln.title or ""
            acc = accession_from_header(hit_id)
            # SwissProt/nr titles often: sp|P12345|NAME or gi|...|sp|P12345|NAME
            if "|" in (aln.title or ""):
                acc = accession_from_header(aln.title or hit_id)
            best = aln.hsps[0] if aln.hsps else None
            if best is None:
                continue
            identity = 100.0 * best.identities / best.align_length if best.align_length else 0.0
            hits.append(
                {
                    "query": qid,
                    "hit_id": hit_id,
                    "hit_def": aln.hit_def or aln.title or "",
                    "accession": acc,
                    "evalue": float(best.expect),
                    "bitscore": float(best.bits),
                    "identity": round(identity, 2),
                    "align_length": int(best.align_length),
                }
            )
    return hits


def blast_remote_one(
    sequence: str,
    *,
    program: str = "blastp",
    database: str = "swissprot",
    evalue: float = 1e-5,
    hitlist_size: int = 50,
    expect_wait: int = 10,
) -> dict[str, Any]:
    from Bio.Blast import NCBIWWW

    try:
        handle = NCBIWWW.qblast(
            program,
            database,
            sequence,
            expect=evalue,
            hitlist_size=hitlist_size,
            format_type="XML",
        )
        xml_text = handle.read()
        handle.close()
        # NCBIWWW 有时返回空或 HTML 错误页
        if not xml_text or "<HTML>" in xml_text[:200].upper():
            return result_fail("NCBI 返回非 XML 结果", xml_len=len(xml_text or ""))
        hits = _parse_blast_xml(xml_text)
        return result_ok(xml=xml_text, hits=hits, count=len(hits))
    except Exception as e:
        logger.exception("远程 BLAST 失败")
        return result_fail(f"远程 BLAST 失败: {e}")


def blast_local_fasta(
    query_fasta: Path,
    out_tsv: Path,
    *,
    database: str,
    evalue: float = 1e-5,
    max_target_seqs: int = 100,
    num_threads: int = 1,
    mt_mode: int = 0,
) -> dict[str, Any]:
    blastp = _which("blastp")
    if not blastp:
        return result_fail("未找到本地 blastp，请 conda install blast 或改用 remote 模式")
    db = Path(str(database))
    if not (
        db.exists()
        or Path(str(database) + ".pin").exists()
        or Path(str(database) + ".pal").exists()  # NR 等 alias
        or Path(str(database) + ".nal").exists()
    ):
        logger.warning("本地库路径可能不存在: %s", database)

    ensure_dir(out_tsv.parent)
    threads = max(1, int(num_threads))
    mode = int(mt_mode)
    cmd = [
        blastp,
        "-query",
        str(query_fasta),
        "-db",
        str(database),
        "-evalue",
        str(evalue),
        "-max_target_seqs",
        str(max_target_seqs),
        "-num_threads",
        str(threads),
        "-mt_mode",
        str(mode),
        "-outfmt",
        "6 qseqid sseqid sacc pident length evalue bitscore stitle",
        "-out",
        str(out_tsv),
    ]
    logger.info("本地 blastp threads=%s mt_mode=%s db=%s", threads, mode, database)
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        logger.exception("本地 blastp 失败")
        return result_fail(f"本地 blastp 失败: {e.stderr or e}")
    return result_ok(path=str(out_tsv))


def _diamond_db_path(database: str) -> Path | None:
    """解析 DIAMOND 库前缀；接受 path、path.dmnd。"""
    p = Path(str(database))
    if p.suffix == ".dmnd" and p.exists():
        return p.with_suffix("")
    if Path(str(p) + ".dmnd").exists():
        return p
    if p.exists() and p.is_file():
        return p.with_suffix("") if p.suffix == ".dmnd" else p
    return None


def diamond_local_fasta(
    query_fasta: Path,
    out_tsv: Path,
    *,
    database: str,
    evalue: float = 1e-5,
    max_target_seqs: int = 100,
    num_threads: int = 1,
    sensitivity: str = "",
    force: bool = False,
) -> dict[str, Any]:
    """DIAMOND blastp；写出与 blastp 同列的 tabular（含 sacc）。"""
    diamond = _which("diamond")
    if not diamond:
        return result_fail("未找到 diamond，请 conda install -c bioconda diamond")

    db_prefix = _diamond_db_path(database)
    if db_prefix is None:
        return result_fail(
            f"DIAMOND 库不存在: {database}.dmnd（请运行 scripts/build_diamond_nr.sh）"
        )

    ensure_dir(out_tsv.parent)
    threads = max(1, int(num_threads))
    raw_tsv = out_tsv.parent / (out_tsv.stem + ".diamond_raw.tsv")

    if not force and out_tsv.exists() and out_tsv.stat().st_size > 0:
        logger.info("复用已有规范化 DIAMOND 结果: %s", out_tsv)
        return result_ok(path=str(out_tsv), engine="diamond", skipped=True)

    reuse_raw = not force and raw_tsv.exists() and raw_tsv.stat().st_size > 0
    if reuse_raw:
        logger.info("复用已有 DIAMOND 原始结果: %s (%s bytes)", raw_tsv, raw_tsv.stat().st_size)
    else:
        cmd = [
            diamond,
            "blastp",
            "-q",
            str(query_fasta),
            "-d",
            str(db_prefix),
            "-o",
            str(raw_tsv),
            "-e",
            str(evalue),
            "-k",
            str(max_target_seqs),
            "--threads",
            str(threads),
            "--outfmt",
            "6",
            "qseqid",
            "sseqid",
            "pident",
            "length",
            "evalue",
            "bitscore",
            "stitle",
        ]
        sens = (sensitivity or "").strip()
        if sens:
            # e.g. --sensitive / --very-sensitive / --ultra-sensitive
            if not sens.startswith("-"):
                sens = f"--{sens.lstrip('-')}"
            cmd.append(sens)

        logger.info("本地 diamond blastp threads=%s db=%s", threads, db_prefix)
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            logger.exception("本地 diamond 失败")
            return result_fail(f"本地 diamond 失败: {e.stderr or e}")

    # 规范为 blast 列: qseqid sseqid sacc pident length evalue bitscore stitle
    lines_out: list[str] = []
    if raw_tsv.exists():
        with raw_tsv.open(encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 6:
                    continue
                q, sseqid, pident, length, ev, bits = parts[:6]
                title = parts[6] if len(parts) > 6 else ""
                try:
                    acc = accession_from_header(sseqid) or accession_from_header(title)
                except Exception:
                    logger.warning("accession 解析失败 sseqid=%r title=%r", sseqid[:80], title[:80])
                    acc = (sseqid or "").split()[0].split(".")[0] if sseqid else ""
                lines_out.append(
                    "\t".join([q, sseqid, acc, pident, length, ev, bits, title])
                )
    atomic_write_text(out_tsv, "\n".join(lines_out) + ("\n" if lines_out else ""))
    # 保留 raw 便于排错；成功写出正规表后可删（可选）
    logger.info("DIAMOND 规范化完成: hits=%s → %s", len(lines_out), out_tsv)
    return result_ok(path=str(out_tsv), engine="diamond", hit_rows=len(lines_out))


def _hits_from_local_tsv(tsv_path: Path) -> list[dict[str, Any]]:
    hits: list[dict[str, Any]] = []
    if not tsv_path.exists():
        return hits
    for line in tsv_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 7:
            continue
        q, sseqid, sacc, pident, length, evalue, bits = parts[:7]
        title = parts[7] if len(parts) > 7 else ""
        acc = sacc or accession_from_header(sseqid)
        hits.append(
            {
                "query": q,
                "hit_id": sseqid,
                "hit_def": title,
                "accession": acc,
                "evalue": float(evalue),
                "bitscore": float(bits),
                "identity": float(pident),
                "align_length": int(float(length)),
            }
        )
    return hits


def fetch_sequences_from_blastdb(
    entries: list[str],
    out_fasta: Path,
    *,
    database: str,
    batch_size: int = 2000,
) -> dict[str, Any]:
    """用本地 blastdbcmd 按 entry（含版本号）批量取 FASTA。"""
    blastdbcmd = _which("blastdbcmd")
    if not blastdbcmd:
        return result_fail("未找到 blastdbcmd")
    db = str(database)
    if not (
        Path(db).exists()
        or Path(db + ".pin").exists()
        or Path(db + ".pal").exists()
        or Path(db + ".nal").exists()
    ):
        return result_fail(f"本地 BLAST 库不可用: {db}")

    uniq: list[str] = []
    seen: set[str] = set()
    for e in entries:
        e = (e or "").strip()
        if not e:
            continue
        key = e.upper()
        if key in seen:
            continue
        seen.add(key)
        uniq.append(e)
    if not uniq:
        write_fasta(out_fasta, [])
        return result_ok(path=str(out_fasta), count=0, method="blastdbcmd")

    ensure_dir(out_fasta.parent)
    logger.info(
        "blastdbcmd 开始回拉: entries=%s batch_size=%s db=%s",
        len(uniq),
        batch_size,
        db,
    )
    # 追加写入单一 FASTA，避免上万临时文件
    if out_fasta.exists():
        out_fasta.unlink()
    got = 0
    n_batches = (len(uniq) + batch_size - 1) // batch_size
    for i in range(0, len(uniq), batch_size):
        batch = uniq[i : i + batch_size]
        batch_idx = i // batch_size
        id_file = out_fasta.parent / f".blastdbcmd_ids_{batch_idx}.txt"
        part_fa = out_fasta.parent / f".blastdbcmd_part_{batch_idx}.fasta"
        atomic_write_text(id_file, "\n".join(batch) + "\n")
        cmd = [
            blastdbcmd,
            "-db",
            db,
            "-entry_batch",
            str(id_file),
            "-out",
            str(part_fa),
            "-outfmt",
            "%f",
        ]
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            logger.warning(
                "blastdbcmd 批次 %s/%s 警告: %s",
                batch_idx + 1,
                n_batches,
                (e.stderr or str(e))[:300],
            )
        if part_fa.exists() and part_fa.stat().st_size > 0:
            # append
            with out_fasta.open("a", encoding="utf-8") as out_f, part_fa.open(
                encoding="utf-8", errors="ignore"
            ) as in_f:
                out_f.write(in_f.read())
            got += len(parse_fasta(part_fa))
        try:
            id_file.unlink(missing_ok=True)
            part_fa.unlink(missing_ok=True)
        except OSError:
            pass
        if (batch_idx + 1) % 20 == 0 or batch_idx + 1 == n_batches:
            logger.info(
                "blastdbcmd 进度: %s/%s batches, ~%s seqs",
                batch_idx + 1,
                n_batches,
                got,
            )

    # 按 accession 去重
    records: list[tuple[str, str]] = []
    seen_acc: set[str] = set()
    if out_fasta.exists():
        for header, seq in parse_fasta(out_fasta):
            acc = accession_from_header(header).upper() or header.split()[0].upper()
            if not acc or acc in seen_acc or not seq:
                continue
            seen_acc.add(acc)
            records.append((header, seq))
    write_fasta(out_fasta, records)
    logger.info(
        "blastdbcmd 回拉完成: requested=%s got=%s db=%s",
        len(uniq),
        len(records),
        db,
    )
    return result_ok(
        path=str(out_fasta),
        count=len(records),
        requested=len(uniq),
        method="blastdbcmd",
    )

def _looks_like_uniprot_acc(acc: str) -> bool:
    return looks_like_uniprot_acc(acc)


def run_blast(
    query_fasta: Path,
    out_dir: Path,
    *,
    program: str = "blastp",
    database: str = "swissprot",
    evalue: float = 1e-5,
    max_target_seqs: int = 100,
    hitlist_size: int = 50,
    mode: str = "remote",
    sleep_sec: float = 3.0,
    force: bool = False,
    source_tag: str = "blast_r1",
    num_threads: int = 1,
    mt_mode: int = 0,
    engine: str = "blastp",
    diamond_sensitivity: str = "",
) -> dict[str, Any]:
    """对 query FASTA 跑 BLAST/DIAMOND，写出 hits TSV，并回拉命中序列 FASTA。"""
    ensure_dir(out_dir)
    hits_tsv = out_dir / "blast_hits.tsv"
    hits_fasta = out_dir / "blast_hits.fasta"
    xml_dir = out_dir / "blast_xml"
    ensure_dir(xml_dir)

    if should_skip(hits_tsv, force) and should_skip(hits_fasta, force):
        logger.info("跳过 BLAST（结果已存在）: %s", out_dir)
        return result_ok(
            hits_tsv=str(hits_tsv),
            hits_fasta=str(hits_fasta),
            skipped=True,
            source=source_tag,
        )

    records = parse_fasta(query_fasta)
    if not records:
        return result_fail("BLAST query FASTA 为空", fasta=str(query_fasta))

    all_hits: list[dict[str, Any]] = []
    engine_name = (engine or program or "blastp").strip().lower()
    if engine_name in {"diamond", "diamond_blastp"}:
        engine_name = "diamond"
    else:
        engine_name = "blastp"

    if mode == "local":
        local_raw = out_dir / "blast_local.tsv"
        if engine_name == "diamond":
            lr = diamond_local_fasta(
                query_fasta,
                local_raw,
                database=database,
                evalue=evalue,
                max_target_seqs=max_target_seqs,
                num_threads=num_threads,
                sensitivity=diamond_sensitivity,
                force=force,
            )
            if not lr.get("ok"):
                # 大批量 query 不可回退 remote
                return result_fail(lr.get("reason", "DIAMOND 失败"), engine="diamond")
            all_hits = _hits_from_local_tsv(local_raw)
        else:
            lr = blast_local_fasta(
                query_fasta,
                local_raw,
                database=database,
                evalue=evalue,
                max_target_seqs=max_target_seqs,
                num_threads=num_threads,
                mt_mode=mt_mode,
            )
            if not lr.get("ok"):
                logger.warning("本地 BLAST 失败，回退 remote: %s", lr.get("reason"))
                mode = "remote"
            else:
                all_hits = _hits_from_local_tsv(local_raw)

    if mode == "remote":
        if engine_name == "diamond":
            return result_fail("DIAMOND 不支持 remote 模式，请设置 blast.mode: local")
        for idx, (header, seq) in enumerate(records):
            if not seq:
                continue
            logger.info("远程 BLAST %s/%s: %s (len=%s)", idx + 1, len(records), header[:60], len(seq))
            r = blast_remote_one(
                seq,
                program=program,
                database=database,
                evalue=evalue,
                hitlist_size=hitlist_size,
            )
            if not r.get("ok"):
                logger.warning("BLAST 失败 %s: %s", header[:40], r.get("reason"))
                time.sleep(sleep_sec)
                continue
            xml_path = xml_dir / f"q{idx:04d}.xml"
            atomic_write_text(xml_path, r.get("xml") or "")
            for h in r.get("hits") or []:
                h = dict(h)
                h["query"] = header
                all_hits.append(h)
            time.sleep(sleep_sec)

    # 过滤 evalue
    filtered = [h for h in all_hits if float(h.get("evalue", 1)) <= evalue]
    # 写 TSV
    header_line = "query\thit_id\taccession\tevalue\tbitscore\tidentity\talign_length\thit_def\tsource"
    lines = [header_line]
    accessions: list[str] = []
    for h in filtered:
        acc = str(h.get("accession") or "").strip()
        if acc:
            accessions.append(acc)
        lines.append(
            "\t".join(
                [
                    str(h.get("query", "")),
                    str(h.get("hit_id", "")),
                    acc,
                    str(h.get("evalue", "")),
                    str(h.get("bitscore", "")),
                    str(h.get("identity", "")),
                    str(h.get("align_length", "")),
                    str(h.get("hit_def", "")).replace("\t", " "),
                    source_tag,
                ]
            )
        )
    atomic_write_text(hits_tsv, "\n".join(lines) + "\n")

    # 回拉序列：本地 BLAST/NR 用 blastdbcmd；DIAMOND 库本身无序列，改读平行 NR
    hit_ids = [str(h.get("hit_id") or "").strip() for h in filtered if h.get("hit_id")]
    fr: dict[str, Any] = result_fail("未回拉")
    blast_seq_db = str(database)
    if "diamond" in blast_seq_db.lower() or Path(blast_seq_db + ".dmnd").exists():
        for cand in ("/db/ncbi/NR/nr", "/db/ncbi/blastdb/swissprot"):
            if Path(cand + ".pal").exists() or Path(cand + ".pin").exists():
                blast_seq_db = cand
                break

    can_blastdbcmd = bool(hit_ids) and (
        Path(blast_seq_db + ".pal").exists() or Path(blast_seq_db + ".pin").exists()
    )
    if can_blastdbcmd:
        fr = fetch_sequences_from_blastdb(hit_ids, hits_fasta, database=blast_seq_db)
        if not fr.get("ok") or int(fr.get("count") or 0) == 0:
            logger.warning(
                "blastdbcmd 回拉不足，回退 UniProt（仅 UniProt 风格 accession）: %s",
                fr.get("reason"),
            )
            up_accs = [a for a in accessions if _looks_like_uniprot_acc(a)]
            fr = fetch_sequences_by_accessions(up_accs, hits_fasta)
    else:
        up_accs = [a for a in accessions if _looks_like_uniprot_acc(a)]
        if len(up_accs) < len(set(a for a in accessions if a)):
            logger.info(
                "过滤非 UniProt accession: %s → %s",
                len(set(accessions)),
                len(set(up_accs)),
            )
        fr = fetch_sequences_by_accessions(up_accs or accessions, hits_fasta)

    if not fr.get("ok"):
        logger.warning("BLAST 命中序列回拉失败: %s", fr.get("reason"))
        atomic_write_text(hits_fasta, "")

    logger.info(
        "BLAST 完成: engine=%s hits=%s seqs=%s fetch=%s",
        engine_name,
        len(filtered),
        fr.get("count"),
        fr.get("method") or "uniprot",
    )
    return result_ok(
        hits_tsv=str(hits_tsv),
        hits_fasta=str(hits_fasta),
        hit_count=len(filtered),
        seq_count=fr.get("count", 0),
        source=source_tag,
        engine=engine_name,
    )
