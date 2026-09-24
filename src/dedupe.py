"""按 accession 与 CD-HIT 去重。"""

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
    parse_fasta,
    result_fail,
    result_ok,
    should_skip,
    write_fasta,
)

logger = logging.getLogger(__name__)


def merge_fastas(paths: list[Path], out_fasta: Path) -> dict[str, Any]:
    """合并多个 FASTA，按 accession 首次出现保留。"""
    seen: set[str] = set()
    records: list[tuple[str, str]] = []
    for p in paths:
        if not p or not Path(p).exists():
            continue
        for header, seq in parse_fasta(Path(p)):
            acc = accession_from_header(header).upper()
            if not acc or acc in seen or not seq:
                continue
            seen.add(acc)
            records.append((f"{acc} {header}", seq))
    write_fasta(out_fasta, records)
    return result_ok(path=str(out_fasta), count=len(records))


def dedupe_by_accession(fasta_path: Path, out_fasta: Path) -> dict[str, Any]:
    return merge_fastas([fasta_path], out_fasta)


def run_cdhit(
    fasta_path: Path,
    out_fasta: Path,
    cluster_map: Path,
    *,
    identity: float = 0.9,
    word_size: int = 5,
    threads: int = 1,
    force: bool = False,
) -> dict[str, Any]:
    """CD-HIT 聚类；若无 cd-hit 则退化为 accession 去重。"""
    ensure_dir(out_fasta.parent)
    if should_skip(out_fasta, force) and should_skip(cluster_map, force):
        logger.info("跳过去重（结果已存在）")
        return result_ok(
            unique_fasta=str(out_fasta),
            cluster_map=str(cluster_map),
            skipped=True,
        )

    records = parse_fasta(fasta_path)
    if not records:
        write_fasta(out_fasta, [])
        atomic_write_text(cluster_map, "representative\tmember\n")
        return result_fail("待去重 FASTA 为空", unique_fasta=str(out_fasta))

    # 先 accession 去重
    tmp_acc = out_fasta.parent / "merged_by_acc.fasta"
    dedupe_by_accession(fasta_path, tmp_acc)

    cdhit = shutil.which("cd-hit") or shutil.which("cdhit")
    if not cdhit:
        logger.warning("未找到 cd-hit，仅做 accession 去重")
        shutil.copyfile(tmp_acc, out_fasta)
        lines = ["representative\tmember"]
        for header, _ in parse_fasta(out_fasta):
            acc = accession_from_header(header).upper()
            lines.append(f"{acc}\t{acc}")
        atomic_write_text(cluster_map, "\n".join(lines) + "\n")
        return result_ok(
            unique_fasta=str(out_fasta),
            cluster_map=str(cluster_map),
            count=len(parse_fasta(out_fasta)),
            method="accession_only",
        )

    logger.info("CD-HIT threads=%s identity=%s", max(1, int(threads)), identity)
    # CD-HIT 输出；-M 0 取消默认内存上限（多线程时缓冲按线程倍增，易触发 “not enough memory”）
    try:
        subprocess.run(
            [
                cdhit,
                "-i",
                str(tmp_acc),
                "-o",
                str(out_fasta),
                "-c",
                str(identity),
                "-n",
                str(word_size),
                "-T",
                str(max(1, int(threads))),
                "-M",
                "0",
                "-d",
                "0",
            ],
            check=True,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError as e:
        logger.exception("cd-hit 失败，回退 accession 去重")
        shutil.copyfile(tmp_acc, out_fasta)
        lines = ["representative\tmember"]
        for header, _ in parse_fasta(out_fasta):
            acc = accession_from_header(header).upper()
            lines.append(f"{acc}\t{acc}")
        atomic_write_text(cluster_map, "\n".join(lines) + "\n")
        return result_ok(
            unique_fasta=str(out_fasta),
            cluster_map=str(cluster_map),
            count=len(parse_fasta(out_fasta)),
            method="accession_fallback",
            reason=str(e.stderr or e),
        )

    # 解析 .clstr
    clstr = Path(str(out_fasta) + ".clstr")
    lines = ["representative\tmember"]
    if clstr.exists():
        rep = ""
        members: list[str] = []
        for line in clstr.read_text(encoding="utf-8", errors="ignore").splitlines():
            if line.startswith(">"):
                if rep:
                    for m in members:
                        lines.append(f"{rep}\t{m}")
                rep = ""
                members = []
                continue
            # 0	123aa, >P12345 ... *
            if ">" in line:
                mid = line.split(">")[1]
                acc = mid.split("...")[0].split()[0]
                acc = accession_from_header(acc).upper()
                members.append(acc)
                if line.strip().endswith("*"):
                    rep = acc
        if rep:
            for m in members:
                lines.append(f"{rep}\t{m}")
    else:
        for header, _ in parse_fasta(out_fasta):
            acc = accession_from_header(header).upper()
            lines.append(f"{acc}\t{acc}")

    atomic_write_text(cluster_map, "\n".join(lines) + "\n")
    n = len(parse_fasta(out_fasta))
    logger.info("去重完成: 代表序列 %s (identity=%s)", n, identity)
    return result_ok(
        unique_fasta=str(out_fasta),
        cluster_map=str(cluster_map),
        count=n,
        method="cd-hit",
    )
