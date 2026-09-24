"""AlphaFold DB 结构下载。"""

from __future__ import annotations

import logging
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import requests

from .io_utils import (
    accession_from_header,
    atomic_write_bytes,
    atomic_write_text,
    ensure_dir,
    looks_like_uniprot_acc,
    parse_fasta,
    result_fail,
    result_ok,
    should_skip,
)

logger = logging.getLogger(__name__)

AFDB_BASE = "https://alphafold.ebi.ac.uk/files"


def _fetch_afdb_pdb(url: str, timeout: float) -> bytes | None:
    """拉取 AFDB PDB；404 立即返回 None（不重试），网络错误有限重试。"""
    for attempt in range(3):
        try:
            r = requests.get(url, timeout=timeout)
            if r.status_code == 404:
                return None
            if r.status_code == 429 or r.status_code >= 500:
                time.sleep(min(2**attempt, 30))
                continue
            if r.status_code == 200 and len(r.content) > 1000:
                return r.content
            return None
        except requests.RequestException as e:
            logger.debug("AFDB 请求异常 %s: %s", url, e)
            time.sleep(min(2**attempt, 30))
    return None


def download_alphafold_pdb(
    uniprot_id: str,
    out_pdb: Path,
    *,
    versions: list[str] | None = None,
    timeout: float = 90,
    force: bool = False,
) -> dict[str, Any]:
    """下载 AFDB 模型 PDB。"""
    acc = uniprot_id.strip().upper().split("-")[0]
    if should_skip(out_pdb, force):
        return result_ok(path=str(out_pdb), accession=acc, skipped=True)

    versions = versions or ["v4", "v6", "v3", "v2"]
    for ver in versions:
        url = f"{AFDB_BASE}/AF-{acc}-F1-model_{ver}.pdb"
        content = _fetch_afdb_pdb(url, timeout=timeout)
        if content:
            atomic_write_bytes(out_pdb, content)
            return result_ok(path=str(out_pdb), accession=acc, version=ver)
    return result_fail(f"未找到 AlphaFold 结构: {acc}", accession=acc, path=str(out_pdb))


def _try_reuse_local(acc: str, out_pdb: Path, out_dir: Path) -> bool:
    if out_pdb.exists() and out_pdb.stat().st_size > 1000:
        return True
    if not out_dir.parent.exists():
        return False
    for other in out_dir.parent.glob(f"*/AF-{acc}-F1.pdb"):
        if other.resolve() == out_pdb.resolve():
            continue
        if other.exists() and other.stat().st_size > 1000:
            try:
                out_pdb.write_bytes(other.read_bytes())
                return True
            except OSError:
                pass
    return False


def download_structures_for_fasta(
    fasta_path: Path,
    out_dir: Path,
    *,
    versions: list[str] | None = None,
    timeout: float = 90,
    force: bool = False,
    workers: int = 16,
    uniprot_only: bool = True,
) -> dict[str, Any]:
    """为 FASTA 中每条序列尝试下载 AFDB 结构（可并行；默认仅 UniProt accession）。"""
    ensure_dir(out_dir)
    records = parse_fasta(fasta_path)
    if not records:
        return result_fail("FASTA 为空，无法下载结构", fasta=str(fasta_path))

    versions = versions or ["v4", "v6", "v3", "v2"]
    n_workers = max(1, int(workers))

    candidates: list[str] = []
    seen: set[str] = set()
    skipped_non_uniprot = 0
    for header, _seq in records:
        acc = accession_from_header(header).upper()
        if not acc or acc in seen:
            continue
        seen.add(acc)
        if uniprot_only and not looks_like_uniprot_acc(acc):
            skipped_non_uniprot += 1
            continue
        candidates.append(acc)

    logger.info(
        "AFDB 待下载: unique=%s uniprot_like=%s skipped_non_uniprot=%s workers=%s",
        len(seen),
        len(candidates),
        skipped_non_uniprot,
        n_workers,
    )

    ok_list: list[str] = []
    fail_list: list[str] = []
    index_lines = ["accession\tpdb_path\tstatus"]
    # 已存在的直接记成功
    todo: list[str] = []
    for acc in candidates:
        out_pdb = out_dir / f"AF-{acc}-F1.pdb"
        if not force and _try_reuse_local(acc, out_pdb, out_dir):
            ok_list.append(acc)
            index_lines.append(f"{acc}\t{out_pdb}\tok")
        else:
            todo.append(acc)

    logger.info("AFDB 已有/复用 %s，需请求 %s", len(ok_list), len(todo))

    def _one(acc: str) -> tuple[str, bool, str]:
        out_pdb = out_dir / f"AF-{acc}-F1.pdb"
        r = download_alphafold_pdb(
            acc, out_pdb, versions=versions, timeout=timeout, force=force
        )
        if r.get("ok"):
            return acc, True, str(out_pdb)
        return acc, False, ""

    done = 0
    if todo:
        with ThreadPoolExecutor(max_workers=n_workers) as pool:
            futs = {pool.submit(_one, acc): acc for acc in todo}
            for fut in as_completed(futs):
                acc, ok, path = fut.result()
                done += 1
                if ok:
                    ok_list.append(acc)
                    index_lines.append(f"{acc}\t{path}\tok")
                else:
                    fail_list.append(acc)
                    index_lines.append(f"{acc}\t\tfail")
                if done % 200 == 0 or done == len(todo):
                    logger.info(
                        "AFDB 下载进度: %s/%s ok=%s fail=%s",
                        done,
                        len(todo),
                        len(ok_list),
                        len(fail_list),
                    )

    index_path = out_dir / "structure_index.tsv"
    atomic_write_text(index_path, "\n".join(index_lines) + "\n")
    logger.info(
        "AFDB 下载完成: 成功 %s / 失败 %s / 跳过非UniProt %s",
        len(ok_list),
        len(fail_list),
        skipped_non_uniprot,
    )
    return result_ok(
        index=str(index_path),
        ok_count=len(ok_list),
        fail_count=len(fail_list),
        skipped_non_uniprot=skipped_non_uniprot,
        ok_accessions=ok_list,
        fail_accessions=fail_list,
        structure_dir=str(out_dir),
    )


def collect_structure_accessions(*dirs: Path) -> set[str]:
    """扫描目录中已有 AFDB PDB 的 accession。"""
    out: set[str] = set()
    for sdir in dirs:
        if not sdir or not Path(sdir).exists():
            continue
        for pdb in list(Path(sdir).glob("AF-*-F1.pdb")) + list(Path(sdir).glob("AF-*-F1.cif")):
            parts = pdb.stem.split("-")
            if len(parts) >= 2:
                out.add(parts[1].upper())
    return out
