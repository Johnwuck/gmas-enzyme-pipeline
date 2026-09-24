"""路径、配置、重试与原子写入工具。"""

from __future__ import annotations

import logging
import os
import tempfile
import time
from pathlib import Path
from typing import Any

import requests
import yaml

logger = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[1]


def setup_logging(level: int = logging.INFO) -> None:
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def load_config(path: str | Path | None = None) -> dict[str, Any]:
    cfg_path = Path(path) if path else ROOT / "config.yaml"
    with cfg_path.open(encoding="utf-8") as f:
        cfg = yaml.safe_load(f) or {}
    if not isinstance(cfg, dict):
        raise ValueError(f"配置文件格式错误: {cfg_path}")
    return cfg


def resolve_path(rel: str | Path, base: Path | None = None) -> Path:
    p = Path(rel)
    if p.is_absolute():
        return p
    return (base or ROOT) / p


def resolve_threads(configured: Any = 0, *, cap: int = 48, reserve: int = 8) -> int:
    """配置 >0 用配置值；0/缺省则按 CPU 与 1 分钟负载自动选取。

    超线程对 BLAST 收益有限，且机器上常有其他高负载，默认上限 48。
    """
    try:
        n_cfg = int(configured) if configured is not None else 0
    except (TypeError, ValueError):
        n_cfg = 0
    if n_cfg > 0:
        return n_cfg
    ncpu = os.cpu_count() or 4
    try:
        load1 = float(os.getloadavg()[0])
    except (OSError, AttributeError):
        load1 = 0.0
    busy = max(0, int(load1) - 1)
    auto = ncpu - max(reserve, busy)
    floor = 8 if ncpu >= 16 else 1
    return max(floor, min(int(cap), auto))


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def atomic_write_text(path: Path, text: str, encoding: str = "utf-8") -> None:
    ensure_dir(path.parent)
    fd, tmp_name = tempfile.mkstemp(dir=str(path.parent), prefix=f".{path.name}.", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding=encoding) as f:
            f.write(text)
        os.replace(tmp_name, path)
    except Exception:
        try:
            os.unlink(tmp_name)
        except OSError:
            pass
        raise


def atomic_write_bytes(path: Path, data: bytes) -> None:
    ensure_dir(path.parent)
    fd, tmp_name = tempfile.mkstemp(dir=str(path.parent), prefix=f".{path.name}.", suffix=".tmp")
    try:
        with os.fdopen(fd, "wb") as f:
            f.write(data)
        os.replace(tmp_name, path)
    except Exception:
        try:
            os.unlink(tmp_name)
        except OSError:
            pass
        raise


def request_get(
    url: str,
    *,
    params: dict[str, Any] | None = None,
    timeout: float = 60,
    max_retries: int = 5,
    headers: dict[str, str] | None = None,
) -> requests.Response:
    """带指数退避的 GET，处理 429/5xx。"""
    last_exc: Exception | None = None
    for attempt in range(max_retries):
        try:
            r = requests.get(url, params=params, timeout=timeout, headers=headers)
            if r.status_code == 429 or r.status_code >= 500:
                wait = min(2**attempt, 60)
                retry_after = r.headers.get("Retry-After")
                if retry_after and retry_after.isdigit():
                    wait = max(wait, int(retry_after))
                logger.warning(
                    "请求 %s 返回 %s，%ss 后重试 (%s/%s)",
                    url,
                    r.status_code,
                    wait,
                    attempt + 1,
                    max_retries,
                )
                time.sleep(wait)
                continue
            r.raise_for_status()
            return r
        except requests.RequestException as e:
            last_exc = e
            wait = min(2**attempt, 60)
            logger.warning(
                "请求失败 %s: %s，%ss 后重试 (%s/%s)",
                url,
                e,
                wait,
                attempt + 1,
                max_retries,
            )
            time.sleep(wait)
    raise RuntimeError(f"请求最终失败: {url}") from last_exc


def parse_link_next(link_header: str | None) -> str | None:
    """从 Link 头解析 rel=next URL。"""
    if not link_header:
        return None
    for part in link_header.split(","):
        piece = part.strip()
        if 'rel="next"' in piece or "rel=next" in piece:
            start = piece.find("<")
            end = piece.find(">")
            if start >= 0 and end > start:
                return piece[start + 1 : end]
    return None


def should_skip(path: Path, force: bool) -> bool:
    if force:
        return False
    return path.exists() and path.stat().st_size > 0


def result_ok(**kwargs: Any) -> dict[str, Any]:
    out: dict[str, Any] = {"ok": True}
    out.update(kwargs)
    return out


def result_fail(reason: str, **kwargs: Any) -> dict[str, Any]:
    out: dict[str, Any] = {"ok": False, "reason": reason}
    out.update(kwargs)
    return out


def parse_fasta(path: Path) -> list[tuple[str, str]]:
    """返回 [(header_without_gt, sequence), ...]。"""
    records: list[tuple[str, str]] = []
    if not path.exists():
        return records
    header: str | None = None
    seq_parts: list[str] = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq_parts)))
            header = line[1:].strip()
            seq_parts = []
        else:
            seq_parts.append(line.replace(" ", "").upper())
    if header is not None:
        records.append((header, "".join(seq_parts)))
    return records


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    lines: list[str] = []
    for header, seq in records:
        lines.append(f">{header}")
        for i in range(0, len(seq), 60):
            lines.append(seq[i : i + 60])
    atomic_write_text(path, "\n".join(lines) + ("\n" if lines else ""))


def accession_from_header(header: str) -> str:
    """从 FASTA header / BLAST/DIAMOND hit id 提取 accession（容错空段与尾部 |）。"""
    h = (header or "").strip().lstrip(">")
    if not h:
        return ""

    def _first_token(s: str) -> str:
        s = (s or "").strip()
        if not s:
            return ""
        return s.split()[0].split("-")[0].strip()

    if "|" in h:
        parts = [p.strip() for p in h.split("|")]
        # sp|P12345|NAME / pir|T29551| / gb|AAA12345.1|
        if len(parts) >= 2 and parts[0].lower() in {
            "sp",
            "tr",
            "dbj",
            "emb",
            "gb",
            "ref",
            "pdb",
            "pir",
            "prf",
            "tpg",
            "tpe",
            "tpd",
        }:
            acc = _first_token(parts[1])
            if acc:
                return acc
        # 取第一个非空、非 db 标签的字段
        skip = {"sp", "tr", "dbj", "emb", "gb", "ref", "pdb", "pir", "prf", "gi", ""}
        for p in parts:
            if p.lower() in skip:
                continue
            acc = _first_token(p)
            if acc and not acc.lower().startswith("eco:"):
                return acc

    tokens = h.split()
    if not tokens:
        return ""
    token = tokens[0]
    if token.startswith("AF-"):
        mid = token.split("-")
        if len(mid) >= 2 and mid[1]:
            return mid[1]
    return token.split("/")[0].split(".")[0]


def looks_like_uniprot_acc(acc: str) -> bool:
    """粗判是否像 UniProtKB accession（用于 AFDB / UniProt REST）。"""
    import re

    a = (acc or "").strip().upper()
    if not a or "_" in a or "." in a:
        return False
    return bool(
        re.fullmatch(r"[OPQ][0-9][A-Z0-9]{3}[0-9]", a)
        or re.fullmatch(r"[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]", a)
        or re.fullmatch(r"A0A[A-Z0-9]{7,}", a)
    )
