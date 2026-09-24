"""UniProtKB 搜索与序列下载。"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

from .io_utils import (
    accession_from_header,
    atomic_write_text,
    ensure_dir,
    parse_fasta,
    parse_link_next,
    request_get,
    result_fail,
    result_ok,
    should_skip,
    write_fasta,
)

logger = logging.getLogger(__name__)

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_STREAM = "https://rest.uniprot.org/uniprotkb/stream"
UNIPROT_ENTRY = "https://rest.uniprot.org/uniprotkb"


def _extract_gene(entry: dict[str, Any]) -> str:
    genes = entry.get("genes") or []
    names: list[str] = []
    for g in genes:
        if isinstance(g, dict):
            gn = g.get("geneName") or {}
            if isinstance(gn, dict) and gn.get("value"):
                names.append(str(gn["value"]))
            for syn in g.get("synonyms") or []:
                if isinstance(syn, dict) and syn.get("value"):
                    names.append(str(syn["value"]))
    return ";".join(names)


def _extract_protein_name(entry: dict[str, Any]) -> str:
    pd = entry.get("proteinDescription") or {}
    rec = pd.get("recommendedName") or {}
    full = rec.get("fullName") or {}
    if isinstance(full, dict) and full.get("value"):
        return str(full["value"])
    alts = pd.get("submissionNames") or pd.get("alternativeNames") or []
    for alt in alts:
        if isinstance(alt, dict):
            fn = (alt.get("fullName") or {}).get("value")
            if fn:
                return str(fn)
    return ""


def _extract_ec(entry: dict[str, Any]) -> str:
    ecs: list[str] = []
    pd = entry.get("proteinDescription") or {}
    rec = pd.get("recommendedName") or {}
    for ec in rec.get("ecNumbers") or []:
        if isinstance(ec, dict) and ec.get("value"):
            ecs.append(str(ec["value"]))
    return ";".join(ecs)


def _extract_organism(entry: dict[str, Any]) -> str:
    org = entry.get("organism") or {}
    return str(org.get("scientificName") or org.get("commonName") or "")


def _extract_sequence(entry: dict[str, Any]) -> str:
    seq = entry.get("sequence") or {}
    return str(seq.get("value") or "").replace(" ", "").upper()


def _entry_to_row(entry: dict[str, Any]) -> dict[str, str]:
    acc = str(entry.get("primaryAccession") or "")
    seq = _extract_sequence(entry)
    return {
        "accession": acc,
        "gene": _extract_gene(entry),
        "protein_name": _extract_protein_name(entry),
        "organism": _extract_organism(entry),
        "length": str(len(seq) if seq else entry.get("sequence", {}).get("length") or ""),
        "ec": _extract_ec(entry),
        "reviewed": "reviewed" if entry.get("entryType") == "UniProtKB reviewed (Swiss-Prot)" else "unreviewed",
        "sequence": seq,
    }


def search_uniprotkb(
    query: str,
    *,
    page_size: int = 500,
    timeout: float = 60,
    max_retries: int = 5,
) -> dict[str, Any]:
    """分页搜索 UniProtKB，返回全部 JSON entries。"""
    rows: list[dict[str, Any]] = []
    url: str | None = UNIPROT_SEARCH
    params: dict[str, Any] | None = {
        "query": query,
        "format": "json",
        "size": page_size,
        "fields": "accession,id,gene_names,protein_name,organism_name,length,ec,reviewed,sequence",
    }
    # fields + format=json 可能不返回完整 sequence 对象；改用不限 fields 的完整 entry
    params = {"query": query, "format": "json", "size": page_size}

    while url:
        try:
            r = request_get(url, params=params, timeout=timeout, max_retries=max_retries)
        except Exception as e:
            logger.exception("UniProt 搜索失败")
            return result_fail(f"UniProt 搜索失败: {e}", count=len(rows))
        data = r.json()
        batch = data.get("results") or []
        rows.extend(batch)
        logger.info("UniProt 已获取 %s 条", len(rows))
        next_url = parse_link_next(r.headers.get("Link"))
        url = next_url
        params = None  # next URL 已含参数
    return result_ok(entries=rows, count=len(rows))


def stream_fasta(
    query: str,
    out_fasta: Path,
    *,
    timeout: float = 120,
    max_retries: int = 5,
) -> dict[str, Any]:
    """用 stream 端点下载 FASTA。"""
    try:
        r = request_get(
            UNIPROT_STREAM,
            params={"query": query, "format": "fasta"},
            timeout=timeout,
            max_retries=max_retries,
        )
        text = r.text
        if not text.strip():
            return result_fail("UniProt stream 返回空 FASTA", path=str(out_fasta))
        atomic_write_text(out_fasta, text if text.endswith("\n") else text + "\n")
        n = len(parse_fasta(out_fasta))
        return result_ok(path=str(out_fasta), count=n)
    except Exception as e:
        logger.exception("UniProt FASTA stream 失败")
        return result_fail(f"FASTA 下载失败: {e}", path=str(out_fasta))


_TSV_HEADER_MAP = {
    "Entry": "accession",
    "Gene Names": "gene_names",
    "Protein names": "protein_name",
    "Organism": "organism_name",
    "Length": "length",
    "EC number": "ec",
    "Reviewed": "reviewed",
}


def _normalize_uniprot_tsv(text: str) -> str:
    """将 UniProt stream 表头规范为稳定英文列名。"""
    lines = text.splitlines()
    if not lines:
        return text
    cols = lines[0].split("\t")
    mapped = [_TSV_HEADER_MAP.get(c, c.strip().lower().replace(" ", "_")) for c in cols]
    lines[0] = "\t".join(mapped)
    return "\n".join(lines) + "\n"


def stream_tsv_metadata(
    query: str,
    out_tsv: Path,
    *,
    timeout: float = 120,
    max_retries: int = 5,
) -> dict[str, Any]:
    fields = "accession,gene_names,protein_name,organism_name,length,ec,reviewed"
    try:
        r = request_get(
            UNIPROT_STREAM,
            params={"query": query, "format": "tsv", "fields": fields},
            timeout=timeout,
            max_retries=max_retries,
        )
        text = r.text
        if not text.strip():
            return result_fail("UniProt stream 返回空 TSV", path=str(out_tsv))
        normalized = _normalize_uniprot_tsv(text)
        atomic_write_text(out_tsv, normalized)
        lines = [ln for ln in normalized.splitlines() if ln.strip()]
        count = max(0, len(lines) - 1)
        return result_ok(path=str(out_tsv), count=count)
    except Exception as e:
        logger.exception("UniProt TSV stream 失败")
        return result_fail(f"TSV 下载失败: {e}", path=str(out_tsv))


def fetch_sequences_by_accessions(
    accessions: list[str],
    out_fasta: Path,
    *,
    timeout: float = 60,
    max_retries: int = 5,
    batch_size: int = 100,
) -> dict[str, Any]:
    """按 accession 批量拉取 FASTA（去重）。"""
    uniq = sorted({a.strip().upper() for a in accessions if a and a.strip()})
    if not uniq:
        write_fasta(out_fasta, [])
        return result_ok(path=str(out_fasta), count=0)

    records: list[tuple[str, str]] = []
    for i in range(0, len(uniq), batch_size):
        batch = uniq[i : i + batch_size]
        query = " OR ".join(f"accession:{a}" for a in batch)
        try:
            r = request_get(
                UNIPROT_STREAM,
                params={"query": query, "format": "fasta"},
                timeout=timeout,
                max_retries=max_retries,
            )
            tmp = out_fasta.parent / f".batch_{i}.fasta"
            atomic_write_text(tmp, r.text)
            for header, seq in parse_fasta(tmp):
                records.append((header, seq))
            tmp.unlink(missing_ok=True)
        except Exception as e:
            logger.warning("批量拉取 accession 失败 (%s): %s", batch[:3], e)
            # 回退单条
            for acc in batch:
                try:
                    r = request_get(
                        f"{UNIPROT_ENTRY}/{acc}.fasta",
                        timeout=timeout,
                        max_retries=max_retries,
                    )
                    for header, seq in _parse_fasta_text(r.text):
                        records.append((header, seq))
                except Exception as e2:
                    logger.warning("单条拉取失败 %s: %s", acc, e2)

    # 按 accession 去重
    seen: set[str] = set()
    deduped: list[tuple[str, str]] = []
    for header, seq in records:
        acc = accession_from_header(header).upper()
        if acc in seen:
            continue
        seen.add(acc)
        deduped.append((header, seq))
    write_fasta(out_fasta, deduped)
    return result_ok(path=str(out_fasta), count=len(deduped), requested=len(uniq))


def _parse_fasta_text(text: str) -> list[tuple[str, str]]:
    from tempfile import NamedTemporaryFile

    with NamedTemporaryFile("w", suffix=".fasta", delete=False, encoding="utf-8") as f:
        f.write(text)
        name = f.name
    try:
        return parse_fasta(Path(name))
    finally:
        Path(name).unlink(missing_ok=True)


def download_seed(
    query: str,
    raw_dir: Path,
    *,
    force: bool = False,
    timeout: float = 120,
    max_retries: int = 5,
    page_size: int = 500,
    seed_prefix: str = "seed",
) -> dict[str, Any]:
    """搜索并下载 seed FASTA + metadata TSV。"""
    ensure_dir(raw_dir)
    fasta_path = raw_dir / f"{seed_prefix}_seed.fasta"
    tsv_path = raw_dir / f"{seed_prefix}_uniprot.tsv"

    if should_skip(fasta_path, force) and should_skip(tsv_path, force):
        n = len(parse_fasta(fasta_path))
        logger.info("跳过 UniProt 下载（已存在），seed=%s", n)
        if n == 0:
            return result_fail("已有 seed FASTA 为空", fasta=str(fasta_path), tsv=str(tsv_path))
        return result_ok(fasta=str(fasta_path), tsv=str(tsv_path), count=n, skipped=True)

    # 先 stream FASTA（全量序列）
    fr = stream_fasta(query, fasta_path, timeout=timeout, max_retries=max_retries)
    if not fr.get("ok"):
        # 回退分页 JSON
        logger.warning("stream FASTA 失败，回退分页 JSON: %s", fr.get("reason"))
        sr = search_uniprotkb(
            query, page_size=page_size, timeout=timeout, max_retries=max_retries
        )
        if not sr.get("ok"):
            return sr
        entries = sr.get("entries") or []
        if not entries:
            return result_fail("UniProt 搜索结果为空", query=query)
        records = []
        tsv_lines = ["accession\tgene_names\tprotein_name\torganism_name\tlength\tec\treviewed"]
        for e in entries:
            row = _entry_to_row(e)
            if not row["accession"] or not row["sequence"]:
                continue
            records.append((f"{row['accession']}|{row['gene']}|{row['organism']}", row["sequence"]))
            tsv_lines.append(
                "\t".join(
                    [
                        row["accession"],
                        row["gene"],
                        row["protein_name"],
                        row["organism"],
                        row["length"],
                        row["ec"],
                        row["reviewed"],
                    ]
                )
            )
        write_fasta(fasta_path, records)
        atomic_write_text(tsv_path, "\n".join(tsv_lines) + "\n")
        if not records:
            return result_fail("UniProt 搜索结果为空", query=query)
        return result_ok(fasta=str(fasta_path), tsv=str(tsv_path), count=len(records))

    tr = stream_tsv_metadata(query, tsv_path, timeout=timeout, max_retries=max_retries)
    if not tr.get("ok"):
        logger.warning("TSV metadata 下载失败: %s（FASTA 已就绪）", tr.get("reason"))

    count = int(fr.get("count") or 0)
    if count == 0:
        return result_fail("UniProt 搜索结果为空", query=query, fasta=str(fasta_path))
    return result_ok(fasta=str(fasta_path), tsv=str(tsv_path), count=count)
