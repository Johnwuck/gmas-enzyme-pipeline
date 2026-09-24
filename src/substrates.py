"""底物下载与缓存（PubChem）。"""

from __future__ import annotations

import json
import logging
import re
from pathlib import Path
from typing import Any

from .io_utils import atomic_write_text, ensure_dir, request_get, result_fail, result_ok

logger = logging.getLogger(__name__)

PUBCHEM_PROP = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/"
    "ConnectivitySMILES,CanonicalSMILES,IUPACName/JSON"
)


def _slug(name: str) -> str:
    s = re.sub(r"[^A-Za-z0-9._-]+", "_", name.strip())
    return s.strip("_") or "substrate"


def fetch_smiles_from_pubchem(name: str, *, timeout: float = 30) -> dict[str, Any]:
    """按化合物名称从 PubChem 取 ConnectivitySMILES。"""
    url = PUBCHEM_PROP.format(name=requests_quote(name))
    try:
        r = request_get(url, timeout=timeout, max_retries=3)
        props = (r.json().get("PropertyTable") or {}).get("Properties") or []
        if not props:
            return result_fail(f"PubChem 无结果: {name}")
        p0 = props[0]
        smiles = (
            p0.get("ConnectivitySMILES")
            or p0.get("SMILES")
            or p0.get("CanonicalSMILES")
            or ""
        )
        if not smiles:
            return result_fail(f"PubChem 无 SMILES: {name}", raw=p0)
        return result_ok(
            name=name,
            smiles=str(smiles),
            cid=p0.get("CID"),
            iupac=p0.get("IUPACName") or "",
            source="pubchem",
        )
    except Exception as e:
        logger.exception("PubChem 查询失败: %s", name)
        return result_fail(f"PubChem 查询失败: {e}", name=name)


def requests_quote(name: str) -> str:
    from urllib.parse import quote

    return quote(name, safe="")


def resolve_substrates(
    cfg_list: list[dict[str, Any]] | None,
    cache_dir: Path,
    *,
    force: bool = False,
) -> dict[str, Any]:
    """解析配置中的底物列表，缺 SMILES 时下载并缓存 JSON。"""
    if not cfg_list:
        return result_fail("未配置 substrates，请在 config.yaml 中添加底物")

    ensure_dir(cache_dir)
    resolved: list[dict[str, Any]] = []
    for item in cfg_list:
        if not isinstance(item, dict):
            continue
        name = str(item.get("name") or "").strip()
        if not name:
            continue
        smiles = str(item.get("smiles") or "").strip()
        cid = item.get("cid")
        cache_path = cache_dir / f"{_slug(name)}.json"

        if smiles and not force:
            rec = {
                "name": name,
                "smiles": smiles,
                "cid": cid,
                "source": "config",
            }
            atomic_write_text(cache_path, json.dumps(rec, ensure_ascii=False, indent=2) + "\n")
            resolved.append(rec)
            continue

        if cache_path.exists() and not force and not smiles:
            try:
                rec = json.loads(cache_path.read_text(encoding="utf-8"))
                if rec.get("smiles"):
                    resolved.append(rec)
                    continue
            except Exception:
                pass

        fr = fetch_smiles_from_pubchem(name)
        if not fr.get("ok"):
            logger.warning("底物 %s 下载失败: %s", name, fr.get("reason"))
            continue
        rec = {
            "name": name,
            "smiles": fr["smiles"],
            "cid": fr.get("cid") or cid,
            "iupac": fr.get("iupac") or "",
            "source": "pubchem",
        }
        atomic_write_text(cache_path, json.dumps(rec, ensure_ascii=False, indent=2) + "\n")
        resolved.append(rec)
        logger.info("已缓存底物 %s → %s (%s)", name, rec["smiles"], cache_path)

    if not resolved:
        return result_fail("无可用底物（全部解析失败）")
    return result_ok(substrates=resolved, count=len(resolved), cache_dir=str(cache_dir))
