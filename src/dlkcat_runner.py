"""DLKcat kcat 预测封装。"""

from __future__ import annotations

import logging
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path
from typing import Any

from .io_utils import (
    ROOT,
    accession_from_header,
    atomic_write_text,
    ensure_dir,
    parse_fasta,
    resolve_path,
    result_fail,
    result_ok,
    should_skip,
)
from .substrates import resolve_substrates

logger = logging.getLogger(__name__)

DEFAULT_DLKCAT_ROOT = "third_party/DLKcat/DeeplearningApproach"


def _resolve_dlkcat_root(raw: str | None) -> Path:
    p = Path(raw or DEFAULT_DLKCAT_ROOT)
    if not p.is_absolute():
        p = resolve_path(p)
    return p


def _resolve_python(raw: str | None) -> str:
    """配置为 python/空 时使用当前解释器，便于 Docker/裸机统一。"""
    s = (raw or "").strip()
    if not s or s in {"python", "python3"}:
        return sys.executable
    path = Path(s)
    if path.is_file():
        return str(path)
    which = shutil.which(s)
    return which or s


def ensure_dlkcat_input_dicts(dlkcat_root: Path) -> dict[str, Any]:
    """解压 Data/input.zip（若字典不存在）。"""
    data_dir = dlkcat_root / "Data"
    input_dir = data_dir / "input"
    marker = input_dir / "fingerprint_dict.pickle"
    if marker.exists():
        return result_ok(input_dir=str(input_dir), skipped=True)
    zpath = data_dir / "input.zip"
    if not zpath.exists():
        return result_fail(f"缺少 DLKcat 字典包: {zpath}")
    try:
        ensure_dir(input_dir)
        with zipfile.ZipFile(zpath, "r") as zf:
            zf.extractall(data_dir)
        if not marker.exists():
            return result_fail("解压后未找到 fingerprint_dict.pickle")
        return result_ok(input_dir=str(input_dir))
    except Exception as e:
        logger.exception("解压 DLKcat input.zip 失败")
        return result_fail(f"解压 input.zip 失败: {e}")


def build_dlkcat_input_tsv(
    enzyme_fasta: Path,
    substrates: list[dict[str, Any]],
    out_tsv: Path,
    *,
    max_enzymes: int = 0,
) -> dict[str, Any]:
    """生成官方格式 TSV（另含 Accession 列便于回写）。"""
    records = parse_fasta(enzyme_fasta)
    if not records:
        return result_fail("酶 FASTA 为空", fasta=str(enzyme_fasta))
    if not substrates:
        return result_fail("底物列表为空")

    if max_enzymes and max_enzymes > 0:
        records = records[:max_enzymes]

    lines = ["Accession\tSubstrate Name\tSubstrate SMILES\tProtein Sequence"]
    n = 0
    for header, seq in records:
        acc = accession_from_header(header).upper()
        if not seq:
            continue
        for sub in substrates:
            name = str(sub.get("name") or "")
            smiles = str(sub.get("smiles") or "")
            if not smiles or "." in smiles:
                continue
            lines.append(f"{acc}\t{name}\t{smiles}\t{seq}")
            n += 1
    if n == 0:
        return result_fail("未生成任何酶×底物行")
    atomic_write_text(out_tsv, "\n".join(lines) + "\n")
    return result_ok(path=str(out_tsv), rows=n, enzymes=len(records), substrates=len(substrates))


def _official_input_from_ours(ours: Path, official: Path) -> None:
    """官方脚本只认 3 列：Name / SMILES / Sequence。"""
    lines_out = ["Substrate Name\tSubstrate SMILES\tProtein Sequence"]
    rows = ours.read_text(encoding="utf-8").splitlines()
    for line in rows[1:]:
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        lines_out.append("\t".join(parts[1:4]))
    atomic_write_text(official, "\n".join(lines_out) + "\n")


def run_dlkcat_prediction(
    enzyme_fasta: Path,
    out_dir: Path,
    *,
    cfg: dict[str, Any],
    force: bool = False,
) -> dict[str, Any]:
    """对最终酶集 × 配置底物跑 DLKcat。"""
    ensure_dir(out_dir)
    pred_tsv = out_dir / "dlkcat_predictions.tsv"
    input_tsv = out_dir / "dlkcat_input.tsv"
    if should_skip(pred_tsv, force):
        logger.info("跳过 DLKcat（结果已存在）: %s", pred_tsv)
        return result_ok(predictions=str(pred_tsv), input=str(input_tsv), skipped=True)

    dl_cfg = cfg.get("dlkcat") or {}
    dlkcat_root = _resolve_dlkcat_root(dl_cfg.get("root"))
    py = _resolve_python(dl_cfg.get("python"))
    max_enzymes = int(dl_cfg.get("max_enzymes") or 0)

    if not dlkcat_root.exists():
        return result_fail(
            f"DLKcat 根目录不存在: {dlkcat_root}；请运行 bash scripts/vendor_dlkcat.sh"
        )

    er = ensure_dlkcat_input_dicts(dlkcat_root)
    if not er.get("ok"):
        return er

    cache_dir = Path(
        (cfg.get("paths") or {}).get("substrates")
        or (Path(cfg.get("_project_root") or ROOT) / "data" / "substrates")
    )
    cache_dir = resolve_path(cache_dir) if not Path(cache_dir).is_absolute() else Path(cache_dir)

    sr = resolve_substrates(cfg.get("substrates"), cache_dir, force=False)
    if not sr.get("ok"):
        return sr
    substrates = sr["substrates"]

    br = build_dlkcat_input_tsv(
        enzyme_fasta, substrates, input_tsv, max_enzymes=max_enzymes
    )
    if not br.get("ok"):
        return br

    example_dir = dlkcat_root / "Code" / "example"
    script = example_dir / "prediction_for_input.py"
    if not script.exists():
        return result_fail(f"缺少预测脚本: {script}")

    work_in = example_dir / "_pipeline_input.tsv"
    work_out = example_dir / "output.tsv"
    _official_input_from_ours(input_tsv, work_in)

    bak = example_dir / "output.tsv.bak_pipeline"
    if work_out.exists():
        try:
            shutil.copy2(work_out, bak)
        except OSError:
            pass

    logger.info("运行 DLKcat: rows=%s python=%s root=%s", br.get("rows"), py, dlkcat_root)
    try:
        proc = subprocess.run(
            [py, str(script.name), str(work_in.name)],
            cwd=str(example_dir),
            capture_output=True,
            text=True,
            timeout=int(dl_cfg.get("timeout_sec") or 7200),
        )
    except subprocess.TimeoutExpired:
        return result_fail("DLKcat 预测超时")
    except Exception as e:
        logger.exception("DLKcat 子进程失败")
        return result_fail(f"DLKcat 子进程失败: {e}")

    if proc.returncode != 0:
        logger.error("DLKcat stderr: %s", (proc.stderr or "")[-2000:])
        return result_fail(
            f"DLKcat 退出码 {proc.returncode}",
            stderr=(proc.stderr or "")[-500:],
            stdout=(proc.stdout or "")[-500:],
        )

    if not work_out.exists():
        return result_fail("DLKcat 未生成 output.tsv")

    our_rows = [
        ln.split("\t")
        for ln in input_tsv.read_text(encoding="utf-8").splitlines()[1:]
        if ln.strip()
    ]
    out_lines = work_out.read_text(encoding="utf-8").splitlines()
    pred_lines = [
        "accession\tsubstrate_name\tsubstrate_smiles\tseq_length\tkcat_per_s\tstatus"
    ]
    out_data = [ln.split("\t") for ln in out_lines[1:] if ln.strip()]

    for i, parts in enumerate(our_rows):
        acc = parts[0]
        name = parts[1]
        smiles = parts[2]
        seq = parts[3]
        seq_len = str(len(seq))
        if i < len(out_data) and len(out_data[i]) >= 4:
            kcat = out_data[i][3]
            status = "ok" if kcat not in {"", "None", "none"} else "failed"
        else:
            kcat = "None"
            status = "missing_output"
        pred_lines.append(f"{acc}\t{name}\t{smiles}\t{seq_len}\t{kcat}\t{status}")

    atomic_write_text(pred_tsv, "\n".join(pred_lines) + "\n")
    shutil.copy2(work_out, out_dir / "dlkcat_output_raw.tsv")

    ok_n = sum(1 for ln in pred_lines[1:] if ln.endswith("\tok"))
    logger.info("DLKcat 完成: 总行 %s，成功 %s → %s", len(pred_lines) - 1, ok_n, pred_tsv)
    return result_ok(
        predictions=str(pred_tsv),
        input=str(input_tsv),
        rows=len(pred_lines) - 1,
        ok_rows=ok_n,
        substrates=len(substrates),
    )
