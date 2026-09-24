#!/usr/bin/env bash
# Docker 容器入口：检查依赖后转发到 pipeline CLI
set -euo pipefail

ROOT="${PIPELINE_ROOT:-/app}"
cd "${ROOT}"

if [[ ! -d "${ROOT}/third_party/DLKcat/DeeplearningApproach" ]]; then
  echo "[entrypoint] 拉取 DLKcat ..."
  bash "${ROOT}/scripts/vendor_dlkcat.sh"
fi

# 从 config 读取 blast.database；并探测常见库前缀
blast_db=""
if [[ -f "${ROOT}/config.yaml" ]] && command -v python >/dev/null 2>&1; then
  blast_db="$(python - <<'PY' 2>/dev/null || true
import yaml
from pathlib import Path
p = Path("/app/config.yaml")
cfg = yaml.safe_load(p.read_text(encoding="utf-8")) or {}
print((cfg.get("blast") or {}).get("database") or "")
PY
)"
fi

candidates=()
if [[ -n "${blast_db}" ]]; then
  candidates+=("${blast_db}.dmnd" "${blast_db}.pal" "${blast_db}.pin" "${blast_db}.nal" "${blast_db}")
fi
candidates+=(
  "${NCBI_DB_ROOT:-/db/ncbi}/diamond/nr.dmnd"
  "${NCBI_DB_ROOT:-/db/ncbi}/NR/nr.pal"
  "${NCBI_DB_ROOT:-/db/ncbi}/blastdb/swissprot.pin"
  "${ROOT}/data/blastdb/swissprot.pin"
)

found=""
for p in "${candidates[@]}"; do
  if [[ -e "${p}" ]]; then
    found="${p}"
    break
  fi
done

if [[ -z "${found}" ]]; then
  echo "[entrypoint] 提示: 未检测到本地比对库"
  echo "  DIAMOND NR: docker compose --profile setup run --rm build_diamond_nr"
  echo "  Swiss-Prot: docker compose --profile setup run --rm download_blastdb"
  echo "  或将 blast.mode 改为 remote（仅 blastp）"
else
  echo "[entrypoint] 比对库就绪: ${found}${blast_db:+ (config: ${blast_db})}"
fi

if [[ $# -eq 0 ]]; then
  echo "用法示例:"
  echo "  docker compose run --rm pipeline --keyword gmas --step all"
  echo "  docker compose run --rm pipeline --keyword gmas --step dlkcat"
  python -m src.pipeline --help
  exit 0
fi

exec python -m src.pipeline --config "${ROOT}/config.yaml" "$@"
