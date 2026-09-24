#!/usr/bin/env bash
# 裸机一键环境：conda 环境 + DLKcat vendor
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_NAME="${ENV_NAME:-gmas_pipeline}"
YML="${ROOT}/environment.yml"

usage() {
  cat <<'EOF'
用法: bash scripts/setup_env.sh [--skip-dlkcat] [--skip-blastdb] [--help]

创建/更新 conda 环境 gmas_pipeline，并拉取 third_party/DLKcat。
可选随后下载 Swiss-Prot BLAST 库。

环境变量:
  ENV_NAME          conda 环境名（默认 gmas_pipeline）
  FORCE_VENDOR=1    强制重拉 DLKcat
  FORCE_BLASTDB=1   强制重下 BLAST 库
EOF
}

SKIP_DLKCAT=0
SKIP_BLASTDB=0
for arg in "$@"; do
  case "$arg" in
    --help|-h) usage; exit 0 ;;
    --skip-dlkcat) SKIP_DLKCAT=1 ;;
    --skip-blastdb) SKIP_BLASTDB=1 ;;
    *) echo "未知参数: $arg"; usage; exit 1 ;;
  esac
done

if ! command -v conda >/dev/null 2>&1; then
  echo "[setup_env] 错误: 未找到 conda，请先安装 Miniconda/Anaconda" >&2
  exit 1
fi

# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"

if conda env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
  echo "[setup_env] 更新已有环境 ${ENV_NAME} ..."
  conda env update -n "${ENV_NAME}" -f "${YML}" --prune
else
  echo "[setup_env] 创建环境 ${ENV_NAME} ..."
  conda env create -n "${ENV_NAME}" -f "${YML}"
fi

conda activate "${ENV_NAME}"
echo "[setup_env] Python: $(which python) ($(python -V 2>&1))"

if [[ "${SKIP_DLKCAT}" != "1" ]]; then
  bash "${ROOT}/scripts/vendor_dlkcat.sh"
fi

if [[ "${SKIP_BLASTDB}" != "1" ]]; then
  echo "[setup_env] 下载 BLAST 库（可 --skip-blastdb 跳过）..."
  bash "${ROOT}/scripts/download_blastdb.sh" || {
    echo "[setup_env] BLAST 库下载失败，可稍后手动: bash scripts/download_blastdb.sh" >&2
  }
fi

echo
echo "[setup_env] 完成。使用方式:"
echo "  conda activate ${ENV_NAME}"
echo "  cd ${ROOT}"
echo "  python -m src.pipeline --keyword gmas --config config.yaml"
echo "文档: docs/DEPLOY.md  docs/USAGE.md"
