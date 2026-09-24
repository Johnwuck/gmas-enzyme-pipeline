#!/usr/bin/env bash
# 将 SysBioChalmers/DLKcat 拉取到 third_party/DLKcat
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEST="${ROOT}/third_party/DLKcat"
REPO_URL="${DLKCAT_REPO_URL:-https://github.com/SysBioChalmers/DLKcat.git}"
REF="${DLKCAT_REF:-master}"
# 若本机已有完整 DLKcat: LOCAL_DLKCAT=/path/to/DLKcat VENDOR_MODE=link bash scripts/vendor_dlkcat.sh
LOCAL_DLKCAT="${LOCAL_DLKCAT:-}"

mkdir -p "${ROOT}/third_party"

copy_or_link_local() {
  local src="$1"
  echo "[vendor_dlkcat] 使用本地副本: ${src} → ${DEST}"
  rm -rf "${DEST}"
  if [[ "${VENDOR_MODE:-copy}" == "link" ]]; then
    ln -sfn "${src}" "${DEST}"
  else
    mkdir -p "${DEST}"
    cp -a "${src}/." "${DEST}/"
  fi
}

if [[ -n "${LOCAL_DLKCAT}" && -d "${LOCAL_DLKCAT}/DeeplearningApproach" ]]; then
  if [[ "${FORCE_VENDOR:-0}" == "1" || ! -e "${DEST}" ]]; then
    copy_or_link_local "${LOCAL_DLKCAT}"
  fi
elif [[ -d "${DEST}/.git" || -L "${DEST}" || -d "${DEST}/DeeplearningApproach" ]]; then
  echo "[vendor_dlkcat] 已存在 ${DEST}（FORCE_VENDOR=1 可强制重拉）"
  if [[ "${FORCE_VENDOR:-0}" == "1" ]]; then
    rm -rf "${DEST}"
  fi
fi

if [[ ! -d "${DEST}/DeeplearningApproach" ]]; then
  echo "[vendor_dlkcat] clone ${REPO_URL} @ ${REF} → ${DEST}"
  git clone --depth 1 --branch "${REF}" "${REPO_URL}" "${DEST}" 2>/dev/null || \
    git clone --depth 1 "${REPO_URL}" "${DEST}"
fi

INPUT_ZIP="${DEST}/DeeplearningApproach/Data/input.zip"
INPUT_DIR="${DEST}/DeeplearningApproach/Data/input"
MARKER="${INPUT_DIR}/fingerprint_dict.pickle"

if [[ ! -f "${MARKER}" ]]; then
  if [[ -f "${INPUT_ZIP}" ]]; then
    echo "[vendor_dlkcat] 解压 input.zip ..."
    if ! unzip -qo "${INPUT_ZIP}" -d "${DEST}/DeeplearningApproach/Data"; then
      echo "[vendor_dlkcat] 警告: input.zip 解压失败，请检查文件完整性" >&2
    fi
  else
    echo "[vendor_dlkcat] 警告: 缺少 ${INPUT_ZIP}" >&2
  fi
fi

WEIGHT="${DEST}/DeeplearningApproach/Results/output/all--radius2--ngram3--dim20--layer_gnn3--window11--layer_cnn3--layer_output3--lr1e-3--lr_decay0.5--decay_interval10--weight_decay1e-6--iteration50"
if [[ ! -f "${WEIGHT}" ]]; then
  echo "[vendor_dlkcat] 警告: 未找到模型权重: ${WEIGHT}" >&2
else
  echo "[vendor_dlkcat] 权重就绪"
fi

if [[ -f "${MARKER}" ]]; then
  echo "[vendor_dlkcat] 字典就绪"
else
  echo "[vendor_dlkcat] 警告: 字典未就绪 ${MARKER}" >&2
fi

echo "[vendor_dlkcat] 完成 → ${DEST}"
