#!/usr/bin/env bash
# 下载 / 更新 NCBI Swiss-Prot BLAST 本地库
# 默认写入项目 data/blastdb；可通过 BLASTDB_DIR 指向共享库，例如：
#   BLASTDB_DIR=/mnt/sdb/tmp/Databases/NCBI/blastdb bash scripts/download_blastdb.sh
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DB_DIR="${BLASTDB_DIR:-${ROOT}/data/blastdb}"
DB_NAME="${BLASTDB_NAME:-swissprot}"
NCBI_BASE="${BLASTDB_URL_BASE:-https://ftp.ncbi.nlm.nih.gov/blast/db}"

mkdir -p "${DB_DIR}"
cd "${DB_DIR}"

if [[ -f "${DB_DIR}/${DB_NAME}.pin" ]]; then
  echo "[download_blastdb] 已存在 ${DB_DIR}/${DB_NAME}.pin"
  if [[ "${FORCE_BLASTDB:-0}" != "1" ]]; then
    echo "[download_blastdb] 跳过（FORCE_BLASTDB=1 可强制重下）"
    exit 0
  fi
  echo "[download_blastdb] FORCE_BLASTDB=1，清理旧 ${DB_NAME}.* 后重下"
  rm -f "${DB_DIR}/${DB_NAME}".*
fi

# update_blastdb.pl 仅查找 /usr/bin/curl 与 /usr/local/bin/curl（不认 conda PATH）
ensure_system_curl() {
  if [[ -x /usr/bin/curl || -x /usr/local/bin/curl ]]; then
    return 0
  fi
  local curl_bin
  curl_bin="$(command -v curl || true)"
  if [[ -z "${curl_bin}" ]]; then
    return 1
  fi
  if ln -sf "${curl_bin}" /usr/local/bin/curl 2>/dev/null \
    || ln -sf "${curl_bin}" /usr/bin/curl 2>/dev/null; then
    echo "[download_blastdb] 已链接 ${curl_bin} → 系统 curl 路径"
    return 0
  fi
  return 1
}

fetch_file() {
  local url="$1"
  local out="$2"
  if command -v aria2c >/dev/null 2>&1; then
    aria2c -x 16 -s 16 -k 1M --max-tries=8 --retry-wait=3 \
      --file-allocation=none --allow-overwrite=true --auto-file-renaming=false \
      --console-log-level=notice -o "${out}" "${url}"
  elif command -v curl >/dev/null 2>&1; then
    curl -fL --retry 8 --retry-delay 3 --connect-timeout 30 -o "${out}" "${url}"
  else
    echo "[download_blastdb] 错误: 需要 aria2c 或 curl" >&2
    return 1
  fi
}

download_direct() {
  local archive="${DB_NAME}.tar.gz"
  local md5file="${archive}.md5"
  echo "[download_blastdb] 直接下载 ${NCBI_BASE}/${archive} → ${DB_DIR}"
  echo "[download_blastdb] 约 215 MB，视网络可能需数十分钟"
  fetch_file "${NCBI_BASE}/${archive}" "${archive}"
  fetch_file "${NCBI_BASE}/${md5file}" "${md5file}"
  if command -v md5sum >/dev/null 2>&1; then
    local expected actual
    expected="$(awk '{print $1}' "${md5file}")"
    actual="$(md5sum "${archive}" | awk '{print $1}')"
    if [[ "${expected}" != "${actual}" ]]; then
      echo "[download_blastdb] MD5 校验失败: expected=${expected} actual=${actual}" >&2
      rm -f "${archive}" "${md5file}"
      return 1
    fi
    echo "[download_blastdb] MD5 校验通过"
  fi
  echo "[download_blastdb] 解压 ${archive} ..."
  tar -xzf "${archive}"
  rm -f "${archive}" "${md5file}"
}

if [[ "${USE_UPDATE_BLASTDB:-0}" == "1" ]]; then
  ensure_system_curl || true
  if ! command -v update_blastdb.pl >/dev/null 2>&1; then
    echo "[download_blastdb] 错误: 未找到 update_blastdb.pl" >&2
    exit 1
  fi
  if [[ ! -x /usr/bin/curl && ! -x /usr/local/bin/curl ]]; then
    echo "[download_blastdb] 错误: update_blastdb.pl 需要 /usr/bin/curl 或 /usr/local/bin/curl" >&2
    exit 1
  fi
  echo "[download_blastdb] 经 update_blastdb.pl 下载并解压 ${DB_NAME} → ${DB_DIR}"
  update_blastdb.pl --decompress "${DB_NAME}"
else
  download_direct
fi

if [[ ! -f "${DB_DIR}/${DB_NAME}.pin" ]]; then
  echo "[download_blastdb] 失败: 未生成 ${DB_NAME}.pin" >&2
  exit 1
fi

echo "[download_blastdb] 完成 → ${DB_DIR}/${DB_NAME}"
echo "[download_blastdb] 请将 config 中 blast.database 设为该前缀路径（Docker 默认挂载为 /db/ncbi/blastdb/${DB_NAME}）"
