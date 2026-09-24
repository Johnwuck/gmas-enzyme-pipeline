#!/usr/bin/env bash
# 下载 Foldseek 格式的 AlphaFold DB，供「序列搜 AFDB」使用。
#
# 默认走 https://opendata.mmseqs.org/foldseek/（官方 workers.dev 在部分网络不可达）。
# 默认库: Alphafold/Swiss-Prot（~1.5 GiB，覆盖较好、下载快）
# 可选: Alphafold/UniProt50（afdb50，~114 GiB）| Alphafold/Proteome
#
# 用法:
#   bash scripts/download_foldseek_afdb.sh
#   FOLDSEEK_AFDB_NAME=Alphafold/UniProt50 bash scripts/download_foldseek_afdb.sh
#   docker compose --profile setup run --rm download_foldseek_afdb
set -euo pipefail

if [[ -z "${NCBI_DB_ROOT:-}" ]]; then
  if [[ -d /db/ncbi ]]; then
    NCBI_DB_ROOT=/db/ncbi
  else
    NCBI_DB_ROOT=/mnt/sdb/tmp/Databases/NCBI
  fi
fi

OUT_DIR="${FOLDSEEK_AFDB_DIR:-${NCBI_DB_ROOT}/foldseek}"
DB_NAME="${FOLDSEEK_AFDB_NAME:-Alphafold/Swiss-Prot}"
MIRROR="${FOLDSEEK_AFDB_MIRROR:-https://opendata.mmseqs.org/foldseek}"
THREADS="${FOLDSEEK_THREADS:-${NUM_THREADS:-48}}"
TMP_DIR="${OUT_DIR}/tmp_download"

case "${DB_NAME}" in
  *Swiss-Prot*)
    PREFIX="${OUT_DIR}/afdb_swissprot"
    TARBALL="afdb_swissprot.tar.gz"
    ;;
  *UniProt50-minimal*|*UniProt50*)
    PREFIX="${OUT_DIR}/afdb_uniprot50"
    TARBALL="afdb50.tar.gz"
    ;;
  *Proteome*)
    PREFIX="${OUT_DIR}/afdb_proteome"
    TARBALL="afdb_proteome.tar.gz"
    ;;
  *UniProt*)
    PREFIX="${OUT_DIR}/afdb_uniprot"
    TARBALL="afdb.tar.gz"
    ;;
  *)
    PREFIX="${OUT_DIR}/afdb"
    TARBALL="afdb_swissprot.tar.gz"
    ;;
esac

mkdir -p "${OUT_DIR}" "${TMP_DIR}"

if [[ -f "${PREFIX}.dbtype" || -f "${PREFIX}" ]]; then
  echo "[download_foldseek_afdb] 已存在 ${PREFIX}"
  if [[ "${FORCE_FOLDSEEK_AFDB:-0}" != "1" ]]; then
    echo "[download_foldseek_afdb] 跳过（FORCE_FOLDSEEK_AFDB=1 可强制重下）"
    ls -lah "${PREFIX}"* 2>/dev/null | head -20
    exit 0
  fi
  echo "[download_foldseek_afdb] FORCE=1，清理旧库"
  rm -f "${PREFIX}"*
fi

URL="${MIRROR}/${TARBALL}"
ARCHIVE="${TMP_DIR}/${TARBALL}"
echo "[download_foldseek_afdb] name=${DB_NAME} url=${URL} → ${PREFIX}"

# 1) 优先镜像直链（aria2c / wget / curl）
download_ok=0
if command -v aria2c >/dev/null 2>&1; then
  if aria2c -c -x 8 -s 8 -d "${TMP_DIR}" -o "${TARBALL}" "${URL}"; then
    download_ok=1
  fi
fi
if [[ "${download_ok}" != "1" ]] && command -v wget >/dev/null 2>&1; then
  if wget -c -O "${ARCHIVE}" "${URL}"; then
    download_ok=1
  fi
fi
if [[ "${download_ok}" != "1" ]] && command -v curl >/dev/null 2>&1; then
  if curl -L --retry 5 --retry-delay 5 -C - -o "${ARCHIVE}" "${URL}"; then
    download_ok=1
  fi
fi

if [[ "${download_ok}" == "1" && -s "${ARCHIVE}" ]]; then
  echo "[download_foldseek_afdb] 解压 ${ARCHIVE}"
  EXTRACT_DIR="${TMP_DIR}/extract_$$"
  mkdir -p "${EXTRACT_DIR}"
  tar -xzf "${ARCHIVE}" -C "${EXTRACT_DIR}"
  # tarball 内通常已是带前缀的 DB 文件，或单一子目录
  shopt -s nullglob
  # 找 .dbtype
  mapfile -t DBTYPES < <(find "${EXTRACT_DIR}" -name "*.dbtype" | head -20)
  if [[ ${#DBTYPES[@]} -eq 0 ]]; then
    # 无扩展名源文件
    mapfile -t DBTYPES < <(find "${EXTRACT_DIR}" -type f -name "*.index" | head -5)
  fi
  if [[ ${#DBTYPES[@]} -eq 0 ]]; then
    echo "[download_foldseek_afdb] 错误: 解压后未找到 Foldseek DB 文件" >&2
    ls -laR "${EXTRACT_DIR}" | head -50 >&2
    exit 1
  fi
  # 取第一个 .dbtype 的 stem 作为源前缀
  sample="${DBTYPES[0]}"
  src_prefix="${sample%.dbtype}"
  src_dir=$(dirname "${src_prefix}")
  src_base=$(basename "${src_prefix}")
  echo "[download_foldseek_afdb] 源前缀: ${src_prefix}"
  for f in "${src_dir}/${src_base}"*; do
    bn=$(basename "$f")
    suffix="${bn#${src_base}}"
    dest="${PREFIX}${suffix}"
    echo "  mv ${bn} → $(basename "${dest}")"
    mv -f "$f" "${dest}"
  done
  rm -rf "${EXTRACT_DIR}"
  # 可选保留 tarball 以便校验；默认删释放空间
  if [[ "${KEEP_FOLDSEEK_TARBALL:-0}" != "1" ]]; then
    rm -f "${ARCHIVE}"
  fi
else
  echo "[download_foldseek_afdb] 镜像直链失败，回退 foldseek databases（需可达 workers.dev）"
  if ! command -v foldseek >/dev/null 2>&1; then
    echo "[download_foldseek_afdb] 错误: 未找到 foldseek" >&2
    exit 1
  fi
  foldseek databases "${DB_NAME}" "${PREFIX}" "${TMP_DIR}" --threads "${THREADS}"
fi

echo "[download_foldseek_afdb] 完成:"
ls -lah "${PREFIX}"* 2>/dev/null | head -30
echo "[download_foldseek_afdb] config: foldseek.afdb_database: ${PREFIX}"
