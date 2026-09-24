#!/usr/bin/env bash
# 增量更新本地 NCBI NR / NT BLAST 库
# 用法:
#   bash scripts/update_ncbi_nr_nt.sh              # 更新 nr + nt
#   bash scripts/update_ncbi_nr_nt.sh nr           # 仅 nr
#   bash scripts/update_ncbi_nr_nt.sh nt           # 仅 nt
# 环境变量:
#   NCBI_DB_ROOT   默认 /mnt/sdb/tmp/Databases/NCBI
#   SOURCE         aws|ncbi|aspera|auto（默认 auto）
#   NUM_THREADS    update_blastdb 并行数（默认 4）
#   ASPERA_RATE    ascp -l 限速，默认 200m
set -euo pipefail

NCBI_ENV="${NCBI_ENV:-/home/wangwenhao/.conda/envs/ncbi}"
export PATH="${NCBI_ENV}/bin:/usr/local/bin:/usr/bin:${PATH}"

NCBI_DB_ROOT="${NCBI_DB_ROOT:-/mnt/sdb/tmp/Databases/NCBI}"
SOURCE="${SOURCE:-auto}"
NUM_THREADS="${NUM_THREADS:-4}"
ASPERA_RATE="${ASPERA_RATE:-200m}"
ASCP_KEY="${ASCP_KEY:-${NCBI_ENV}/etc/asperaweb_id_dsa.openssh}"
LOG_DIR="${LOG_DIR:-/tmp}"
if [[ $# -eq 0 ]]; then
  TARGETS=(nr nt)
else
  TARGETS=("$@")
fi

UPDATE_BLASTDB="$(command -v update_blastdb.pl || true)"
if [[ -z "${UPDATE_BLASTDB}" ]]; then
  echo "[error] 未找到 update_blastdb.pl，请: conda activate ncbi" >&2
  exit 1
fi

if ! command -v aws >/dev/null 2>&1; then
  echo "[warn] 未找到 awscli；建议: conda install -n ncbi -c conda-forge awscli" >&2
fi

# update_blastdb.pl 只认系统路径 curl（无 awscli 时会走 curl）
if [[ ! -x /usr/bin/curl && ! -x /usr/local/bin/curl ]]; then
  curl_bin="$(command -v curl || true)"
  if [[ -n "${curl_bin}" ]]; then
    ln -sf "${curl_bin}" /usr/local/bin/curl 2>/dev/null \
      || sudo ln -sf "${curl_bin}" /usr/local/bin/curl || true
  fi
fi

extract_pending_tarballs() {
  local db_dir="$1"
  local kind="$2" # nr|nt
  local marker="pin"
  [[ "${kind}" == "nt" ]] && marker="nin"
  local f stem
  shopt -s nullglob
  for f in "${db_dir}/${kind}".*.tar.gz; do
    stem="$(basename "${f}" .tar.gz)"
    if [[ -f "${db_dir}/${stem}.${marker}" ]]; then
      continue
    fi
    echo "[extract] ${f}"
    tar -xzf "${f}" -C "${db_dir}"
  done
  shopt -u nullglob
}

pick_source() {
  if [[ "${SOURCE}" != "auto" ]]; then
    echo "${SOURCE}"
    return
  fi
  if command -v aws >/dev/null 2>&1 \
    && aws s3 ls s3://ncbi-blast-databases/ --no-sign-request >/dev/null 2>&1; then
    echo aws
  else
    echo ncbi
  fi
}

update_via_update_blastdb() {
  local db_dir="$1"
  local db_name="$2"
  local src="$3"
  cd "${db_dir}"
  echo "[update_blastdb] db=${db_name} source=${src} threads=${NUM_THREADS} dir=${db_dir}"
  if [[ "${src}" == "aws" ]]; then
    if ! command -v aws >/dev/null 2>&1; then
      echo "[error] source=aws 需要 awscli（conda install -n ncbi -c conda-forge awscli）" >&2
      return 1
    fi
    # AWS 桶为解压后的分卷；--decompress 对 aws 通常可忽略，保留以兼容
    "${UPDATE_BLASTDB}" --source aws --num_threads "${NUM_THREADS}" "${db_name}"
  else
    "${UPDATE_BLASTDB}" --source ncbi --decompress "${db_name}"
  fi
}

update_one() {
  local db_name="$1"
  local db_dir
  case "${db_name}" in
    nr) db_dir="${NCBI_DB_ROOT}/NR" ;;
    nt) db_dir="${NCBI_DB_ROOT}/NT" ;;
    *) echo "[error] 未知库: ${db_name}（仅支持 nr|nt）" >&2; return 1 ;;
  esac
  mkdir -p "${db_dir}"
  local log="${LOG_DIR}/update_${db_name}_$(date +%Y%m%d_%H%M%S).log"
  echo "[start] ${db_name} → ${db_dir}  log=${log}"
  {
    echo "==== $(date) update ${db_name} ===="
    # 仅当走 NCBI tar 源时需要先解开本地未解压包；AWS 直接覆盖分卷文件
    if [[ "$(pick_source)" == "ncbi" ]]; then
      extract_pending_tarballs "${db_dir}" "${db_name}"
    fi
    local src
    src="$(pick_source)"
    echo "[source] ${src}"
    update_via_update_blastdb "${db_dir}" "${db_name}" "${src}"
    echo "==== $(date) done ${db_name} ===="
  } 2>&1 | tee -a "${log}"
  echo "[done] ${db_name} 详见 ${log}"
}

echo "[info] UPDATE_BLASTDB=${UPDATE_BLASTDB}"
echo "[info] NCBI_DB_ROOT=${NCBI_DB_ROOT}"
echo "[info] targets=${TARGETS[*]} SOURCE=${SOURCE}"
echo "[info] 远端 NR/NT 为 TB 级，增量更新可能需数小时到数天；用 tail -f 看日志"

for t in "${TARGETS[@]}"; do
  update_one "${t}"
done
