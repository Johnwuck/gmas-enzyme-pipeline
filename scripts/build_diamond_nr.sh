#!/usr/bin/env bash
# 构建 DIAMOND NR 库（.dmnd）。DIAMOND 不能直接读 BLAST .pal/.psq。
#
# 默认输出: ${NCBI_DB_ROOT}/diamond/nr.dmnd
# 数据来源优先级:
#   1) DIAMOND_NR_FASTA 指向的 faa / faa.gz / nr.gz
#   2) ${DIAMOND_DB_DIR}/nr.faa[.gz] 或 ${NCBI_DB_ROOT}/FASTA/nr.gz
#   3) 从已有 BLAST NR 用 blastdbcmd 导出并管道进 makedb（推荐本机）
#   4) 从 NCBI FTP 下载 blast/db/FASTA/nr.gz
#
# 用法:
#   bash scripts/build_diamond_nr.sh
#   NCBI_DB_ROOT=/mnt/sdb/tmp/Databases/NCBI bash scripts/build_diamond_nr.sh
#   docker compose --profile setup run --rm build_diamond_nr
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# 宿主机默认共享库；容器内 compose 会设 NCBI_DB_ROOT=/db/ncbi
if [[ -z "${NCBI_DB_ROOT:-}" ]]; then
  if [[ -d /db/ncbi ]]; then
    NCBI_DB_ROOT=/db/ncbi
  else
    NCBI_DB_ROOT=/mnt/sdb/tmp/Databases/NCBI
  fi
fi
DIAMOND_DB_DIR="${DIAMOND_DB_DIR:-${NCBI_DB_ROOT}/diamond}"
PREFIX="${DIAMOND_DB_DIR}/nr"
THREADS="${DIAMOND_THREADS:-${NUM_THREADS:-48}}"
BLAST_NR="${BLAST_NR_DB:-${NCBI_DB_ROOT}/NR/nr}"
NCBI_FASTA_URL="${NCBI_NR_FASTA_URL:-https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz}"

mkdir -p "${DIAMOND_DB_DIR}"

if [[ -f "${PREFIX}.dmnd" ]]; then
  echo "[build_diamond_nr] 已存在 ${PREFIX}.dmnd"
  if [[ "${FORCE_DIAMOND_NR:-0}" != "1" ]]; then
    echo "[build_diamond_nr] 跳过（FORCE_DIAMOND_NR=1 可强制重建）"
    ls -lah "${PREFIX}.dmnd"
    exit 0
  fi
  echo "[build_diamond_nr] FORCE_DIAMOND_NR=1，删除旧库后重建"
  rm -f "${PREFIX}.dmnd"
fi

if ! command -v diamond >/dev/null 2>&1; then
  echo "[build_diamond_nr] 错误: 未找到 diamond，请 conda install -c bioconda diamond" >&2
  exit 1
fi

echo "[build_diamond_nr] diamond=$(command -v diamond) threads=${THREADS}"
echo "[build_diamond_nr] 输出前缀: ${PREFIX}"

resolve_fasta() {
  if [[ -n "${DIAMOND_NR_FASTA:-}" && -f "${DIAMOND_NR_FASTA}" ]]; then
    echo "${DIAMOND_NR_FASTA}"
    return 0
  fi
  local cand
  for cand in \
    "${DIAMOND_DB_DIR}/nr.faa.gz" \
    "${DIAMOND_DB_DIR}/nr.faa" \
    "${DIAMOND_DB_DIR}/nr.gz" \
    "${NCBI_DB_ROOT}/FASTA/nr.gz" \
    "${NCBI_DB_ROOT}/nr.gz" \
    "${NCBI_DB_ROOT}/nr.faa.gz"
  do
    if [[ -f "${cand}" ]]; then
      echo "${cand}"
      return 0
    fi
  done
  return 1
}

makedb_from_fasta() {
  local fasta="$1"
  echo "[build_diamond_nr] diamond makedb --in ${fasta} -d ${PREFIX}"
  diamond makedb --in "${fasta}" -d "${PREFIX}" --threads "${THREADS}"
}

makedb_from_blast_nr() {
  if ! command -v blastdbcmd >/dev/null 2>&1; then
    return 1
  fi
  if [[ ! -f "${BLAST_NR}.pal" && ! -f "${BLAST_NR}.pin" ]]; then
    echo "[build_diamond_nr] BLAST NR 不可用: ${BLAST_NR}"
    return 1
  fi
  echo "[build_diamond_nr] 从 BLAST NR 管道导出 → diamond makedb（耗时长，无需另下 nr.gz）"
  echo "[build_diamond_nr] BLAST_NR=${BLAST_NR}"
  # 部分 diamond 版本对 stdin 支持不稳；写临时 faa.gz 更稳妥但占盘。优先管道。
  if blastdbcmd -db "${BLAST_NR}" -entry all -outfmt %f 2>/tmp/blastdbcmd_nr.err \
    | diamond makedb --in - -d "${PREFIX}" --threads "${THREADS}"
  then
    return 0
  fi
  echo "[build_diamond_nr] 管道建库失败，见 /tmp/blastdbcmd_nr.err；尝试落盘再 makedb" >&2
  local tmp_faa="${DIAMOND_DB_DIR}/nr_from_blast.faa"
  echo "[build_diamond_nr] blastdbcmd → ${tmp_faa}（体积可能数百 GB）"
  blastdbcmd -db "${BLAST_NR}" -entry all -outfmt %f -out "${tmp_faa}"
  diamond makedb --in "${tmp_faa}" -d "${PREFIX}" --threads "${THREADS}"
  if [[ "${KEEP_NR_FASTA:-0}" != "1" ]]; then
    rm -f "${tmp_faa}"
  fi
}

download_nr_gz() {
  local out="${DIAMOND_DB_DIR}/nr.gz"
  echo "[build_diamond_nr] 下载 ${NCBI_FASTA_URL} → ${out}"
  if command -v aria2c >/dev/null 2>&1; then
    aria2c -x 16 -s 16 -k 1M --max-tries=8 --retry-wait=3 \
      --file-allocation=none --allow-overwrite=true --auto-file-renaming=false \
      --console-log-level=notice -d "${DIAMOND_DB_DIR}" -o nr.gz "${NCBI_FASTA_URL}"
  elif command -v curl >/dev/null 2>&1; then
    curl -fL --retry 8 --retry-delay 3 --connect-timeout 30 -o "${out}" "${NCBI_FASTA_URL}"
  else
    echo "[build_diamond_nr] 错误: 需要 aria2c 或 curl" >&2
    return 1
  fi
  echo "${out}"
}

FASTA=""
if FASTA="$(resolve_fasta)"; then
  echo "[build_diamond_nr] 使用已有 FASTA: ${FASTA}"
  makedb_from_fasta "${FASTA}"
elif makedb_from_blast_nr; then
  :
else
  echo "[build_diamond_nr] 本地无 FASTA / BLAST 管道失败，改为下载 NCBI nr.gz"
  FASTA="$(download_nr_gz)"
  makedb_from_fasta "${FASTA}"
fi

if [[ ! -f "${PREFIX}.dmnd" ]]; then
  echo "[build_diamond_nr] 错误: 未生成 ${PREFIX}.dmnd" >&2
  exit 1
fi
ls -lah "${PREFIX}.dmnd"
echo "[build_diamond_nr] 完成。config 中设置 blast.database: /db/ncbi/diamond/nr 与 blast.engine: diamond"
