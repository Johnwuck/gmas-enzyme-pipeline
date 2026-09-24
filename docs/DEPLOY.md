# 远程部署指南

关键词 → UniProtKB → NCBI BLAST → AlphaFold DB → 交集去重 → DLKcat 的整合包部署说明。

## 1. 机器要求

| 项 | 建议 |
|----|------|
| 系统 | Linux x86_64（含 WSL2） |
| 磁盘 | ≥ 20 GB 可用（镜像 + Swiss-Prot + 结果/结构） |
| 内存 | ≥ 8 GB |
| 网络 | 可访问 UniProt、NCBI、AlphaFold DB、PubChem、GitHub |
| GPU | **不需要**（默认 CPU PyTorch）；有 GPU 时可自行改装 CUDA 版 pytorch |

## 2. 方式 A：Docker（推荐）

### 2.1 安装 Docker / Compose

按官方文档安装 Docker Engine 与 Docker Compose 插件。

### 2.2 获取代码

将本仓库拷到服务器（git clone / scp / rsync），进入项目根目录。

### 2.3 构建镜像

```bash
docker compose build
```

构建期会：创建 conda 环境 `gmas_pipeline`、安装 blast/cd-hit/foldseek/pytorch/rdkit，并 `vendor` DLKcat。

### 2.4 BLAST / DIAMOND 库（本地比对）

本机推荐复用共享库目录（默认）：

`/mnt/sdb/tmp/Databases/NCBI` → 容器内 `/db/ncbi`

```bash
# 下载 / 更新 Swiss-Prot（blastp 小库，约 215 MB）
docker compose --profile setup run --rm download_blastdb

# 构建 DIAMOND NR（推荐搜全量 NR；可从已有 BLAST NR 导出，耗时长）
docker compose --profile setup run --rm build_diamond_nr
# 或: NCBI_DB_ROOT=/mnt/sdb/tmp/Databases/NCBI bash scripts/build_diamond_nr.sh

# Foldseek-AFDB（序列搜 AFDB；默认 UniProt50-minimal）
docker compose --profile setup run --rm download_foldseek_afdb
# 可选更小库: FOLDSEEK_AFDB_NAME=Alphafold/Swiss-Prot docker compose --profile setup run --rm download_foldseek_afdb

# 增量更新 NR / NT BLAST 分卷（需 conda 环境 ncbi + awscli；TB 级，耗时长）
conda activate ncbi
bash scripts/update_ncbi_nr_nt.sh          # 默认 SOURCE=aws
# bash scripts/update_ncbi_nr_nt.sh nr     # 仅 NR
# SOURCE=ncbi NUM_THREADS=2 bash scripts/update_ncbi_nr_nt.sh nt
```

日志：`/tmp/update_nr_aws.log`、`/tmp/update_nt_aws.log`（或脚本生成的带时间戳日志）。

`config.yaml` 中推荐：

- `blast.engine: diamond` + `blast.database: /db/ncbi/diamond/nr`（NR，快）
- 或 `blast.engine: blastp` + `blast.database: /db/ncbi/blastdb/swissprot`（小库）

亦可 `blastp` + `/db/ncbi/NR/nr`（极慢，不推荐大批量 query）。裸机请改成宿主机绝对路径。  
若暂不下库，可将 `blast.mode` 改为 `remote`（慢、受 NCBI 限流；DIAMOND 不支持 remote）。

### 2.5 运行流水线

```bash
# 全流程
docker compose run --rm pipeline --keyword gmas --step all

# 仅 DLKcat（需已有 results/gmas/final/complete_enzymes.fasta）
docker compose run --rm pipeline --keyword gmas --step dlkcat

# 查看帮助
docker compose run --rm pipeline
```

数据与结果通过卷挂载持久化在宿主机 `./data`、`./results`。

## 3. 方式 B：裸机 Conda

### 3.1 安装 Miniconda

略。确保 `conda` 在 PATH 中。

### 3.2 一键环境

```bash
bash scripts/setup_env.sh
# 仅环境 + DLKcat，稍后下库:
# bash scripts/setup_env.sh --skip-blastdb
# bash scripts/download_blastdb.sh
```

`setup_env.sh` 会：

1. 按 `environment.yml` 创建/更新 **`gmas_pipeline`**
2. 执行 `scripts/vendor_dlkcat.sh` → `third_party/DLKcat`
3. （默认）下载 Swiss-Prot BLAST 库

### 3.3 运行

```bash
conda activate gmas_pipeline
cd /path/to/gmas_enzyme_pipeline
cp config.example.yaml config.yaml   # 首次
python -m src.pipeline --keyword gmas --config config.yaml
```

## 4. 关键配置项

模板见 [`config.example.yaml`](../config.example.yaml)；本地复制为 `config.yaml` 后修改（`config.yaml` 不入库）。

| 键 | 含义 |
|----|------|
| `keyword` | 搜索关键词；CLI `--keyword` 优先 |
| `query` | 空则自动 `(gene:{kw}) OR (protein_name:{kw}) OR {kw}` |
| `substrates` | DLKcat 底物（名称 + SMILES） |
| `blast.mode` | `local` / `remote` |
| `blast.engine` | `diamond`（推荐 NR）或 `blastp` |
| `blast.database` | 如 `/db/ncbi/diamond/nr` 或 `/db/ncbi/blastdb/swissprot` |
| `dlkcat.root` | 默认 `third_party/DLKcat/DeeplearningApproach` |
| `dlkcat.python` | 默认 `python`（当前环境解释器） |

**禁止**在远程配置中写死某台机器的绝对路径（如 `/home/xxx/...`）。

## 5. 环境依赖一览

定义文件：[`environment.yml`](../environment.yml)

| 组件 | 来源 |
|------|------|
| Python 3.12、requests、biopython、pandas、pyyaml | conda-forge |
| blast、diamond、cd-hit、foldseek | bioconda |
| pytorch（CPU）+ rdkit + scikit-learn | pytorch / conda-forge |
| DLKcat 代码与权重 | `scripts/vendor_dlkcat.sh` → GitHub SysBioChalmers/DLKcat |
| Swiss-Prot BLAST DB | `scripts/download_blastdb.sh` |
| DIAMOND NR（.dmnd） | `scripts/build_diamond_nr.sh` |
| Foldseek AFDB | `scripts/download_foldseek_afdb.sh` |

可选 GPU：自行在环境中改用 `pytorch-cuda` 渠道包，并调整 Docker base；本仓库默认 CPU 以保证无 GPU 远程机可构建。

## 6. 脚本速查

| 脚本 | 作用 |
|------|------|
| `scripts/setup_env.sh` | 裸机一键环境 |
| `scripts/vendor_dlkcat.sh` | 拉取/解压 DLKcat |
| `scripts/download_blastdb.sh` | 下载 Swiss-Prot |
| `scripts/build_diamond_nr.sh` | 构建 DIAMOND NR（.dmnd） |
| `scripts/download_foldseek_afdb.sh` | 下载 Foldseek 格式 AFDB |
| `scripts/entrypoint.sh` | 容器入口 |

```bash
bash scripts/setup_env.sh --help
FORCE_VENDOR=1 bash scripts/vendor_dlkcat.sh
# 本机已有 DLKcat 时可:
# LOCAL_DLKCAT=/path/to/DLKcat VENDOR_MODE=link bash scripts/vendor_dlkcat.sh
FORCE_BLASTDB=1 bash scripts/download_blastdb.sh
```

## 7. 常见问题

**UniProt / AFDB 429 或超时**  
网络不稳定时重试；流水线已有指数退避。可单步重跑：`--step afdb` / `--step finalize`。

**最终酶集少于 BLAST 命中**  
设计如此：最终只保留 **NCBI 序列 ∩ AlphaFold DB 结构** 再去重。无结构的 accession 见 `dropped_no_afdb.tsv`。

**DLKcat：根目录不存在**  
运行 `bash scripts/vendor_dlkcat.sh`，或重建 Docker 镜像。

**DLKcat：缺少 fingerprint_dict.pickle**  
确认 `third_party/DLKcat/DeeplearningApproach/Data/input.zip` 已解压（vendor 脚本会处理）。

**blastp 找不到库**  
检查 `/db/ncbi/blastdb/swissprot.pin`（宿主机对应 `.../NCBI/blastdb/swissprot.pin`）；或改 `blast.mode: remote`。

**容器内改了 config 不生效**  
`config.yaml` 以只读挂载；请改宿主机文件后重跑。

## 8. 验收清单

- [ ] 已从 `config.example.yaml` 复制并改好 `config.yaml`（无误写的他人机器路径）  
- [ ] `docker compose build` 成功，或 `setup_env.sh` 成功  
- [ ] local 模式下 DIAMOND/BLAST 库路径可读（如 `/db/ncbi/diamond/nr.dmnd` 或 Swiss-Prot）  
- [ ] `third_party/DLKcat/DeeplearningApproach` 存在（`bash scripts/vendor_dlkcat.sh`）  
- [ ] `python -m src.pipeline --keyword gmas --step dlkcat`（或 Docker 等价命令）可产出预测表  

使用细节见 [USAGE.md](USAGE.md)。
