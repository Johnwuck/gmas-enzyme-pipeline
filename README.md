# GMAS Enzyme Pipeline

关键词驱动的酶同源搜索与 DLKcat kcat 预测流水线。

```
关键词
  → UniProtKB 种子
  → CD-HIT 去重
  → DIAMOND/BLAST（NR）扩同源
  → AlphaFold DB 结构
  → Foldseek / 序列映 AFDB
  → NCBI∩AFDB 交集再去重
  → DLKcat kcat 预测
```

> 本流水线使用 **AlphaFold DB** 公开模型 + Foldseek，不跑本地 AlphaFold3。

## 流水线步骤

| 步骤 | 作用 |
|------|------|
| `uniprot` | 按关键词搜 UniProtKB，下载种子序列 |
| `seed_dedupe` | 种子 FASTA 做 CD-HIT，缩小后续比对量 |
| `blast1` / `blast2` | 本地 DIAMOND（NR）或 blastp 扩同源 |
| `afdb` | 按 UniProt accession 下载 AlphaFold DB PDB |
| `foldseek_afdb` | 序列映射 AFDB（DIAMOND→Swiss-Prot 或 Foldseek-AFDB） |
| `foldseek1` / `foldseek2` | 本地结构 all-vs-all |
| `dedupe` | 中间 CD-HIT 去重 |
| `finalize` | NCBI∩AFDB 交集后再去重，写出最终酶集 |
| `dlkcat` | 最终酶集 × 配置底物 → kcat 预测 |

逐步说明与产物路径见 [docs/USAGE.md](docs/USAGE.md)。

## 快速开始

### 1. 配置

```bash
cp config.example.yaml config.yaml
# 按需改 keyword、substrates、数据库路径
```

### 2. Docker（推荐）

```bash
docker compose build
docker compose --profile setup run --rm build_diamond_nr   # 一次性构建 DIAMOND NR
# 或小库: docker compose --profile setup run --rm download_blastdb
docker compose run --rm pipeline --keyword Phytase --step all
```

### 3. 裸机 Conda

```bash
bash scripts/setup_env.sh
conda activate gmas_pipeline
python -m src.pipeline --keyword Phytase --config config.yaml --step all
```

## 文档

| 文档 | 内容 |
|------|------|
| [docs/USAGE.md](docs/USAGE.md) | 步骤、命令、产物路径、换底物/关键词 |
| [docs/DEPLOY.md](docs/DEPLOY.md) | Docker / 裸机部署、库下载、FAQ |
| [config.example.yaml](config.example.yaml) | 远程部署配置模板 |

## 仓库结构

```
├── src/                  # 流水线源码
│   ├── pipeline.py       # CLI 与步骤编排
│   ├── uniprot_client.py
│   ├── blast_runner.py
│   ├── alphafold_db.py
│   ├── foldseek_runner.py
│   ├── dedupe.py
│   ├── dlkcat_runner.py
│   └── ...
├── scripts/              # 环境、库下载、容器入口
├── docs/                 # 使用与部署文档
├── data/substrates/      # 底物 SMILES 缓存（示例）
├── config.example.yaml
├── docker-compose.yml
├── Dockerfile
├── environment.yml
└── requirements.txt
```

`third_party/DLKcat`、`results/`、大体量结构与 BLAST/DIAMOND 库不入库，部署时按文档拉取/挂载。

## 最终酶集含义

**最终酶 = 有序列（NCBI/UniProt）且有 AlphaFold DB 结构，再经 CD-HIT 去重。**  
DLKcat 只对该集合预测；预测值为模型估计，需实验验证。
