# 使用指南

## 流水线步骤

```
uniprot → seed_dedupe → blast1 → afdb → foldseek_afdb → foldseek1 → dedupe → blast2 → foldseek2 → finalize → dlkcat
```

| 步骤 | 作用 |
|------|------|
| `uniprot` | 按关键词搜 UniProtKB 并下载序列 |
| `seed_dedupe` | 对种子 FASTA 做 CD-HIT，缩小后续比对查询量 |
| `blast1` / `blast2` | 本地 DIAMOND（NR）或 blastp 扩同源 |
| `afdb` | 按 UniProt accession 直接下 AlphaFold DB 结构 |
| `foldseek_afdb` | 用序列映射 AFDB 结构（优先 Foldseek-AFDB；否则 DIAMOND→Swiss-Prot→下 PDB） |
| `foldseek1` / `foldseek2` | 本地结构 all-vs-all |
| `dedupe` | 中间 CD-HIT 去重 |
| `finalize` | NCBI∩AFDB 交集后再去重，写最终酶集与 PDB |
| `dlkcat` | 对最终酶集 × 配置底物预测 kcat |

## 常用命令

裸机：

```bash
conda activate gmas_pipeline
python -m src.pipeline --keyword gmas --config config.yaml
python -m src.pipeline --keyword gmas --step uniprot,finalize,dlkcat
python -m src.pipeline --keyword gmas --step dlkcat --force
```

Docker：

```bash
docker compose run --rm pipeline --keyword gmas --step all
docker compose run --rm pipeline --keyword gmas --step dlkcat
```

`--keyword` 优先于配置文件中的 `keyword`；输出目录按关键词隔离。

## 产物路径

| 路径 | 说明 |
|------|------|
| `data/{keyword}/raw/` | UniProt 种子 FASTA/TSV |
| `data/{keyword}/structures/` | AFDB PDB |
| `results/{keyword}/round1|round2/` | BLAST / Foldseek / foldseek_afdb |
| `results/{keyword}/round1/foldseek_afdb.tsv` | 序列搜 AFDB 命中 |
| `results/{keyword}/final/complete_enzymes.fasta` | 最终酶集 |
| `results/{keyword}/final/complete_summary.tsv` | 汇总（均有 NCBI+AFDB） |
| `results/{keyword}/final/structures/` | 最终 PDB |
| `results/{keyword}/final/dlkcat_predictions.tsv` | kcat 预测 |
| `results/{keyword}/final/dropped_no_afdb.tsv` | 无 AFDB 的丢弃项 |
| `data/substrates/` | 底物缓存 JSON |

## 更换底物

编辑本地 `config.yaml`（由 `config.example.yaml` 复制）：

```yaml
substrates:
  - name: L-glutamate
    smiles: "C(CC(=O)O)C(C(=O)O)N"
  - name: methylamine
    smiles: "CN"
```

仅名称、无 SMILES 时会尝试 PubChem 补全。改完后：

```bash
python -m src.pipeline --keyword gmas --step dlkcat --force
```

## 更换关键词

```bash
python -m src.pipeline --keyword your_gene --step all
```

新关键词使用独立目录 `data/your_gene/`、`results/your_gene/`，不覆盖已有 `gmas` 结果。

## 最终集含义

**最终酶 = 有序列（NCBI/UniProt）且有 AlphaFold DB 结构，再经 CD-HIT 去重。**  
DLKcat 只对该集合预测；预测值为模型估计，需实验验证。

部署与环境见 [DEPLOY.md](DEPLOY.md)。
