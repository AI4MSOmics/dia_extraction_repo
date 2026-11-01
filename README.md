# DIA Spectrum Evidence Extraction

[![License: Apache-2.0](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](LICENSE)
[![MSHub Plugin](https://img.shields.io/badge/MSHub-ms.s1.dia--extraction-teal)](plugin.yaml)

---

## English

**MSHub capability pack** for DIA spectrum evidence extraction.
Part of the [AI4MSOmics](https://github.com/AI4MSOmics) ecosystem.
See [MSHub](https://github.com/AI4MSOmics/MSHub) for the hub specification.

Extracts observed MS2 spectra from DIA raw files using RT-anchored targeted extraction, and renders per-peptide A4 evidence pages with isoPSV-style layout: info lines, ion ladder, observed spectrum, inverted AlphaPeptDeep prediction, ppm panel, and RT heatmap.

### Installation

```bash
git clone https://github.com/AI4MSOmics/dia_extraction_repo.git
cd dia_extraction_repo
# Dependencies: Python >=3.8 (numpy/pandas/matplotlib/pymupdf/openpyxl/h5py),
#               R >=4.0 (base), ThermoFisher RawFileReader DLLs, AlphaPeptDeep
```

### Usage — three integration paths

#### 1. MSHub skill card

This repository follows the [MSHub plugin specification](https://github.com/AI4MSOmics/MSHub/blob/main/docs/03-plugin-spec.md). The [SKILL.md](SKILL.md) at the root is the MSHub-format skill card.

#### 2. ZCode / Claude / DeepSeek skill

Copy `zcode/` to your agent's skills directory, or configure `agent/tool_manifest.json` for function calling (see [agent/deepseek_guide.md](agent/deepseek_guide.md)).

#### 3. Standalone CLI

```bash
# Full chain (7 scripts)
python scripts/70_build_s14_targets.py       # Parse search results → targets
python scripts/71_extract_targeted.py        # RT-anchored extraction
python scripts/73_profiles.py                # Full-RT XIC profiles
python scripts/76_build_isopsv_inputs.py     # Build isoPSV R inputs
cd risoPSV && env -u LC_CTYPE Rscript 78_s14_export.R  # Render top blocks
python scripts/79_compose_a4.py              # Compose final A4 PDFs
python scripts/85_audit.py                   # Automated audit
```

### Key features

| Feature | Method | Verified on |
|---------|--------|-------------|
| RT-anchored extraction | Sequence-ladder apex scoring in ±0.6 min window | 3,319/3,319 rows, median ΔRT 0.035 min |
| isoPSV evidence layout | R engine verbatim (ion ladder, spectrum, ppm) | 112-peptide PDF, 3-round visual acceptance |
| Inverted AlphaPeptDeep | Mirror prediction panel + cosine annotation | 112/112 coverage (11 on-the-fly) |
| RT heatmap | Spectrum-correlation XIC + 3-tier background removal | 19→0 wide-coverage peptides |
| Table 1 cross-reference | Dual-chain (S14 + v9-master fallback) | 159/159 peptides, 0 placeholders |
| Automated audit | 3-layer (data/T1/presentation) programmatic check | errors=0, warnings=1 |

### Documentation

| Document | Description |
|----------|-------------|
| [SKILL.md](SKILL.md) | MSHub skill card (agent-facing) |
| [docs/architecture.md](docs/architecture.md) | Script chain and data flow |
| [docs/configuration.md](docs/configuration.md) | Env vars and parameters |
| [docs/troubleshooting.md](docs/troubleshooting.md) | 15+ known issues |
| [plugin.yaml](plugin.yaml) | MSHub capability pack manifest |
| [schemas/](schemas/) | Input/output JSON Schemas |

### Related projects

- [MSHub](https://github.com/AI4MSOmics/MSHub) — Agent-native capability hub
- [eval_ms_pipeline](https://github.com/AI4MSOmics/eval_ms_pipeline_repo) — DDA parameter evaluation (companion)
- [isoPSV](https://github.com/Songhao-Jiang/isoPSV) — Spectrum visualization engine
- [AlphaPeptDeep](https://github.com/AI4MSOmics/AlphaPeptDeep) — Theoretical spectrum prediction

---

## 中文说明

**MSHub 能力包** — DIA 谱图证据提取流水线。
隶属于 [AI4MSOmics](https://github.com/AI4MSOmics) 生态。

从 DIA raw 文件中用 RT 锚定靶向提取实测 MS2 谱，以 isoPSV 版式渲染逐肽 A4 证据页：信息行、离子梯、观测谱图、倒置 AlphaPeptDeep 理论谱、质量差面板、RT 热图。

### 三种使用方式

#### 1. MSHub 能力包
根目录 [SKILL.md](SKILL.md) 为 MSHub 格式技能卡。

#### 2. ZCode / DeepSeek 技能
```bash
cp -r zcode/ ~/.zcode/skills/dia_extraction/
```

#### 3. 命令行直接运行
```bash
python scripts/70_build_s14_targets.py
python scripts/71_extract_targeted.py
python scripts/73_profiles.py
python scripts/76_build_isopsv_inputs.py
cd risoPSV && env -u LC_CTYPE Rscript 78_s14_export.R
python scripts/79_compose_a4.py
python scripts/85_audit.py
```

### 核心功能

| 功能 | 方法 | 验证结果 |
|------|------|---------|
| RT 锚定提取 | 序列梯顶点打分（±0.6 min 窗内） | 3,319/3,319 全部成功，ΔRT 中位 0.035 min |
| isoPSV 版式 | R 引擎原样渲染（离子梯/谱图/ppm） | 112 肽 PDF，三轮视觉验收 |
| 倒置理论谱 | AlphaPeptDeep 镜像 + cosine 标注 | 112/112 覆盖（11 条按需补算） |
| RT 热图 | 谱型相关 XIC + 三层背景去除 | 大面积显色肽 19→0 |
| Table 1 对应 | 双链（S14 + v9 主表兜底） | 159/159 全有谱、零占位 |
| 自动审计 | 三层程序化检查 | 0 错误、1 警告 |

---

## License / 许可证

Apache-2.0. See [LICENSE](LICENSE).
