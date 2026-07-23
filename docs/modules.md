# 模块清单

流程由小型模块组成。推荐在 `config/config.yaml` 中统一选择模块，让一次运行的选择可以
被报告和复现。

```yaml
modules:
  enabled:
    - qc
    - core
    - report
  auto_dependencies: true
```

默认只启用 `qc`。`auto_dependencies: true` 会补齐上游模块，但不会替用户创造缺失的
ROI、实验设计或 gene-set 资源。

## 稳定模块

| 模块/target | 做什么 | 依赖 | 典型输出 |
|---|---|---|---|
| `qc` | 输入检查、标准化、六项 QC 和综合评分 | Space Ranger 输入 | `results/input/`、`results/qc/`、`work/ingested/` |
| `core` | eligibility、HVG/PCA、表达图、UMAP/Leiden、空间 domain | `qc` | `results/preprocessing/`、`results/embeddings/`、`results/spatial/` |
| `roi` | ROI coverage、raw-count pseudobulk、ROI-vs-rest effect | `qc`、ROI 表 | `results/roi/` |
| `svg` | 每样本 ROI 图和空间变异基因候选 | `qc`、ROI 表 | `results/svg/<sample>/` |
| `condition_2x2` | 两因素 main/simple/interaction effect | `roi`、2×2 样本信息 | `results/condition/` |
| `pathway` | 对描述性 condition ranking 做可恢复的 prerank 富集 | 描述性 `condition_2x2`、GMT 清单 | `results/pathway/` |
| `figures` | 为已完成模块生成通用图及 source tables | 相应分析模块 | `results/figures/` |
| `resource_report` | 汇总资源监控日志 | 资源日志 | `results/reporting/` |
| `report` | 汇总 QC、图、模块状态、运行来源和文件索引 | 已启用模块 | `results/report/` |

`full` 是一个方便目标，会请求自包含的稳定分析模块。它要求 ROI 和 2×2 设计已配置；
`pathway` 依赖外部 GMT，`resource_report` 依赖预先生成的监控日志，因此两者必须显式
加入，不会被 `full` 暗中启用。

正式使用时在 `modules.enabled` 中写 `full` 并运行默认入口。named `full` 会检查配置是否
已选中同一组模块，不允许绕过配置，否则 run report 无法准确描述本次 DAG。

## 数据流与依赖

`roi` 和 `svg` 都从标准化输入、QC/eligibility 与 ROI 标注开始，二者彼此独立。
因此可以只运行 ROI 汇总，也可以只运行每样本 SVG。

`condition_2x2` 依赖 ROI pseudobulk。其内部先审计四个设计 cell：

- 每 cell 一个独立样本时，输出描述性 effect 和方向；
- 有足够的独立生物学重复时，进入 replicated 分支；
- 不完整或身份不清时，保留审计结果并停止推断。

样本数审计由 `workflow/scripts/metadata/audit_sample_design.py` 单独完成，同时输出每个
cell 的样本行数和唯一生物学单位数。当前版本不把同一动物的多张切片当作独立重复，也
不支持 paired/repeated-measure 模型；这两种情况会明确标记或停止。

v0.1 的 `pathway` 只接受描述性分支的 effect ranking。replicated 分支包含
p 值和 FDR，必须先确定并记录适合该推断结果的排名策略；当前流程会明确停止，而不会
把两种结果静默混用。

## 综合 QC

QC 总览由六个分项组成：

1. `in_tissue` 数据完整性和一致性；
2. 每个 spot 的 total counts；
3. 每个 spot 的 detected genes；
4. mitochondrial fraction；
5. H&E 与 spot 配准；
6. 空间伪影。

报告同时展示 0–100 分数和 evidence coverage。正式 `qc_score` 只有在六项证据完整时
才出现；达到最低覆盖度但仍缺项时，只显示明确命名的 `provisional_score`。缺失、关闭、
无法计算或尚未审核的分项不会被当作零分。综合评分用于快速发现需要复核的样本，不会
自动删除 spot 或样本。

## 使用 named target

配置文件适合正式运行；named target 适合单独检查一个模块：

```bash
snakemake \
  --snakefile workflow/Snakefile \
  --directory . \
  --cores 1 \
  --dry-run \
  roi
```

正式执行时移除 `--dry-run`，并按需加入 `--sdm conda`。

查看当前版本可点名的目标：

```bash
snakemake \
  --snakefile workflow/Snakefile \
  --directory . \
  --list-target-rules
```

## Figures 与 report

`figures` 只绘制已有分析结果，不重新定义分析 population。图表应同时保留 source table，
方便追踪和重画。

`report` 汇总当前选择，不会为未启用模块编造结果。模块状态应区分：

- completed；
- not_requested；
- review_required；
- completed_with_qc_flags；
- completed_no_eligible_results；
- completed_with_model_failures。

若某条 Snakemake rule 失败，本次 DAG 会直接以非零状态停止，不会生成一份把失败模块
写成 completed 的最终报告。具体失败原因保留在相应日志中。

H5AD、完整矩阵和原始图像只出现在 artifact manifest 中，不嵌入 HTML。
run manifest 分别记录 defaults、活动 override、samples 和合并后 effective config 的
校验值；effective config 会把仓库内路径转为相对路径，并脱敏外部绝对路径。

## 实验性外部比较器

仓库保留一个旧格式 external validation comparator，供兼容特定历史目录结构。它具有以下
边界：

- experimental、specialized，不属于通用 v0.1 主流程；
- 不由 `full` 自动运行；
- 需要用户明确提供匹配格式的外部目录；
- 外部聚类或报告只是 comparator，不是真值；
- 在公开使用前应改造成带显式 adapter schema 的通用模块。
