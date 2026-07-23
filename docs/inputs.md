# 输入准备

v0.1 将一个 Space Ranger `outs` 目录视为一个空间切片或建库样本。仓库不包含真实数据
或测试数据，用户需要在仓库外准备输入，并在样本表中引用。

## 最低输入

运行默认 `qc` 模块只需要：

1. 一个或多个可读取的 Space Ranger `outs` 目录；
2. `config/config.yaml`；
3. `config/samples.tsv`。

Space Ranger 目录通常应包含 filtered feature-barcode matrix、空间坐标和
scalefactors。raw matrix 和已配准组织图像不是构建标准化表达矩阵的硬性要求，但缺少
它们时，相应 background QC 或 H&E 配准证据会显示为不可用。

v0.1 不接受只有表达矩阵和 `array_row/array_col` 的通用输入。后续版本会通过独立
adapter 支持这类数据。

## 创建配置文件

从模板开始：

```bash
cp config/config.template.yaml config/config.yaml
cp config/samples.template.tsv config/samples.tsv
cp config/qc_reviews.template.tsv config/qc_reviews.tsv
```

`config/defaults.yaml` 保存仓库默认值，不建议为单个项目直接修改它。项目差异写入
`config/config.yaml`；未覆盖的内容会继续使用默认值。

## 样本表

`samples.tsv` 是制表符分隔文件，一行代表一个空间切片。最低必填列是：

| 列 | 含义 |
|---|---|
| `sample_id` | 脱敏且唯一的样本编号，也会用于结果目录名 |
| `input_type` | v0.1 固定填写 `spaceranger` |
| `input_path` | 对应 Space Ranger `outs` 目录 |

优先使用相对于 `samples.tsv` 所在目录的路径。例如：

```text
sample_id	input_type	input_path
sample_01	spaceranger	../data/spaceranger/sample_01/outs
sample_02	spaceranger	../data/spaceranger/sample_02/outs
```

`sample_id` 只能使用字母、数字、点、下划线和连字符，并应以字母或数字开头。不要使用
姓名、病历号、送样编号或能直接追溯到个体的标识。

## 推荐的样本信息

基础 QC 不要求实验设计字段，但为后续模块保留以下信息会更方便：

- `animal_id` 或其他脱敏后的独立生物学个体编号；
- `biological_replicate`；
- `technical_batch`、`slide_id`、`capture_area` 和 `library_id`；
- `assay`、`species`、`genome_reference` 和 probe set 信息；
- `section_level` 和 `orientation`；
- `genotype`、`treatment` 和 `condition`。

spot 或 ROI 不能代替生物学重复。`condition_2x2` 会根据独立样本数选择分析路径。

## H&E 图像

如果 Space Ranger 输出中存在与空间坐标匹配的 registered image 和 scalefactors，
`qc` 会生成 H&E/spot overlay。流程不会：

- 把完整 H&E 图像嵌入大型 AnnData；
- 猜测图像与 scalefactor 的配对；
- 自动修正旋转、平移或形变。

缺少图像时，表达量和坐标 QC 仍可运行，配准分项会标记为缺少证据。

## 综合 QC 评分输入

仓库默认使用 `config/qc_profiles/unconfigured_v1.yaml`。其中数值阈值故意为空，因此
第一次运行会保留 `UNCALIBRATED`，不会把未知质量写成通过。正式评分前：

1. 复制该文件并创建新的、带版本号的 assay profile；
2. 根据明确的实验平台、组织类型和实验 SOP 填写 counts、detected genes 和
   mitochondrial fraction 阈值，并在 profile 的 `assays` 中列出兼容的 assay；
3. 在 `config/config.yaml` 中把 `qc.score.profile` 指向新 profile；
4. 查看每个样本的 H&E overlay 和 spatial QC 图；
5. 从 `config/qc_reviews.template.tsv` 创建审核表，并为
   `image_alignment`、`spatial_artifacts` 记录决定和证据。

完成的 PASS/WARN/FAIL 审核必须同时填写证据路径、脱敏 reviewer ID 和 ISO-8601
审核时间。若对应图实际不可用或检查被关闭，审核表不能把该分项覆盖为 PASS。
已校准 profile 的 `assays` 必须与 `samples.tsv` 的 `assay` 精确匹配，避免把另一平台的
阈值误用于当前数据。

固定的 v1 权重与 PASS/WARN/FAIL 分值属于评分方法版本，不能在同一版本中任意修改。
需要改变方法时应升级 `method_version`、schema、文档和测试。综合分只用于复核与报告，
不会自动过滤 spot 或样本。

## ROI 模块输入

启用 `roi` 或 `svg` 时，在对应样本行填写：

| 列 | 含义 |
|---|---|
| `roi_path` | ROI 导出表 |
| `roi_barcode_column` | barcode 列名 |
| `roi_label_column` | ROI 标签列名 |

ROI 表至少应提供能与该样本 spot 对应的 barcode 以及一个区域标签。流程只使用明确的
barcode 连接规则，不做模糊字符串匹配。

当 `roi`、`svg`、`condition_2x2` 或 `pathway` 由配置启用时，当前版本要求每个样本都有
`roi_path`；未启用这些模块时，该列可以省略或留空。

如果不同项目中的 ROI 名称需要合并，复制别名模板：

```bash
cp config/roi_label_aliases.template.tsv config/roi_label_aliases.tsv
```

再记录经过确认的 source label 与 canonical label。未经确认的别名不应当作解剖学真值。

## 2×2 条件比较输入

启用 `condition_2x2` 时，每个样本必须填写两个因素及其独立样本信息。配置中还要明确：

- factor A 的 reference 和 alternative 水平；
- factor B 的 reference 和 alternative 水平；
- 分析模式，通常使用 `auto`。

`auto` 的行为是：

- 每个设计 cell 一个独立样本：生成描述性 effect，不生成样本层面的 p 值或 FDR；
- 每个设计 cell 有足够的独立生物学重复：使用重复样本分析；
- 设计不完整或 replicate 身份不清：先报告 eligibility，再停止条件比较。

两条分支都会逐 ROI 检查每个独立 section/unit 的 spot 数；不足
`analysis.condition.min_roi_spots_per_unit` 的 ROI 不进入 effect ranking 或模型。

## Pathway 输入

启用 `pathway` 时，复制 gene-set 清单：

```bash
cp config/pathway_gene_sets.template.tsv config/pathway_gene_sets.tsv
```

为每个启用的 GMT 资源填写相对路径、真实 SHA-256、来源和版本说明。模板中的占位
checksum 不能用于正式分析。v0.1 只把描述性 2×2 effect 送入 pathway；replicated
PyDESeq2 结果不会自动转换为 pathway ranking。

## 输入检查

编辑完成后先运行：

```bash
snakemake \
  --snakefile workflow/Snakefile \
  --directory . \
  --cores 1 \
  --dry-run
```

如果输入路径、schema 或模块依赖有问题，应先修正这些问题，再启动长时间任务。
