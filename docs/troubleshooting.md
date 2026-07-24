# 故障排查

遇到问题时先运行 dry-run。它通常能在长任务开始前发现配置、输入和依赖错误。

```bash
snakemake \
  --snakefile workflow/Snakefile \
  --directory . \
  --cores 1 \
  --dry-run
```

## 找不到 config 文件

确认已经从模板创建活动配置和样本表：

```bash
cp config/config.template.yaml config/config.yaml
cp config/samples.template.tsv config/samples.tsv
cp config/qc_reviews.template.tsv config/qc_reviews.tsv
```

不要把包含真实路径或样本信息的活动文件提交到仓库。

## 样本路径不存在

`input_path` 和 `roi_path` 的相对路径从 `samples.tsv` 所在目录解析，而不是从终端当前目录
解析。优先使用相对路径，并确认它们指向实际文件。

v0.1 的 `input_path` 应指向 Space Ranger `outs` 目录，而不是任意 H5AD 或表达矩阵。

## Config 或 samples schema 失败

常见原因包括：

- 列名拼写错误；
- YAML 缩进错误；
- 使用了当前版本不支持的配置项；
- `sample_id` 重复或含空格；
- `input_type` 不是 `spaceranger`；
- 数值参数超出允许范围。

根据报错中的字段名修正模板副本，不要通过关闭 schema 绕过错误。

## 线粒体比例显示不可用

流程需要能够识别 gene symbol 并匹配配置中的线粒体前缀。如果输入只有无法映射的
feature ID，该分项会保持 NA。它不是零分，也不会由流程猜测。

检查 feature metadata、物种和 `qc.mitochondrial` 配置。

## QC 总分显示 UNCALIBRATED 或 PENDING

这是安全状态，不是程序失败。`UNCALIBRATED` 表示当前 profile 没有经过确认的数值阈值；
`PENDING` 表示配准或空间伪影还没有人工审核。按
[输入准备](inputs.md#综合-qc-评分输入)完成 profile 和 review 表，不要为了得到总分而
随意填写阈值或 PASS。

## 没有 H&E overlay

确认 Space Ranger 目录中存在 registered image、spatial coordinates 和相匹配的
scalefactors。流程不会尝试错误配对，也不会自动进行图像变换。

图像缺失时，其他 QC 可以继续，配准分项会降低 evidence coverage。

## ROI barcode 对不上

依次检查：

1. ROI 表是否属于当前样本；
2. `roi_barcode_column` 和 `roi_label_column` 是否正确；
3. barcode 是否保留或移除了 10x suffix；
4. 是否存在重复 barcode；
5. ROI 标签是否需要经过明确审核的 alias。

不要用模糊匹配强行连接 barcode。匹配失败应保留为可审计问题。

## 2×2 模块提示 design not eligible

四个 factor 组合必须完整，并且 biological replicate 身份必须可信。

- 每 cell 一个独立样本只能进入 descriptive 分支；
- replicated 分支需要每 cell 有足够的独立样本；
- spot 数、ROI 数或同一样本中的切片行不能冒充生物学重复。

检查样本表中的 factor 水平、`animal_id` 和 `biological_replicate`，再重新 dry-run。

## Pathway 模块拒绝运行

检查：

- `condition_2x2` 是否成功产生 ranking；
- gene-set manifest 是否存在；
- 当前 condition 是否为描述性分支；v0.1 不会把 replicated 结果自动排名；
- 需要使用的行是否设为 enabled；
- GMT 路径是否可读取；
- SHA-256 是否为真实校验值而非模板占位符。

资源来源和版本不明确时，应先冻结资源，而不是跳过校验。

## 内存不足或任务太慢

先降低 `--cores`，因为并发任务会叠加内存。然后查看 `logs/` 中对应 rule 的日志和
benchmark。不要仅因机器核数较多就把 `--cores` 设为全部核数。

大型中间矩阵位于 `work/`。确认可用磁盘空间，并避免在同一目录中保留多个未使用的历史
run。

## 运行被中断

再次执行原命令即可。Snakemake 会继续缺失或过期的任务。

如果上一次任务留下了 Snakemake 锁，先确认没有其他 Snakemake 进程正在使用同一工作
目录，再执行：

```bash
snakemake \
  --snakefile workflow/Snakefile \
  --directory . \
  --unlock
```

不要删除已经完成的大文件来“解除锁定”。

## 输出存在，但 Snakemake 想重跑

常见原因是 config、脚本、输入文件时间或 rule 定义发生变化。先查看文件状态摘要：

```bash
snakemake \
  --snakefile workflow/Snakefile \
  --directory . \
  --summary
```

再运行一次 dry-run；Snakemake 默认会在计划中说明任务需要执行的原因。如果文件来自旧的
standalone 运行，它没有当前 DAG 的 provenance，重新生成通常是正确行为。不要用空文件、
`--touch` 或伪造时间戳掩盖来源差异。

## 报告没有包含大文件

这是预期行为。读者报告不会直接嵌入 H5AD、矩阵和原始图像，应在
`results/report/artifact_manifest.tsv` 中以相对路径、大小和 checksum 引用。

如果链接失效，确认 `results/report/report.html` 没有脱离完整 `results/` 目录单独移动。
Snakemake 技术报告应使用不同文件名，推荐
`results/report/snakemake_report.html`。

## 仍然无法定位问题

收集以下信息再复核：

- 完整命令；
- Snakemake 版本；
- dry-run 输出；
- 失败 rule 的日志；
- 脱敏后的 config 和样本表片段；
- 可用 CPU、内存和磁盘。

不要在 issue 或聊天中粘贴真实样本标识、原始数据路径或未脱敏报告。
