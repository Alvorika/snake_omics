# HTML 报告

流程提供两种用途不同的 HTML。默认的读者报告适合快速查看一次运行完成了什么；Snakemake
技术报告用于检查 rule、provenance 和被标记为 `report()` 的技术附件。

## 生成方式

完成所选模块并生成读者报告：

```bash
snakemake \
  --snakefile workflow/Snakefile \
  --directory . \
  --cores 8 \
  --sdm conda \
  report
```

输出是 `results/report/report.html`。它包含模块状态、综合 QC、声明过的小型摘要和图片、
provenance，以及所有结果文件的相对链接。H5AD、完整矩阵和原始 H&E 不会写进 HTML。

需要排查 DAG 或留存 Snakemake 自身记录时，可另外生成技术报告：

```bash
snakemake \
  --snakefile workflow/Snakefile \
  --directory . \
  --report results/report/snakemake_report.html \
  report
```

第二条命令不会替代读者报告。两个 HTML 的目标文件名必须不同。技术报告可能包含 rule
provenance 和技术附件，默认只用于内部诊断；除非另行审计，不应把它作为公开交付物。

## 结果目录与体积

`reporting.report.inline_image_max_mb` 控制单张图片的内联上限，
`reporting.report.inline_image_total_max_mb` 控制整份报告内联图片的总上限，
`reporting.report.max_table_preview_rows` 控制表格预览行数。超过上限的图片和所有大文件只
保留相对链接；artifact manifest 记录大小和有上限的 checksum 计算状态。

在受控环境中整体移动运行结果时，应保留整个 `results/` 目录结构。单独移动
`results/report/report.html` 会使 `../qc/`、`../figures/` 等链接失效。

完整 `results/` 不是自动脱敏的公开报告包：部分科学 sidecar 可能保留输入 provenance。
公开前先把准备发布的文件复制到单独的 staging 目录，再扫描其中所有文本结果，并按需
加入本项目已知标识。若生成过 Snakemake 技术报告，不要把它放进公开 staging：

```bash
python scripts/audit_run_outputs.py PATH_TO_PUBLIC_STAGING \
  --project-root . \
  --forbid REPLACE_WITH_PROJECT_IDENTIFIER
```

扫描通过后仍需人工决定哪些矩阵、图片和 sidecar 可以公开。若只发布 HTML，未内联图片和
大文件的相对链接将不可用，但 artifact 名称、大小和 checksum 仍保留在报告中。

## 新模块如何进入报告

新增分析模块时：

1. 在 `workflow/module_registry.py` 注册模块、依赖和说明；
2. 在 rule 中声明输出，并在 `workflow/rules/common.smk` 的 `MODULE_OUTPUTS` 登记；
3. 通用模块状态和 artifact 索引会自动纳入该模块，不需要修改 HTML 生成器；
4. 只有需要精选摘要、表格或图片时，才在
   `workflow/report/report_sections.json` 增加声明。

精选内容只能来自 artifact manifest 已登记的小型 sidecar。JSON 只展示声明过的字段，
TSV/CSV 只展示声明过的列，PNG 受内联大小限制。不要把 H5AD、完整表达矩阵、原始图像或
未清理的日志配置为预览。

模块若能给出比通用判定更准确的结束态，可将唯一的
`results/.../report_summary.json` 登记为该模块 artifact。它必须符合
[sidecar schema](../workflow/report/module_report_summary.schema.json)：版本为 `1.0.0`，
`module` 与登记模块一致，`report_status` 使用现有公开状态；除 `completed` 外均需填写
`status_detail`。报告会优先采用该 sidecar；没有时继续使用现有 QC、2×2 或通用完成态。

## 隐私边界

两个 HTML 都是运行产物，不属于可公开的源码。它们可能显示脱敏后的 `sample_id`、分组和
ROI 名、图中文字、module status、effective config、相对文件名和 checksum。流程会把外部
绝对路径及其文件名统一隐藏为 `<external>/REDACTED`，但无法判断普通文本或图中标签是否
敏感。

公开读者报告前必须人工打开检查；技术报告默认不公开，如确需发布也要独立审计。具体清单
见[数据与隐私](privacy.md)。自动测试和输出扫描不能替代人工脱敏审核。
