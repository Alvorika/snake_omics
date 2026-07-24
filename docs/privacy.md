# 数据与隐私

本仓库只用于保存可复用流程。仓库当前不包含真实数据，也不包含测试数据；测试数据将在
后续选择独立、允许使用的数据集后另行接入。

## 不应进入仓库的内容

- 原始测序数据、Space Ranger 输出和 H&E 原图；
- 真实样本表、ROI 导出和 active config；
- 患者、动物、送样、病历或公司项目标识；
- 带有真实样本名称的图、表、HTML 和日志；
- `work/`、`results/`、`logs/` 和 `.snakemake/`；
- 本机用户名、私有目录、访问令牌和环境密钥；
- 能把脱敏编号反向映射到真实身份的 key。

仓库模板只使用 `sample_01`、`group_a`、`treatment_0` 等中性示例。
项目运行时使用的 `config/config.yaml`、`samples.tsv`、QC 审核表、ROI alias 和 pathway
manifest 已列入 `.gitignore`；应提交对应的 defaults/template，而不是活动项目副本。

## 样本编号

在数据进入流程前完成脱敏。推荐创建项目内部随机 ID，并把真实 ID 的映射表保存在仓库
之外、访问受控的位置。

一个合适的 `sample_id` 应：

- 在当前项目内唯一且稳定；
- 不包含姓名、出生日期、病历号或送样编号；
- 不编码不必要的敏感分组信息；
- 可以安全地出现在文件名、图标题和 HTML 中。

流程会把 `sample_id` 写入结果路径和报告，所以运行后再改名会破坏 provenance。

## 路径与 metadata

模板和文档使用相对路径。活动项目也应尽量使用相对路径，以便 rsync 或移动整个运行目录
后仍然可用。

即使表达矩阵已经脱敏，以下 metadata 仍可能泄露来源：

- slide、library、capture area 和实验日期；
- free-text note；
- 原始图像文件名；
- 外部报告标题；
- 软件日志中的输入目录；
- probe set 或定制 panel 名称。

只保留分析和追溯真正需要的字段。需要内部保存的详细 provenance 不应默认进入公开
HTML。

## 报告与大文件索引

读者报告和 Snakemake 技术报告都可能显示样本 ID、分组/ROI 名、模块状态、图中文字、
effective config、相对文件名和部分结果。生成后必须分别打开检查，不能因为路径扫描通过
就默认内容已经脱敏。外部绝对路径和 basename 会被统一替换，但普通文本和图中标签无法
由流程自动判断。

artifact manifest 只应记录运行目录内的相对结果路径。不要把私有数据源的完整本机路径
复制到公开报告。大型文件可以在受控存储中交付，报告只保留文件角色、大小、checksum 和
相对引用。在受控存储中移动时应保留完整 `results/` 目录结构；单独移动 HTML 会破坏
相对链接。但完整结果树不是自动脱敏的公开 bundle，部分科学 summary 仍可能含输入
provenance。公开随附文件前，先把准备发布的内容复制到独立 staging 目录；不要放入
Snakemake 技术报告，然后执行：

```bash
python scripts/audit_run_outputs.py PATH_TO_PUBLIC_STAGING \
  --project-root . \
  --forbid REPLACE_WITH_PROJECT_IDENTIFIER
```

技术报告默认仅用于内部诊断，不属于默认公开交付。输出扫描只能发现已知本机路径和用户
提供的标识，仍不能替代逐文件人工检查。

## 测试数据

仓库暂不携带测试数据。以后新增测试数据时必须同时满足：

- 来源和许可证允许再分发；
- 不含可识别个体信息；
- 文件足够小；
- 只包含测试所需字段和最少行数；
- 在 README 中记录来源、处理方式和许可证；
- 不由当前真实项目简单抽样得到，除非已完成独立的公开许可审核。

单元测试中的名称也应使用中性合成标签，不能复制真实项目的样本名、处理名或 ROI 组合。

## 同步或发布前检查

发布前至少完成：

1. 查看待同步文件清单；
2. 确认 active config、真实 samples/ROI 表和运行目录未被包含；
3. 搜索本机用户名、项目编号、真实样本名和外部报告标题；
4. 删除 bytecode、缓存和 notebook checkpoint；
5. 打开 README、图和 HTML，人工检查可见文字；
6. 确认没有身份映射 key 或 secret。

本仓库提供物理目录审计和 rsync 排除清单。请在仓库根目录先停止运行任务，再执行：

```bash
python scripts/audit_source_tree.py
rsync -a \
  --exclude-from=scripts/rsync-exclude.txt \
  ./ /path/to/clean-destination/
```

审计会把 `.snakemake/`、结果、日志、Python cache、常见组学数据格式、过大文件、本机
home 路径和已知旧项目标识视为错误。排除清单还会去掉 active config 和 active metadata，
只保留 defaults/template。不要省略末尾 `/` 的来源目录语义，也不要在未核对目标时加入
`--delete`。

自动扫描只能发现已知模式，不能代替人工审核。同步完成后还应在目标目录再运行一次同一
审计，并人工打开 README、模板、图和 HTML。
