# HW_MD3：Cu 表面能与低温形貌

参考课程实验指导，使用 LAMMPS 计算 FCC Cu 的 (100)、(110)、(111) 表面能，再通过 Wulff 构造回答低温平衡形貌及晶面截距比问题。结构检查和渲染使用 OVITO Python 接口。

## 报告与结果

- [中文实验报告](report/Cu表面能与低温平衡形貌.md)：Obsidian 兼容 Markdown，标题从二级开始，含 OVITO 图片。
- [离线 HTML 报告](report/Cu表面能与低温平衡形貌.html)：下载后可在浏览器中打开，图片与公式已内嵌。
- [报告配图](report/assets/)：与 Markdown 一起复制进 Obsidian 笔记库。
- [表面能完整数据](results/surface_energies.csv)、[结果汇总](results/summary.json)。
- [OVITO 出图说明](ovito/README.md)、[渲染记录](figures/ovito/render_manifest.json)。
- [计算环境与作业记录](results/run_manifest.json)。

中文报告由用户继续微调，完整英文版待用户提供修改稿后再翻译。`report/Homework_answer_EN.md` 是早期英文简答，仅作参考。

助教可直接在 GitHub 中打开上述 Markdown 阅读正文和图片，无需安装 Python、LAMMPS 或 OVITO。原始结果可通过下面的目录查看，重新计算只用于进一步复核。

基准弛豫表面能为 γ111 = 1.2395、γ100 = 1.3453、γ110 = 1.4755 J/m²。在三晶面族 Wulff 近似中，{111}、{100}、{110} 面积占比分别约为 65.93%、29.32%、4.75%。

(111) 最小化因线搜索精度停止，未达到请求的严格力二范数阈值；报告与结果记录保留了这项限制。实验采用 0 K 分子静力学，Wulff 图是几何平衡构造，不是纳米颗粒 MD 轨迹的最终帧。

## 目录

```text
HW_MD3/
├── README.md
├── .gitignore
├── requirements.txt          # 分析与作图依赖
├── requirements-ovito.txt    # OVITO Python 版本
├── scripts/                 # 输入生成、Slurm 运行、分析、渲染与报告导出
├── inputs/                  # 案例清单和输入脚本；重复势文件不提交
├── reference/README.md      # 教师资料来源、固定版本及势文件校验值
├── results/
│   ├── cluster/             # 正式原始日志、最终构型及执行时脚本快照
│   ├── surface_energies.csv
│   ├── summary.json
│   ├── wulff.json
│   ├── run_manifest.json
│   └── raw_checksums.sha256
├── figures/                 # PNG、SVG 与 OVITO 渲染结果
├── ovito/                   # 出图说明；生成的场景文件仅在本地保存
└── report/                  # 当前报告与 assets 配图
```

本地保留、Git 忽略的内容包括：`example/` 课程下载资料、`reference/mmmsfall/` 完整上游克隆、每案例的重复势文件、`results/first_pass/` 第一轮诊断数据、结果压缩包、旧稿和 OVITO 场景。`local_archive/` 保存旧 README 与最初规划，不作为报告提交内容。未删除任何原始模拟数据。

## 安装 Python 依赖

在项目上级目录使用已有 uv 环境：

```bash
uv pip install --python .venv/bin/python -r HW_MD3/requirements.txt
uv pip install --python .venv/bin/python -r HW_MD3/requirements-ovito.txt
```

普通分析只需第一组依赖，OVITO 自动渲染需要第二组。实验使用 Python 3.12、OVITO 3.17.1。所有 Python 脚本按脚本位置确定项目路径。

## 用保存的数据重新分析

以下命令在项目上级目录执行，不提交集群作业：

```bash
.venv/bin/python HW_MD3/scripts/analyze.py
.venv/bin/python HW_MD3/scripts/render_ovito.py
.venv/bin/python HW_MD3/scripts/build_report.py
```

分析脚本读取 `results/cluster/inputs/` 中的正式数据。渲染脚本显式使用 `lammps/dump` 读取器，输出到 `figures/ovito/`，同时在 `ovito/states/` 保存已配置场景。报告导出脚本读取 Markdown 和 `report/assets/`，生成可离线阅读的 HTML。

`report/assets/` 保留作为报告附件的图片副本，防止重新作图时直接改变正在编辑的报告。需要替换图时，再将选定成图复制到该目录。

在 `HW_MD3` 内可校验原始文件：

```bash
sha256sum --check results/raw_checksums.sha256
```

## 重新准备模拟

教师原始仓库不提交到本项目，先获取本实验使用的固定版本：

```bash
bash HW_MD3/scripts/fetch_reference.sh
.venv/bin/python HW_MD3/scripts/prepare.py
```

来源与校验值见 [reference/README.md](reference/README.md)。`prepare.py` 生成 15 个统一参数案例和 2 个教师复现案例，并复制每案例所需的势文件。现有 `inputs/` 已保留本次输入，重新生成会覆盖其中的生成文件，修改参数前应保留自己的版本。

上传 `inputs/` 和 `scripts/run_cluster.sh` 到集群的独立工作目录后，在该目录提交：

```bash
sbatch scripts/run_cluster.sh
```

脚本对应本次思源环境：分区 `64c512g`、1 个 CPU 核、10 分钟上限，模块为 `lammps/20240829.1-intel-2021.4.0`。其他集群需要调整模块和队列；不要在登录节点执行模拟。

## Git 管理范围

正式输入、原始日志、构型、结果表、分析脚本、报告和配图均纳入版本控制。没有全局忽略 `.log`、`.dump`、`.data`，以免误排除实验依据。完整上游克隆和机器本地环境不提交。

SSH 密码、密钥、令牌及本地环境文件不得写入脚本、报告或仓库。根目录 `.gitignore` 已忽略虚拟环境、常见 Python 缓存、本地凭据目录及 `.env` 文件。

本次整理未执行 `git add`、提交或推送。
