# HW-MD2：一维 Lennard-Jones 分子动力学

将课程的一维简谐振子扩展为 10 粒子的一维 LJ 链，依次模拟自由边界和周期边界。
采用 Velocity Verlet 积分、shifted-force LJ 截断修正和 LJ 约化单位。

计算与可视化独立运行：计算只保存数据，绘图只读取已有结果。
保留原来的模型与默认参数，提供时空图和离线交互联动。可视化不重新积分，也不改写结果文件。

## 文件职责

| 文件 | 职责 |
| --- | --- |
| `example/MD1HO.m`、`example/ForceCal.m` | 课程参考代码，保持原样，不参与 Python 程序运行 |
| `md1d.py` | `SimulationConfig`、LJ 势、成对力计算、Velocity Verlet 和计算命令入口 |
| `trajectory_io.py` | `Trajectory` 数据结构、NPZ 读写、JSON/CSV/XYZ 导出 |
| `analysis.py` | 从轨迹计算能量误差、去质心位移、速度和单粒子内部动能，不积分、不绘图 |
| `visualize.py` | 可视化命令入口，输出时空 PNG 和交互 HTML，不导入或运行计算模块 |
| `interactive.py`、`dashboard/` | 自包含离线交互页的导出、布局与联动逻辑 |
| `report_figures.py` | 基于已有轨迹生成报告静态图和测量汇总，不运行新工况 |
| `REPORT_GUIDE.md` | 报告研究问题、图序、图注草稿、实际发现和解释边界 |
| `requirements.txt` | 当前验证通过的直接依赖版本 |
| `tests/` | 计算独立性、文件读写、周期坐标、分析量物理含义与离线导出检查 |

## uv 与上级共享环境

使用 Python 3.12 和上级目录的共享环境：
`~/workspace/20lesson/homework/mmms_code/.venv`。
不要在 `HW_MD2` 内另建 `.venv`。

```bash
cd ~/workspace/20lesson/homework/mmms_code/HW_MD2

# 首次安装或更新本作业依赖，写入已有共享环境。
uv pip install --python ../.venv/bin/python -r requirements.txt

# 第一步：只计算，默认先 free，后 periodic。
../.venv/bin/python md1d.py

# 第二步：只读取上一步保存的数据，生成时空 PNG 和交互 HTML。
../.venv/bin/python visualize.py
```

只处理一个边界条件：

```bash
../.venv/bin/python md1d.py --boundary free
../.venv/bin/python visualize.py --boundary free
```

自定义输出位置：计算的 `--out` 指定根目录，程序在其中写入 `results/<boundary>/`；
绘图的 `--input` 指向结果根目录、单个边界目录或单个 NPZ，`--out` 直接指定图片目录。

```bash
../.venv/bin/python md1d.py --out /tmp/lj-run
../.venv/bin/python visualize.py --input /tmp/lj-run/results --out /tmp/lj-run/plots

# 单个文件也能独立绘图，无需 metadata.json 或当前计算默认参数。
../.venv/bin/python visualize.py --input results/periodic/trajectory.npz --out /tmp/lj-plots
```

重复执行会覆盖对应输出文件。Matplotlib 使用 Agg 后端，保存图片而不弹出窗口。
交互页包含 Plotly 运行库和全部展示数据，直接用浏览器打开 `plots/dashboard.html`，
不需要服务器、联网或浏览器中的 Python 环境。文件约 7 MB。

```bash
# 可分别生成静态图或交互页；默认 --format both。
../.venv/bin/python visualize.py --format png
../.venv/bin/python visualize.py --format html
```

## 模型、初始化与数值设置

- 粒子数 `n=10`，`epsilon=sigma=mass=1`，`k_B=1`。
- 初始位置 `x[i]=i*2**(1/6)`；第六个粒子速度为 `0.2`，其余速度为零。
- 时间步长 `dt=0.002`，总步数 `10000`，总模拟时间 `20`。
- 每 10 步采样，同时保留初态和终态，默认共 1001 帧。
  自定义步数不能整除采样间隔时，额外保存最终帧。
- 截断半径 `rc=2.5*sigma`，周期盒长 `L=n*spacing≈11.22462`，满足 `L>2*rc`。
- 自由边界下粒子不受盒壁限制；周期边界下使用最小镜像距离并将位置折回 `[0,L)`。
- 遍历所有 `i<j` 的粒子对；每对势能累计一次，对两个粒子施加大小相等、方向相反的力。

默认设置集中在 `md1d.py` 的 `SimulationConfig` 中。计算函数也可直接调用：

```python
from md1d import SimulationConfig, simulate
from trajectory_io import save_trajectory

config = SimulationConfig(nsteps=2000, output_every=10)
trajectory = simulate('periodic', config)  # 只计算并返回数据，不写文件。
save_trajectory(trajectory, 'results/custom-periodic')
```

瞬时位置、速度、力使用一维数组；采样历史使用 `(采样帧数, 粒子数)` 的二维数组。
行表示采样时刻，列表示固定粒子身份。这与参考代码中保存粒子对距离、力的二维矩阵用途不同。

### 势函数

设 `F_LJ(r)=-dU_LJ/dr` 为向外的径向力。`r<rc` 时：

```text
U_LJ(r) = 4*epsilon*((sigma/r)^12 - (sigma/r)^6)
U_SF(r) = U_LJ(r) - U_LJ(rc) + (r-rc)*F_LJ(rc)
F_SF(r) = F_LJ(r) - F_LJ(rc)
```

`r>=rc` 时势能和力均为零，在截断处连续。这是对原始 LJ 势的 shifted-force 修正。

### 初态的物理含义

单粒子速度扰动是本实现的选择，适合观察低能振动、扰动传播与边界效应，不能直接视为热平衡初态。
总动量不为零，预期质心速度为 `0.02`，模拟期间整体漂移 `0.4`。

间距 `2**(1/6)*sigma` 是未修正 LJ 双粒子势的最低点，不是多粒子自由链的精确平衡间距。
自由边缘会松弛；周期等距链的初始合力因左右对称而为零，但不意味着压力为零。

## 结果格式：schema version 1

每种边界在 `results/<boundary>/` 下生成：

| 文件 | 内容 |
| --- | --- |
| `trajectory.npz` | 完整轨迹和内嵌 JSON 元数据，是独立分析的数据来源 |
| `metadata.json` | 同一份元数据的人类可读副本 |
| `thermal.csv` | step、time、势能、总动能、总能量、约化动力学温度 |
| `trajectory.xyz` | 原有多帧普通 XYZ，三维显示坐标为 `(x,0,0)` |

NPZ 使用 `allow_pickle=False` 读取，包含：

- `positions`、`velocities`、`forces`：二维采样历史，默认 `(1001,10)`。
- `positions_unwrapped`：连续实验室坐标，形状相同。
  周期跨界计数在每个积分步更新，不依赖输出帧间的位移推断；自由边界下与 `positions` 相同。
- `thermal`：二维热力学历史，默认 `(1001,6)`。
- `box_length`：盒长标量；自由边界下只是构型参考长度，不是容器。
- `metadata_json`：JSON 字符串，记录完整配置、边界、势函数、积分器、初始化、
  单位、数组轴、坐标约定、采样方式、温度定义和 Python/NumPy 版本。

旧版 NPZ 缺少元数据和连续坐标，读取时会提示重新运行 `md1d.py`，不会猜测历史配置。
普通 XYZ 不携带标准周期盒信息；需要在 OVITO 中设置周期盒，或以后增加 Extended XYZ 导出。

温度为 `T*=mass*sum((v-v_COM)^2)/(n-1)`，扣除一个质心平移自由度。
总动能保留质心平移能，因此始终有 `total_energy=potential_energy+kinetic_energy`。
程序为无恒温器的 NVE 模拟，小体系温度有明显波动。

`analysis.py` 可直接从保存的轨迹计算去质心位移、去质心速度、单粒子内部动能和能量误差。
去质心位移以初态为参考，不代表相对于已松弛平衡构型的位移。

## 时空图与交互联动

`plots/` 中的当前输出为：

- `spacetime.png`：全程时空图，行对应边界条件，列对应位移、速度和内部动能。
- `spacetime_early.png`：同一组量在前 2 个约化时间单位的放大图，更容易观察传播。
- `dashboard.html`：离线交互页；只有一组输入时也能独立使用。
- `initial_conditions.png`：共同初态、周期镜像、截断范围与初始受力。
- `conservation.png`：能量/动量误差、质心运动与内部位移尺度。
- `thermodynamics.png`：动势能交换和瞬时温度，标出有限时段平均值。
- `wave_snapshots.png`：四个保存时刻的速度空间分布，直接对比两种边界。
- `particle_responses.png`：被扰动粒子、邻居、ID 1 和 ID N 的时间响应。
- `report_metrics.csv`：守恒误差、质心位移、温度等实际测量结果。

报告配图围绕“验证守恒 → 区分质心与内部运动 → 比较空间传播 → 解释能量和温度”组织。
推荐正文图序及可直接参考的图注见 [REPORT_GUIDE.md](REPORT_GUIDE.md)。
五张新增报告图为 300 dpi；未增加研究工况，原有计算与交互页保持不变。

旧的六张单独能量、绝对轨迹和温度图片已移除，也不再生成。
新的静态图使用共同时间区间和共同色标，按实际采样时间绘制网格，保留末帧非等间隔采样。
输出文件名固定，单独导出某种边界会替换当前对比视图；如需同时保留多个版本，请指定不同的 `--out`。

交互页支持：

- 时间滑块、播放/暂停、速度切换、早期/全程切换及回到当前区间起点。
- 同一时间下的自由/周期粒子动画、热图游标、粒子读数、运动曲线及能量温度曲线联动。
- 点击热图同时选取时间与粒子；点击粒子图或使用下拉框选择粒子；点击曲线定位时间。
- 缩放任意时间图会同步其他时间图及滑块范围；双击时间图恢复全程。
- 热图切换去质心位移 `u`、去质心速度 `v′` 和单粒子内部动能 `K′`。
  每种量在所有边界、所有帧中使用固定的共同色标。
- 粒子位置切换去质心位移 ×1、×10、实验室连续坐标和盒内折回坐标。
  ×10 仅用于显示放大，热图、曲线及读数始终是真实数值；箭头按固定比例显示去质心速度。
- 两组采样时间不同时，在共同时间范围内按最近保存帧对齐，并标出每个面板的实际采样时间。
  播放按墙钟时间推进并显示已有帧，不插值生成物理数据。

单粒子内部动能为 `K′[i]=mass*(v[i]-v_COM)^2/2`，其求和等于扣除质心运动后的总动能。
它不是局部总能量，也不能直接解释为热流。能量曲线显示相对于初态的动势能变化，
总能量误差 `ΔE` 和瞬时温度同时显示在当前读数中。

## 验证

```bash
../.venv/bin/python -m unittest discover -s tests -v
```

默认完整运行的最大能量误差约为：自由边界 `2.192e-6`、周期边界 `2.123e-6`。
比较不同边界时，不应要求初始总能量相同，因为相互作用的粒子对不同。
