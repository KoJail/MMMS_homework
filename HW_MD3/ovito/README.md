## OVITO 自动出图

使用官方 `ovito==3.17.1` Python 包，通过 `input_format='lammps/dump'` 明确指定读取器，不依赖桌面软件保留的 XYZ 导入设置。

在上级工作目录运行：

```bash
.venv/bin/python HW_MD3/scripts/render_ovito.py
```

图片输出到 `HW_MD3/figures/ovito/`；对应的 `.ovito` 场景保存到 `HW_MD3/ovito/states/`。场景包含数据来源、读取格式、修改器及相机。它们是 OVITO Python 渲染的图片，可以在报告中标明“使用 OVITO 渲染”，不能称为 GUI 截图。

已生成的输出：

| 文件 | 内容 |
|---|---|
| `Cu100_bulk_cna.png` | 周期体相，Adaptive CNA 结构识别 |
| `Cu100_slab_cna.png` | 真空间隙与双表面，Adaptive CNA |
| `Cu100_slab_energy.png` | Cu(100) 薄板按逐原子势能着色 |
| `Cu110_slab_energy.png` | Cu(110) 薄板按逐原子势能着色 |
| `Cu111_slab_energy.png` | Cu(111) 薄板按逐原子势能着色 |

所有势能图使用相同的 −3.55 至 −2.80 eV 色标，覆盖三个晶面的实际逐原子势能范围。CNA 的 `Other` 表示局域结构未归入完整 FCC/HCP 等类别，不能直接解释为非晶。

已验证 Cu(100) 体相的 1728 个原子全部为 FCC；表面模型中 1584 个为 FCC，144 个为 Other。图片和场景可由上述命令重新生成。

数据读取自正式结果 `results/cluster/inputs/`，不使用 `results/first_pass/`。

官方接口参考：[安装](https://www.ovito.org/manual/python/introduction/getting_started/installation.html)、[显式指定读取格式](https://www.ovito.org/manual/python/modules/ovito_io.html)、[渲染](https://www.ovito.org/manual/python/modules/ovito_vis.html)。
