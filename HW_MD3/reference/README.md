# 教师原始文件来源

本实验参考课程仓库的 `0Surface` 目录。完整克隆仅供本地查阅，不纳入本项目 Git，避免提交嵌套仓库及其他课程实验。

| 项目 | 记录 |
|---|---|
| 上游地址 | https://gitee.com/sjtu-liugs/mmmsfall.git |
| 固定版本 | `544b5c5a86f9119b3bc77dfd3f06c065c3168743` |
| 实验目录 | `0Surface` |
| 势文件 | `CuYM.eam.alloy`，文件头注明 Mishin 2001 Cu EAM1 |
| 势文件 SHA-256 | `c0ae005c4680dc6110f6b3ec135f27ab8a2031efacb4a4b52c83e4ec9b54fb24` |

从 GitHub 克隆本项目后，在 `HW_MD3` 中运行：

```bash
bash scripts/fetch_reference.sh
```

该脚本获取并核对固定版本；已有目录会先检查，不会覆盖、切换或清理已有文件。原始输入和势函数的作者及出处仍按上游文件注明的信息引用。

`example/` 中用户下载的课程 HTML、网页附件及示例构型也仅在本地保存。指导来源为课程提供的 *Exercise 1: Molecular Statics Calculation of surface energy (MmmsFall)*：https://notes.sjtu.edu.cn/s/IL0qIhw6F 。
