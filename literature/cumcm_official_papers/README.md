# CUMCM 官方优秀论文链接库

这是一个只保存元数据和官方链接的可检索索引。中国大学生在线的论文展示页明确限制未经书面许可转载，因此本目录不镜像 PDF 或逐页图片。

- 官方总入口：[https://dxs.moe.gov.cn/zx/hd/sxjm/sxjmlw/qkt_sxjm_lw_lwzs.shtml](https://dxs.moe.gov.cn/zx/hd/sxjm/sxjmlw/qkt_sxjm_lw_lwzs.shtml)
- 本次更新：`2026-09-12T13:39:45+00:00`
- 官网索引覆盖：`2012—2025`
- 收录条目：`140`
- 数据文件：[`index.csv`](index.csv)、[`index.json`](index.json)、[`crawl_report.json`](crawl_report.json)

## 按年份浏览

| 年份 | 论文数 | 题号 | PDF | 逐页图片 | 仅网页 | 失败 |
|---:|---:|---|---:|---:|---:|---:|
| [2012](by_year/2012.md) | 3 | A, B | 2 | 0 | 1 | 0 |
| [2013](by_year/2013.md) | 8 | — | 2 | 0 | 6 | 0 |
| [2014](by_year/2014.md) | 9 | A, B, C | 0 | 4 | 5 | 0 |
| [2015](by_year/2015.md) | 14 | A, B, C, D | 11 | 0 | 3 | 0 |
| [2016](by_year/2016.md) | 5 | A, C, D | 5 | 0 | 0 | 0 |
| [2017](by_year/2017.md) | 2 | A, D | 2 | 0 | 0 | 0 |
| [2018](by_year/2018.md) | 12 | A, B, C, D | 12 | 0 | 0 | 0 |
| [2019](by_year/2019.md) | 9 | B, C, D, E | 0 | 9 | 0 | 0 |
| [2020](by_year/2020.md) | 13 | A, B, C, D, E | 0 | 13 | 0 | 0 |
| [2021](by_year/2021.md) | 18 | A, B, C, D, E | 0 | 17 | 1 | 0 |
| [2022](by_year/2022.md) | 11 | A, B, C, D, E | 0 | 11 | 0 | 0 |
| [2023](by_year/2023.md) | 13 | A, B, C, D, E | 0 | 13 | 0 | 0 |
| [2024](by_year/2024.md) | 16 | A, B, C, D, E | 0 | 16 | 0 | 0 |
| [2025](by_year/2025.md) | 7 | A, B, C, D, E | 0 | 7 | 0 | 0 |

## 使用方式

- 用表格软件打开 `index.csv`，可按年份、题号、论文编号和载体筛选。
- `index.json` 适合程序检索；每条记录包含官方页面、官方 PDF（若有）和抓取状态。
- 重新更新：在项目根目录运行 `python scripts/update_cumcm_paper_library.py`。

## 边界说明

- “历年”指官网当前总入口列出的年份，不代表竞赛自 1992 年以来每一年均有论文公开展示。
- 只收录全国大学生数学建模竞赛论文展示，不混入深圳杯、研究生数学建模竞赛或第三方转载。
- `official_pdf` 表示官方页给出了 PDF；`official_image_pages` 表示官方仅提供逐页图片展示；两者都不在本仓库复制。
- 学校和队员信息若官方论文展示页未提供，本索引不会从非官方来源补写。
