#!/usr/bin/env python3
"""Build a metadata-only index of officially displayed CUMCM papers.

The official publication pages state that republication requires written
permission.  This tool therefore stores metadata and official links only; it
does not mirror PDFs or page images.
"""

from __future__ import annotations

import argparse
import csv
import html
import json
import re
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass, replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable
from urllib.error import HTTPError, URLError
from urllib.parse import urljoin, urlparse
from urllib.request import Request, urlopen


OFFICIAL_INDEX_URL = (
    "https://dxs.moe.gov.cn/zx/hd/sxjm/sxjmlw/"
    "qkt_sxjm_lw_lwzs.shtml"
)
OFFICIAL_SITE_ROOT = "https://dxs.moe.gov.cn/zx/"
USER_AGENT = (
    "Mozilla/5.0 (compatible; CUMCMOfficialPaperIndex/1.0; "
    "+https://www.mcm.edu.cn/)"
)
DEFAULT_OUTPUT = Path("literature/cumcm_official_papers")
TRANSPORT = "auto"


@dataclass(frozen=True)
class PaperRecord:
    year: int
    problem: str
    paper_code: str
    title: str
    source: str
    author: str
    publish_date: str
    content_type: str
    official_page_url: str
    direct_pdf_url: str
    page_image_count: int
    official_id: str
    catalog_name: str
    metadata_url: str
    alternate_official_page_urls: str
    duplicate_official_ids: str
    duplicate_count: int
    status: str
    error: str


def fetch_once_with_curl(url: str, timeout: float) -> bytes:
    result = subprocess.run(
        [
            "curl",
            "--fail",
            "--location",
            "--silent",
            "--show-error",
            "--max-time",
            str(timeout),
            "--user-agent",
            USER_AGENT,
            url,
        ],
        check=False,
        capture_output=True,
    )
    if result.returncode != 0:
        message = result.stderr.decode("utf-8", errors="replace").strip()
        raise RuntimeError(f"curl 失败（{result.returncode}）：{message}")
    return result.stdout


def fetch_once_with_urllib(url: str, timeout: float) -> bytes:
    request = Request(url, headers={"User-Agent": USER_AGENT})
    with urlopen(request, timeout=timeout) as response:
        return response.read()


def fetch_bytes(url: str, timeout: float = 30.0, retries: int = 3) -> bytes:
    last_error: Exception | None = None
    use_curl = TRANSPORT == "curl" or (TRANSPORT == "auto" and shutil.which("curl"))
    for attempt in range(retries):
        try:
            if use_curl:
                return fetch_once_with_curl(url, timeout)
            return fetch_once_with_urllib(url, timeout)
        except (HTTPError, URLError, TimeoutError, RuntimeError) as exc:
            last_error = exc
            if attempt + 1 < retries:
                time.sleep(0.4 * (2**attempt))
    assert last_error is not None
    raise last_error


def fetch_text(url: str, timeout: float = 30.0) -> str:
    data = fetch_bytes(url, timeout=timeout)
    return data.decode("utf-8", errors="replace")


def fetch_json(url: str, timeout: float = 30.0) -> Any:
    return json.loads(fetch_text(url, timeout=timeout))


def discover_years(index_html: str) -> list[int]:
    years = {
        int(value)
        for value in re.findall(
            r"(20\d{2})全国大学生数学建模竞赛论文展示", index_html
        )
    }
    return sorted(years)


def catalog_json_url(year: int) -> str:
    return (
        f"{OFFICIAL_SITE_ROOT}json/hd/sxjm/sxjmlw/"
        f"{year}qgdxssxjmjslwzs/"
    )


def unwrap_json_object(payload: Any) -> dict[str, Any]:
    if isinstance(payload, list) and payload and isinstance(payload[0], dict):
        return payload[0]
    if isinstance(payload, dict):
        return payload
    raise ValueError("官方 JSON 结构不是对象或非空对象数组")


def load_catalog(year: int, timeout: float) -> list[dict[str, Any]]:
    base_url = catalog_json_url(year)
    next_url = base_url
    seen_pages: set[str] = set()
    records: list[dict[str, Any]] = []

    while next_url:
        if next_url in seen_pages:
            raise ValueError(f"{year} 年目录出现循环分页：{next_url}")
        seen_pages.add(next_url)
        page = unwrap_json_object(fetch_json(next_url, timeout=timeout))
        data = page.get("data") or []
        if not isinstance(data, list):
            raise ValueError(f"{year} 年目录 data 字段不是列表")
        records.extend(item for item in data if isinstance(item, dict))

        next_page = page.get("nextPageUrl")
        if not next_page or next_page == "#":
            next_url = ""
        elif str(next_page).startswith(("http://", "https://")):
            next_url = str(next_page)
        else:
            next_url = urljoin(base_url, str(next_page))

    unique: dict[str, dict[str, Any]] = {}
    for item in records:
        key = str(item.get("id") or item.get("artUrl") or item.get("link") or "")
        if not key:
            continue
        unique[key] = item
    return list(unique.values())


def extract_pdf_urls(content: str, base_url: str) -> list[str]:
    decoded = html.unescape(content or "")
    candidates = re.findall(
        r'''(?:href\s*=\s*["']([^"']+?\.pdf(?:\?[^"']*)?)["'])''',
        decoded,
        flags=re.IGNORECASE,
    )
    if not candidates:
        candidates = re.findall(
            r"https?://[^\s\"'<>]+?\.pdf(?:\?[^\s\"'<>]*)?",
            decoded,
            flags=re.IGNORECASE,
        )
    urls: list[str] = []
    for candidate in candidates:
        url = urljoin(base_url, candidate.strip())
        if url not in urls:
            urls.append(url)
    return urls


def infer_problem(title: str) -> str:
    patterns = (
        r"([A-E])题",
        r"（([A-E])\d+）",
        r"\(([A-E])\d+\)",
        r"^\d([A-E])\d",
        r"(?:^|[-_])([A-E])\d{2,}",
    )
    for pattern in patterns:
        match = re.search(pattern, title, flags=re.IGNORECASE)
        if match:
            return match.group(1).upper()
    return ""


def infer_paper_code(title: str) -> str:
    for pattern in (
        r"（([A-E]\d+)）",
        r"\(([A-E]\d+)\)",
        r"^([0-9][A-E][A-Z0-9_-]+?)(?=[\u4e00-\u9fff“])",
    ):
        match = re.search(pattern, title, flags=re.IGNORECASE)
        if match:
            return match.group(1).upper()
    return ""


def page_image_count(article: dict[str, Any]) -> int:
    images = article.get("list") or article.get("extendjson") or []
    if not isinstance(images, list):
        return 0
    return sum(
        1
        for item in images
        if isinstance(item, dict)
        and any(item.get(key) for key in ("imagepath", "imagePath", "url"))
    )


def make_record(
    year: int,
    item: dict[str, Any],
    timeout: float,
) -> PaperRecord:
    metadata_url = str(item.get("link") or "")
    base_title = str(item.get("title") or "").strip()
    try:
        if not metadata_url:
            raise ValueError("目录条目缺少元数据 URL")
        article = unwrap_json_object(fetch_json(metadata_url, timeout=timeout))
        title = str(article.get("title") or base_title).strip()
        page_url = str(
            article.get("pcLink")
            or article.get("artUrl")
            or item.get("artUrl")
            or ""
        )
        content = str(article.get("content") or "")
        pdf_urls = extract_pdf_urls(content, page_url or OFFICIAL_SITE_ROOT)
        image_count = page_image_count(article)
        content_type = str(article.get("contentType") or item.get("contentType") or "")
        if pdf_urls:
            status = "official_pdf"
        elif image_count:
            status = "official_image_pages"
        elif page_url:
            status = "official_page_only"
        else:
            status = "metadata_only"
        return PaperRecord(
            year=year,
            problem=infer_problem(title),
            paper_code=infer_paper_code(title),
            title=title,
            source=str(article.get("source") or item.get("source") or "").strip(),
            author=str(article.get("author") or item.get("author") or "").strip(),
            publish_date=str(
                article.get("publishDate") or item.get("publishDate") or ""
            ).strip(),
            content_type=content_type,
            official_page_url=page_url,
            direct_pdf_url=" | ".join(pdf_urls),
            page_image_count=image_count,
            official_id=str(article.get("id") or item.get("id") or ""),
            catalog_name=str(
                article.get("catalogName") or item.get("catalogName") or ""
            ).strip(),
            metadata_url=metadata_url,
            alternate_official_page_urls="",
            duplicate_official_ids="",
            duplicate_count=1,
            status=status,
            error="",
        )
    except Exception as exc:  # preserve the catalog entry and its failure reason
        page_url = str(item.get("artUrl") or "")
        return PaperRecord(
            year=year,
            problem=infer_problem(base_title),
            paper_code=infer_paper_code(base_title),
            title=base_title,
            source=str(item.get("source") or "").strip(),
            author=str(item.get("author") or "").strip(),
            publish_date=str(item.get("publishDate") or "").strip(),
            content_type=str(item.get("contentType") or ""),
            official_page_url=page_url,
            direct_pdf_url="",
            page_image_count=0,
            official_id=str(item.get("id") or ""),
            catalog_name=str(item.get("catalogName") or "").strip(),
            metadata_url=metadata_url,
            alternate_official_page_urls="",
            duplicate_official_ids="",
            duplicate_count=1,
            status="metadata_fetch_failed",
            error=f"{type(exc).__name__}: {exc}",
        )


def is_paper_item(year: int, item: dict[str, Any]) -> bool:
    title = re.sub(r"\s+", "", str(item.get("title") or ""))
    generic_titles = {
        f"{year}全国大学生数学建模竞赛论文展示",
        f"{year}高教社杯全国大学生数学建模竞赛论文展示",
    }
    if title in generic_titles and not infer_paper_code(title):
        return False
    return bool(item.get("id") and item.get("link"))


def normalized_title_key(title: str) -> str:
    value = re.sub(r"\s+", "", title).lower()
    value = re.sub(r"[-—_]*(?:优秀)?论文(?:下载)?$", "", value)
    value = re.sub(r"[-—_]+$", "", value)
    return value


def logical_key(record: PaperRecord) -> tuple[int, str, str]:
    if record.paper_code:
        return (record.year, "code", record.paper_code.upper())
    return (record.year, "title", normalized_title_key(record.title))


def collapse_duplicate_records(records: list[PaperRecord]) -> list[PaperRecord]:
    groups: dict[tuple[int, str, str], list[PaperRecord]] = {}
    for record in records:
        groups.setdefault(logical_key(record), []).append(record)

    status_priority = {
        "official_pdf": 4,
        "official_image_pages": 3,
        "official_page_only": 2,
        "metadata_only": 1,
        "metadata_fetch_failed": 0,
    }
    collapsed: list[PaperRecord] = []
    for group in groups.values():
        primary = max(
            group,
            key=lambda record: (
                status_priority.get(record.status, -1),
                record.publish_date,
                record.page_image_count,
            ),
        )
        alternate_pages = sorted(
            {
                record.official_page_url
                for record in group
                if record.official_page_url
                and record.official_page_url != primary.official_page_url
            }
        )
        all_ids = sorted({record.official_id for record in group if record.official_id})
        collapsed.append(
            replace(
                primary,
                alternate_official_page_urls=" | ".join(alternate_pages),
                duplicate_official_ids=" | ".join(all_ids) if len(group) > 1 else "",
                duplicate_count=len(group),
            )
        )
    return collapsed


def sort_key(record: PaperRecord) -> tuple[int, str, str, str]:
    return (record.year, record.problem or "Z", record.paper_code, record.title)


def md_escape(value: str) -> str:
    return value.replace("|", "\\|").replace("\n", " ").strip()


def atomic_write_text(path: Path, text: str, encoding: str = "utf-8") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text, encoding=encoding)
    temporary.replace(path)


def write_csv(path: Path, records: Iterable[PaperRecord]) -> None:
    rows = [asdict(record) for record in records]
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(PaperRecord.__dataclass_fields__))
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def year_summary(records: list[PaperRecord]) -> list[dict[str, Any]]:
    summary: list[dict[str, Any]] = []
    for year in sorted({record.year for record in records}):
        subset = [record for record in records if record.year == year]
        summary.append(
            {
                "year": year,
                "papers": len(subset),
                "problems": sorted({record.problem for record in subset if record.problem}),
                "official_pdf": sum(r.status == "official_pdf" for r in subset),
                "official_image_pages": sum(
                    r.status == "official_image_pages" for r in subset
                ),
                "page_only": sum(r.status == "official_page_only" for r in subset),
                "failed": sum(r.status == "metadata_fetch_failed" for r in subset),
            }
        )
    return summary


def write_year_pages(output: Path, records: list[PaperRecord]) -> None:
    years_dir = output / "by_year"
    for year in sorted({record.year for record in records}):
        subset = [record for record in records if record.year == year]
        lines = [
            f"# {year} 年 CUMCM 官方优秀论文索引",
            "",
            "正文不在本仓库复制；请通过“官方页面”或“官方 PDF”访问。",
            "",
            "| 题号 | 论文编号 | 标题 | 载体 | 官方页面 | 官方 PDF |",
            "|---|---|---|---|---|---|",
        ]
        for record in subset:
            page = (
                f"[打开]({record.official_page_url})"
                if record.official_page_url
                else "—"
            )
            pdf = (
                f"[PDF]({record.direct_pdf_url.split(' | ')[0]})"
                if record.direct_pdf_url
                else "—"
            )
            carrier = {
                "official_pdf": "PDF",
                "official_image_pages": f"逐页图片（{record.page_image_count} 页）",
                "official_page_only": "网页",
                "metadata_fetch_failed": "元数据失败",
                "metadata_only": "仅元数据",
            }.get(record.status, record.status)
            lines.append(
                "| {problem} | {code} | {title} | {carrier} | {page} | {pdf} |".format(
                    problem=record.problem or "—",
                    code=md_escape(record.paper_code) or "—",
                    title=md_escape(record.title),
                    carrier=carrier,
                    page=page,
                    pdf=pdf,
                )
            )
        lines.append("")
        atomic_write_text(years_dir / f"{year}.md", "\n".join(lines))


def write_readme(
    output: Path,
    records: list[PaperRecord],
    retrieved_at: str,
    discovered_years: list[int],
) -> None:
    summary = year_summary(records)
    lines = [
        "# CUMCM 官方优秀论文链接库",
        "",
        "这是一个只保存元数据和官方链接的可检索索引。中国大学生在线的论文展示页明确限制未经书面许可转载，因此本目录不镜像 PDF 或逐页图片。",
        "",
        f"- 官方总入口：[{OFFICIAL_INDEX_URL}]({OFFICIAL_INDEX_URL})",
        f"- 本次更新：`{retrieved_at}`",
        f"- 官网索引覆盖：`{min(discovered_years)}—{max(discovered_years)}`",
        f"- 收录条目：`{len(records)}`",
        "- 数据文件：[`index.csv`](index.csv)、[`index.json`](index.json)、[`crawl_report.json`](crawl_report.json)",
        "",
        "## 按年份浏览",
        "",
        "| 年份 | 论文数 | 题号 | PDF | 逐页图片 | 仅网页 | 失败 |",
        "|---:|---:|---|---:|---:|---:|---:|",
    ]
    for item in summary:
        year = item["year"]
        lines.append(
            f"| [{year}](by_year/{year}.md) | {item['papers']} | "
            f"{', '.join(item['problems']) or '—'} | {item['official_pdf']} | "
            f"{item['official_image_pages']} | {item['page_only']} | {item['failed']} |"
        )
    lines.extend(
        [
            "",
            "## 使用方式",
            "",
            "- 用表格软件打开 `index.csv`，可按年份、题号、论文编号和载体筛选。",
            "- `index.json` 适合程序检索；每条记录包含官方页面、官方 PDF（若有）和抓取状态。",
            "- 重新更新：在项目根目录运行 `python scripts/update_cumcm_paper_library.py`。",
            "",
            "## 边界说明",
            "",
            "- “历年”指官网当前总入口列出的年份，不代表竞赛自 1992 年以来每一年均有论文公开展示。",
            "- 只收录全国大学生数学建模竞赛论文展示，不混入深圳杯、研究生数学建模竞赛或第三方转载。",
            "- `official_pdf` 表示官方页给出了 PDF；`official_image_pages` 表示官方仅提供逐页图片展示；两者都不在本仓库复制。",
            "- 学校和队员信息若官方论文展示页未提供，本索引不会从非官方来源补写。",
            "",
        ]
    )
    atomic_write_text(output / "README.md", "\n".join(lines))


def write_outputs(
    output: Path,
    records: list[PaperRecord],
    discovered_years: list[int],
    retrieved_at: str,
) -> None:
    output.mkdir(parents=True, exist_ok=True)
    records = sorted(records, key=sort_key)
    write_csv(output / "index.csv", records)
    atomic_write_text(
        output / "index.json",
        json.dumps([asdict(record) for record in records], ensure_ascii=False, indent=2)
        + "\n",
    )
    report = {
        "official_index_url": OFFICIAL_INDEX_URL,
        "retrieved_at_utc": retrieved_at,
        "discovered_years": discovered_years,
        "record_count": len(records),
        "duplicate_groups_collapsed": sum(record.duplicate_count > 1 for record in records),
        "duplicate_source_entries_collapsed": sum(
            record.duplicate_count - 1 for record in records
        ),
        "year_summary": year_summary(records),
        "status_counts": {
            status: sum(record.status == status for record in records)
            for status in sorted({record.status for record in records})
        },
        "failed_records": [
            asdict(record) for record in records if record.status.endswith("failed")
        ],
        "copyright_policy": (
            "Metadata and official links only; paper files and page images are not mirrored."
        ),
    }
    atomic_write_text(
        output / "crawl_report.json",
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
    )
    write_year_pages(output, records)
    write_readme(output, records, retrieved_at, discovered_years)


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--year", type=int, action="append", dest="years")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--timeout", type=float, default=30.0)
    parser.add_argument(
        "--transport",
        choices=("auto", "curl", "urllib"),
        default="auto",
        help="联网方式；auto 优先使用 curl，以兼容 macOS 的系统证书链。",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    global TRANSPORT
    args = parse_args(argv or sys.argv[1:])
    TRANSPORT = args.transport
    if TRANSPORT == "curl" and not shutil.which("curl"):
        raise RuntimeError("指定了 curl 传输，但系统中找不到 curl")
    index_html = fetch_text(OFFICIAL_INDEX_URL, timeout=args.timeout)
    discovered_years = discover_years(index_html)
    if not discovered_years:
        raise RuntimeError("未能从官方总入口识别任何年份")

    selected_years = sorted(set(args.years or discovered_years))
    unknown = sorted(set(selected_years) - set(discovered_years))
    if unknown:
        raise ValueError(f"以下年份未出现在官方总入口：{unknown}")

    catalogs: dict[int, list[dict[str, Any]]] = {}
    for year in selected_years:
        catalogs[year] = load_catalog(year, timeout=args.timeout)
        print(f"{year}: 发现 {len(catalogs[year])} 条", flush=True)

    jobs = [
        (year, item)
        for year in selected_years
        for item in catalogs[year]
        if is_paper_item(year, item)
    ]
    records: list[PaperRecord] = []
    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as executor:
        futures = {
            executor.submit(make_record, year, item, args.timeout): (year, item)
            for year, item in jobs
        }
        for index, future in enumerate(as_completed(futures), start=1):
            records.append(future.result())
            if index % 25 == 0 or index == len(futures):
                print(f"元数据：{index}/{len(futures)}", flush=True)

    expected = len(jobs)
    if len(records) != expected:
        raise RuntimeError(f"条目数不一致：目录 {expected}，结果 {len(records)}")

    records = collapse_duplicate_records(records)

    retrieved_at = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    write_outputs(args.output, records, discovered_years, retrieved_at)
    failures = sum(record.status == "metadata_fetch_failed" for record in records)
    print(
        f"完成：{len(records)} 条，{failures} 条失败；输出到 {args.output}",
        flush=True,
    )
    return 0 if failures == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
