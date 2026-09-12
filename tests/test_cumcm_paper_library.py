import importlib.util
import sys
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "scripts" / "update_cumcm_paper_library.py"
SPEC = importlib.util.spec_from_file_location("cumcm_paper_library", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class PaperLibraryParsingTests(unittest.TestCase):
    def test_discovers_and_sorts_unique_years(self):
        source = "2025全国大学生数学建模竞赛论文展示 2012全国大学生数学建模竞赛论文展示 2025全国大学生数学建模竞赛论文展示"
        self.assertEqual(MODULE.discover_years(source), [2012, 2025])

    def test_extracts_absolute_and_relative_pdf_links(self):
        source = '<a href="/zx/files/A001.pdf">A</a><a href="https://dxs.moe.gov.cn/B002.PDF?x=1">B</a>'
        self.assertEqual(
            MODULE.extract_pdf_urls(source, "https://dxs.moe.gov.cn/zx/a/page.shtml"),
            [
                "https://dxs.moe.gov.cn/zx/files/A001.pdf",
                "https://dxs.moe.gov.cn/B002.PDF?x=1",
            ],
        )

    def test_infers_problem_and_code_from_modern_title(self):
        title = "2025高教社杯全国大学生数学建模竞赛A题论文展示（A196）"
        self.assertEqual(MODULE.infer_problem(title), "A")
        self.assertEqual(MODULE.infer_paper_code(title), "A196")

    def test_infers_problem_and_code_from_legacy_title(self):
        title = "1A23003-A048太阳影子定位的多目标优化模型-数学建模竞赛优秀论文下载"
        self.assertEqual(MODULE.infer_problem(title), "A")
        self.assertEqual(MODULE.infer_paper_code(title), "1A23003-A048")

    def test_excludes_year_landing_page(self):
        item = {
            "id": "123",
            "title": "2023高教社杯全国大学生数学建模竞赛论文展示",
            "link": "https://dxs.moe.gov.cn/example.shtml",
        }
        self.assertFalse(MODULE.is_paper_item(2023, item))


if __name__ == "__main__":
    unittest.main()
