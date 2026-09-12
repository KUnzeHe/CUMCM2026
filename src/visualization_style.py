"""CUMCM 论文图的统一视觉规范。"""

from __future__ import annotations

from matplotlib import font_manager
from matplotlib.axes import Axes
import matplotlib.pyplot as plt


COLORS = {
    "navy": "#19324D",
    "blue": "#356A8A",
    "teal": "#2A9D8F",
    "gold": "#E9C46A",
    "vermilion": "#E76F51",
    "charcoal": "#25313C",
    "gray": "#6B7280",
    "light_gray": "#DCE2E8",
    "paper": "#FFFFFF",
}

RADIAL_COLORS = (
    COLORS["navy"],
    COLORS["blue"],
    COLORS["teal"],
    COLORS["gold"],
    COLORS["vermilion"],
)

RADIAL_STYLES = (
    "-",
    (0, (5, 2)),
    (0, (7, 2, 1.5, 2)),
    (0, (2, 2)),
    (0, (9, 2)),
)


def configure_publication_style() -> None:
    """配置适合中文竞赛论文、PDF 和高分辨率位图的 Matplotlib 样式。"""

    installed = {font.name for font in font_manager.fontManager.ttflist}
    latin_candidates = (
        "Times New Roman",
        "Times",
        "Liberation Serif",
        "DejaVu Serif",
    )
    chinese_candidates = (
        "Songti SC",
        "STSong",
        "SimSun",
        "Noto Serif CJK SC",
        "Source Han Serif SC",
    )
    paper_fonts = [family for family in latin_candidates if family in installed]
    paper_fonts.extend(
        family for family in chinese_candidates if family in installed
    )
    if "DejaVu Serif" not in paper_fonts:
        paper_fonts.append("DejaVu Serif")

    plt.rcParams.update(
        {
            "font.family": paper_fonts,
            "mathtext.fontset": "stix",
            "axes.unicode_minus": False,
            "axes.edgecolor": COLORS["charcoal"],
            "axes.labelcolor": COLORS["charcoal"],
            "axes.titlecolor": COLORS["charcoal"],
            "axes.titleweight": "semibold",
            "axes.titlesize": 12.5,
            "axes.labelsize": 10.5,
            "axes.linewidth": 0.9,
            "xtick.color": COLORS["charcoal"],
            "ytick.color": COLORS["charcoal"],
            "xtick.labelsize": 9.5,
            "ytick.labelsize": 9.5,
            "legend.fontsize": 9.0,
            "legend.frameon": False,
            "font.size": 10.0,
            "figure.facecolor": COLORS["paper"],
            "axes.facecolor": COLORS["paper"],
            "savefig.facecolor": COLORS["paper"],
            "savefig.edgecolor": "none",
            "savefig.dpi": 300,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def polish_axis(axis: Axes, *, grid_axis: str = "y") -> None:
    """统一坐标轴、刻度和辅助网格。"""

    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_color(COLORS["charcoal"])
    axis.spines["bottom"].set_color(COLORS["charcoal"])
    axis.tick_params(direction="out", length=4, width=0.8)
    axis.set_axisbelow(True)
    if grid_axis != "none":
        axis.grid(
            True,
            axis=grid_axis,
            color=COLORS["light_gray"],
            linewidth=0.65,
            alpha=0.78,
        )


def add_panel_label(axis: Axes, label: str) -> None:
    """在多面板图中添加稳定、醒目的子图编号。"""

    axis.text(
        -0.12,
        1.04,
        label,
        transform=axis.transAxes,
        ha="left",
        va="bottom",
        fontsize=11.5,
        fontweight="bold",
        color=COLORS["charcoal"],
        clip_on=False,
    )
