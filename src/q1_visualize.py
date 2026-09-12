"""绘制问题 1 温度场和含水率场的论文用图。"""

from __future__ import annotations

import argparse
import os
import tempfile
from pathlib import Path

matplotlib_cache = Path(tempfile.gettempdir()) / "cumcm2026_matplotlib"
matplotlib_cache.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(matplotlib_cache))

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize

from q1_fvm import DEFAULT_ENV_PATH, PROJECT_ROOT, load_environment
from visualization_style import (
    COLORS,
    RADIAL_COLORS,
    RADIAL_STYLES,
    add_panel_label,
    configure_publication_style,
    polish_axis,
)


DEFAULT_RESULTS_PATH = PROJECT_ROOT / "outputs" / "q1" / "data" / "q1_results_full.npz"
DEFAULT_FIGURE_DIR = PROJECT_ROOT / "outputs" / "q1" / "figures" / "fancy"
DEFAULT_RADII_CM = (0.0, 0.5, 1.0, 1.5, 2.0)
DEFAULT_PROFILE_TIMES_S = (300, 600, 900, 1200, 1500, 1800)
PROFILE_CMAP = LinearSegmentedColormap.from_list(
    "cumcm_time",
    (COLORS["navy"], COLORS["teal"], COLORS["gold"], COLORS["vermilion"]),
)
TEMPERATURE_CMAP = "magma_r"


def load_results(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """读取并校验 Q1 的 NPZ 结果。"""

    if not path.is_file():
        raise FileNotFoundError(f"未找到结果文件：{path}")
    with np.load(path) as data:
        required = {"time_s", "radius_cm", "temperature_c", "moisture_kg_kg"}
        missing = required.difference(data.files)
        if missing:
            raise ValueError(f"结果文件缺少字段：{', '.join(sorted(missing))}")
        time_s = np.asarray(data["time_s"], dtype=float)
        radius_cm = np.asarray(data["radius_cm"], dtype=float)
        temperature_c = np.asarray(data["temperature_c"], dtype=float)
        moisture_kg_kg = np.asarray(data["moisture_kg_kg"], dtype=float)

    expected_shape = (time_s.size, radius_cm.size)
    if temperature_c.shape != expected_shape or moisture_kg_kg.shape != expected_shape:
        raise ValueError(
            "温度或含水率矩阵尺寸与时间、半径坐标不一致："
            f"期望 {expected_shape}，实际 {temperature_c.shape}、{moisture_kg_kg.shape}"
        )
    if time_s.ndim != 1 or radius_cm.ndim != 1:
        raise ValueError("时间和半径坐标必须是一维数组")
    if not np.all(np.diff(time_s) > 0) or not np.all(np.diff(radius_cm) > 0):
        raise ValueError("时间和半径坐标必须严格递增")
    if not np.isfinite(temperature_c).all() or not np.isfinite(moisture_kg_kg).all():
        raise ValueError("结果中存在 NaN 或无穷值")
    return time_s, radius_cm, temperature_c, moisture_kg_kg


def nearest_indices(values: np.ndarray, targets: tuple[float, ...]) -> list[int]:
    """查找目标坐标的最近下标，并拒绝超出数据范围的目标。"""

    if min(targets) < values[0] or max(targets) > values[-1]:
        raise ValueError("指定的绘图位置或时刻超出结果范围")
    return [int(np.argmin(np.abs(values - target))) for target in targets]


def save_figure(
    figure: plt.Figure,
    output_directory: Path,
    stem: str,
    formats: tuple[str, ...],
) -> list[Path]:
    output_directory.mkdir(parents=True, exist_ok=True)
    paths = []
    for suffix in formats:
        path = output_directory / f"{stem}.{suffix}"
        figure.savefig(
            path,
            dpi=300,
            bbox_inches="tight",
            metadata={"Creator": "CUMCM2026 visualization pipeline"},
        )
        paths.append(path.resolve())
    plt.close(figure)
    return paths


def radius_label(radius_cm: float, maximum_radius_cm: float) -> str:
    if np.isclose(radius_cm, 0.0):
        return "中心"
    if np.isclose(radius_cm, maximum_radius_cm):
        return "表面"
    return f"r = {radius_cm:g} cm"


def radial_style(radius_cm: float, maximum_radius_cm: float) -> tuple[str, object]:
    """按相对半径固定语义配色，保证删减曲线后表面仍为朱红色。"""

    if maximum_radius_cm <= 0:
        raise ValueError("最大半径必须为正")
    slot = int(round(radius_cm / maximum_radius_cm * (len(RADIAL_COLORS) - 1)))
    slot = min(max(slot, 0), len(RADIAL_COLORS) - 1)
    return RADIAL_COLORS[slot], RADIAL_STYLES[slot]


def draw_temperature_history(
    axis,
    time_s: np.ndarray,
    radius_cm: np.ndarray,
    temperature_c: np.ndarray,
    ambient_temperature,
    selected_radii_cm: tuple[float, ...],
    *,
    compact: bool = False,
) -> None:
    time_min = time_s / 60.0
    indices = nearest_indices(radius_cm, selected_radii_cm)
    for index in indices:
        color, style = radial_style(radius_cm[index], radius_cm[-1])
        axis.plot(
            time_min,
            temperature_c[:, index],
            color=color,
            linestyle=style,
            linewidth=2.15,
            label=radius_label(radius_cm[index], radius_cm[-1]),
        )
    axis.plot(
        time_min,
        np.asarray(ambient_temperature(time_s), dtype=float),
        color=COLORS["gray"],
        linestyle=(0, (3, 2)),
        linewidth=1.55,
        alpha=0.95,
        label="烘房空气",
    )
    axis.set(
        title="径向各位置的升温响应",
        xlabel="时间/min",
        ylabel="温度/°C",
        xlim=(0.0, time_min[-1]),
    )
    polish_axis(axis, grid_axis="y")
    if not compact:
        surface_index = int(np.argmax(radius_cm))
        center_index = int(np.argmin(radius_cm))
        delta = temperature_c[-1, surface_index] - temperature_c[-1, center_index]
        axis.text(
            0.99,
            0.04,
            f"30 min 表面–中心温差  {delta:.1f} °C",
            transform=axis.transAxes,
            ha="right",
            va="bottom",
            fontsize=9.2,
            color=COLORS["vermilion"],
        )
        axis.legend(ncol=3, loc="upper left", handlelength=2.8, columnspacing=1.4)


def draw_moisture_history(
    axis,
    time_s: np.ndarray,
    radius_cm: np.ndarray,
    moisture_kg_kg: np.ndarray,
    selected_radii_cm: tuple[float, ...],
    *,
    compact: bool = False,
) -> None:
    time_min = time_s / 60.0
    indices = nearest_indices(radius_cm, selected_radii_cm)
    for index in indices:
        color, style = radial_style(radius_cm[index], radius_cm[-1])
        axis.plot(
            time_min,
            moisture_kg_kg[:, index],
            color=color,
            linestyle=style,
            linewidth=2.15,
            label=radius_label(radius_cm[index], radius_cm[-1]),
        )
    axis.set(
        title="径向各位置的含水率演化",
        xlabel="时间/min",
        ylabel="干基含水率 C/(kg·kg⁻¹)",
        xlim=(0.0, time_min[-1]),
        ylim=(
            float(np.min(moisture_kg_kg)) - 0.06,
            float(np.max(moisture_kg_kg)) + (0.25 if compact else 0.13),
        ),
    )
    polish_axis(axis, grid_axis="y")
    if not compact:
        surface_index = int(np.argmax(radius_cm))
        center_index = int(np.argmin(radius_cm))
        contrast = moisture_kg_kg[-1, center_index] - moisture_kg_kg[-1, surface_index]
        axis.text(
            0.99,
            0.53,
            f"30 min 表面较中心低  {contrast:.2f} kg/kg",
            transform=axis.transAxes,
            ha="right",
            va="center",
            fontsize=9.2,
            color=COLORS["vermilion"],
        )
        axis.legend(ncol=3, loc="upper left", handlelength=2.8, columnspacing=1.4)


def draw_radial_profiles(
    axes,
    time_s: np.ndarray,
    radius_cm: np.ndarray,
    temperature_c: np.ndarray,
    moisture_kg_kg: np.ndarray,
    profile_times_s: tuple[int, ...],
) -> list:
    indices = nearest_indices(time_s, tuple(float(value) for value in profile_times_s))
    colors = PROFILE_CMAP(np.linspace(0.06, 0.94, len(indices)))
    handles = []
    for index, color in zip(indices, colors, strict=True):
        label = f"{time_s[index] / 60:g} min"
        (line,) = axes[0].plot(
            radius_cm,
            temperature_c[index],
            color=color,
            linewidth=2.05,
            label=label,
        )
        handles.append(line)
        axes[1].plot(
            radius_cm,
            moisture_kg_kg[index],
            color=color,
            linewidth=2.05,
            label=label,
        )

    axes[0].set(
        title="温度从表面向中心逐步响应",
        xlabel="径向位置 r/cm",
        ylabel="温度/°C",
        xlim=(radius_cm[0], radius_cm[-1]),
    )
    axes[1].set(
        title="低含水率层首先形成于表层",
        xlabel="径向位置 r/cm",
        ylabel="干基含水率/(kg·kg⁻¹)",
        xlim=(radius_cm[0], radius_cm[-1]),
    )
    for axis in axes:
        polish_axis(axis, grid_axis="both")
    return handles


def plot_temperature_history(
    time_s: np.ndarray,
    radius_cm: np.ndarray,
    temperature_c: np.ndarray,
    ambient_temperature,
    selected_radii_cm: tuple[float, ...],
) -> plt.Figure:
    figure, axis = plt.subplots(figsize=(8.6, 4.9), layout="constrained")
    draw_temperature_history(
        axis,
        time_s,
        radius_cm,
        temperature_c,
        ambient_temperature,
        selected_radii_cm,
    )
    return figure


def plot_moisture_history(
    time_s: np.ndarray,
    radius_cm: np.ndarray,
    moisture_kg_kg: np.ndarray,
    selected_radii_cm: tuple[float, ...],
) -> plt.Figure:
    figure, axis = plt.subplots(figsize=(8.6, 4.9), layout="constrained")
    draw_moisture_history(
        axis,
        time_s,
        radius_cm,
        moisture_kg_kg,
        selected_radii_cm,
    )
    return figure


def plot_spatiotemporal_heatmaps(
    time_s: np.ndarray,
    radius_cm: np.ndarray,
    temperature_c: np.ndarray,
    moisture_kg_kg: np.ndarray,
) -> plt.Figure:
    figure, axes = plt.subplots(2, 1, figsize=(9.2, 7.4), layout="constrained")
    time_min = time_s / 60.0
    temperature_mesh = axes[0].pcolormesh(
        time_min,
        radius_cm,
        temperature_c.T,
        shading="auto",
        cmap=TEMPERATURE_CMAP,
        rasterized=True,
    )
    axes[0].set(
        title="温度场的时空演化",
        xlabel="时间/min",
        ylabel="径向位置 r/cm",
    )
    temperature_bar = figure.colorbar(temperature_mesh, ax=axes[0], pad=0.02)
    temperature_bar.set_label("温度（°C）")

    moisture_mesh = axes[1].pcolormesh(
        time_min,
        radius_cm,
        moisture_kg_kg.T,
        shading="auto",
        cmap="YlGnBu",
        rasterized=True,
    )
    axes[1].set(
        title="含水率场的时空演化",
        xlabel="时间/min",
        ylabel="径向位置 r/cm",
    )
    moisture_bar = figure.colorbar(moisture_mesh, ax=axes[1], pad=0.02)
    moisture_bar.set_label("干基含水率 C（kg/kg）")
    ticks = np.arange(0.0, time_min[-1] + 0.1, 5.0)
    for axis in axes:
        axis.set_xlim(0.0, time_min[-1])
        axis.set_xticks(ticks)
        axis.tick_params(direction="out", length=4, width=0.8)
    add_panel_label(axes[0], "(a)")
    add_panel_label(axes[1], "(b)")
    return figure


def plot_cross_section_evolution(
    time_s: np.ndarray,
    radius_cm: np.ndarray,
    temperature_c: np.ndarray,
    moisture_kg_kg: np.ndarray,
    profile_times_s: tuple[int, ...],
) -> plt.Figure:
    """用圆形小多图展示轴对称横截面上温度场与含水率场的演化。"""

    indices = nearest_indices(time_s, tuple(float(value) for value in profile_times_s))
    theta = np.linspace(0.0, 2.0 * np.pi, 241)
    theta_mesh, radius_mesh = np.meshgrid(theta, radius_cm)
    temperature_norm = Normalize(
        vmin=float(np.min(temperature_c[indices])),
        vmax=float(np.max(temperature_c[indices])),
    )
    moisture_norm = Normalize(
        vmin=float(np.min(moisture_kg_kg[indices])),
        vmax=float(np.max(moisture_kg_kg[indices])),
    )

    figure = plt.figure(figsize=(13.2, 5.7), layout="constrained")
    grid = figure.add_gridspec(
        2,
        len(indices) + 2,
        width_ratios=[0.18] + [1.0] * len(indices) + [0.08],
        wspace=0.04,
        hspace=0.08,
    )
    outline_theta = np.linspace(0.0, 2.0 * np.pi, 361)
    rows = (
        (temperature_c, TEMPERATURE_CMAP, temperature_norm),
        (moisture_kg_kg, "YlGnBu", moisture_norm),
    )
    for row, (field, cmap, norm) in enumerate(rows):
        label_axis = figure.add_subplot(grid[row, 0])
        label_axis.set_axis_off()
        label_axis.text(
            0.5,
            0.5,
            "温度场 T" if row == 0 else "干基含水率 C",
            rotation=90,
            va="center",
            ha="center",
            fontsize=10.8,
            fontweight="bold",
            color=COLORS["charcoal"],
        )
        for column, index in enumerate(indices):
            axis = figure.add_subplot(grid[row, column + 1], projection="polar")
            field_mesh = np.repeat(field[index, :, np.newaxis], theta.size, axis=1)
            axis.pcolormesh(
                theta_mesh,
                radius_mesh,
                field_mesh,
                shading="gouraud",
                cmap=cmap,
                norm=norm,
                rasterized=True,
            )
            axis.plot(
                outline_theta,
                np.full_like(outline_theta, radius_cm[-1]),
                color=COLORS["navy"],
                linewidth=1.0,
                alpha=0.72,
            )
            axis.scatter(
                [0.0],
                [0.0],
                s=9,
                color=COLORS["paper"],
                edgecolor=COLORS["navy"],
                linewidth=0.5,
                zorder=4,
            )
            axis.set_ylim(0.0, radius_cm[-1])
            axis.set_axis_off()
            if row == 0:
                axis.set_title(
                    f"{time_s[index] / 60:g} min",
                    fontsize=10.4,
                    fontweight="bold",
                    color=COLORS["charcoal"],
                    pad=8,
                )

        colorbar_axis = figure.add_subplot(grid[row, -1])
        colorbar = figure.colorbar(
            ScalarMappable(norm=norm, cmap=cmap),
            cax=colorbar_axis,
        )
        if row == 0:
            colorbar.set_label("温度/°C")
        else:
            colorbar.set_label("C/(kg·kg⁻¹)")

    figure.suptitle(
        "问题1｜圆柱横截面的径向传热与水分迁移",
        x=0.012,
        ha="left",
        fontsize=16.5,
        fontweight="bold",
        color=COLORS["charcoal"],
    )
    figure.text(
        0.5,
        0.012,
        "统一色标便于跨时刻比较；圆心为 r = 0，外缘为物料表面",
        ha="center",
        va="bottom",
        fontsize=9.2,
        color=COLORS["gray"],
    )
    return figure


def plot_state_trajectories(
    time_s: np.ndarray,
    radius_cm: np.ndarray,
    temperature_c: np.ndarray,
    moisture_kg_kg: np.ndarray,
    selected_radii_cm: tuple[float, ...],
    profile_times_s: tuple[int, ...],
) -> plt.Figure:
    """绘制不同径向位置在温度-含水率状态空间中的演化轨迹。"""

    radial_indices = nearest_indices(radius_cm, selected_radii_cm)
    marker_indices = nearest_indices(
        time_s, tuple(float(value) for value in profile_times_s)
    )
    figure, axis = plt.subplots(figsize=(9.0, 5.7), layout="constrained")
    for radial_index in radial_indices:
        radius = radius_cm[radial_index]
        color, style = radial_style(radius, radius_cm[-1])
        axis.plot(
            temperature_c[:, radial_index],
            moisture_kg_kg[:, radial_index],
            color=color,
            linestyle=style,
            linewidth=2.25,
            label=radius_label(radius, radius_cm[-1]),
        )
        axis.scatter(
            temperature_c[marker_indices, radial_index],
            moisture_kg_kg[marker_indices, radial_index],
            s=31,
            facecolors=COLORS["paper"],
            edgecolors=color,
            linewidth=1.15,
            zorder=3,
        )
        axis.scatter(
            temperature_c[marker_indices[-1], radial_index],
            moisture_kg_kg[marker_indices[-1], radial_index],
            s=54,
            marker=">",
            color=color,
            edgecolor=COLORS["paper"],
            linewidth=0.75,
            zorder=4,
        )

    axis.set(
        title="局部状态轨迹：升温与水分迁移同步推进",
        xlabel="温度 T/°C",
        ylabel="干基含水率 C/(kg·kg⁻¹)",
    )
    polish_axis(axis, grid_axis="both")
    axis.legend(ncol=3, loc="lower left", handlelength=2.8, columnspacing=1.4)
    axis.text(
        0.99,
        0.97,
        "空心圆：每 5 min　▶：30 min",
        transform=axis.transAxes,
        ha="right",
        va="top",
        fontsize=9.2,
        color=COLORS["gray"],
    )
    figure.suptitle(
        "问题1｜不同径向位置的温度-含水率状态路径",
        x=0.01,
        ha="left",
        fontsize=16.0,
        fontweight="bold",
        color=COLORS["charcoal"],
    )
    return figure


def plot_field_landscapes(
    time_s: np.ndarray,
    radius_cm: np.ndarray,
    temperature_c: np.ndarray,
    moisture_kg_kg: np.ndarray,
) -> plt.Figure:
    """绘制温度场与含水率场的三维时空地形图。"""

    stride = max(1, time_s.size // 120)
    sample_indices = np.unique(
        np.append(np.arange(0, time_s.size, stride, dtype=int), time_s.size - 1)
    )
    time_mesh, radius_mesh = np.meshgrid(time_s[sample_indices] / 60.0, radius_cm)
    figure = plt.figure(figsize=(12.4, 5.6), layout="constrained")
    fields = (
        (
            temperature_c[sample_indices].T,
            TEMPERATURE_CMAP,
            "温度场 T",
            "温度/°C",
        ),
        (moisture_kg_kg[sample_indices].T, "YlGnBu", "含水率场 C", "C/(kg·kg⁻¹)"),
    )
    for position, (field, cmap, title, zlabel) in enumerate(fields, start=1):
        axis = figure.add_subplot(1, 2, position, projection="3d")
        surface = axis.plot_surface(
            time_mesh,
            radius_mesh,
            field,
            cmap=cmap,
            linewidth=0.0,
            antialiased=True,
            rcount=radius_cm.size,
            ccount=min(121, sample_indices.size),
        )
        surface.set_rasterized(True)
        z_min = float(np.min(field))
        z_max = float(np.max(field))
        z_span = max(z_max - z_min, 1e-9)
        floor = z_min - 0.08 * z_span
        axis.set(
            title=title,
            xlabel="时间/min",
            ylabel="径向位置 r/cm",
            zlim=(floor, z_max + 0.03 * z_span),
        )
        axis.view_init(elev=27, azim=-126)
        axis.set_box_aspect((1.55, 1.0, 0.82))
        axis.grid(False)
        axis.xaxis.pane.set_facecolor((1.0, 1.0, 1.0, 0.0))
        axis.yaxis.pane.set_facecolor((1.0, 1.0, 1.0, 0.0))
        axis.zaxis.pane.set_facecolor((1.0, 1.0, 1.0, 0.0))
        colorbar = figure.colorbar(surface, ax=axis, shrink=0.62, pad=0.04)
        colorbar.set_label(zlabel)
        axis.text2D(
            0.01,
            0.96,
            f"({chr(96 + position)})",
            transform=axis.transAxes,
            fontsize=11.0,
            fontweight="bold",
            color=COLORS["charcoal"],
        )

    figure.suptitle(
        "问题1｜双场时空演化地形",
        x=0.01,
        ha="left",
        fontsize=16.0,
        fontweight="bold",
        color=COLORS["charcoal"],
    )
    return figure


def plot_radial_profiles(
    time_s: np.ndarray,
    radius_cm: np.ndarray,
    temperature_c: np.ndarray,
    moisture_kg_kg: np.ndarray,
    profile_times_s: tuple[int, ...],
) -> plt.Figure:
    figure, axes = plt.subplots(1, 2, figsize=(11.6, 4.8), layout="constrained")
    handles = draw_radial_profiles(
        axes,
        time_s,
        radius_cm,
        temperature_c,
        moisture_kg_kg,
        profile_times_s,
    )
    add_panel_label(axes[0], "(a)")
    add_panel_label(axes[1], "(b)")
    figure.legend(
        handles=handles,
        labels=[handle.get_label() for handle in handles],
        loc="outside upper center",
        ncol=min(6, len(handles)),
        handlelength=2.4,
        columnspacing=1.35,
    )
    return figure


def plot_q1_overview(
    time_s: np.ndarray,
    radius_cm: np.ndarray,
    temperature_c: np.ndarray,
    moisture_kg_kg: np.ndarray,
    ambient_temperature,
    profile_times_s: tuple[int, ...],
) -> plt.Figure:
    """绘制Q1论文主视觉：时间响应与径向传播的四面板叙事图。"""

    figure, axes = plt.subplots(2, 2, figsize=(12.2, 8.0), layout="constrained")
    key_radii = (0.0, 1.0, 2.0)
    draw_temperature_history(
        axes[0, 0],
        time_s,
        radius_cm,
        temperature_c,
        ambient_temperature,
        key_radii,
        compact=True,
    )
    draw_moisture_history(
        axes[0, 1],
        time_s,
        radius_cm,
        moisture_kg_kg,
        key_radii,
        compact=True,
    )
    handles = draw_radial_profiles(
        axes[1],
        time_s,
        radius_cm,
        temperature_c,
        moisture_kg_kg,
        profile_times_s,
    )
    for axis, label in zip(axes.flat, ("(a)", "(b)", "(c)", "(d)"), strict=True):
        add_panel_label(axis, label)

    axes[0, 0].legend(
        ncol=2,
        loc="upper left",
        fontsize=8.4,
        handlelength=2.5,
        columnspacing=1.1,
    )
    axes[0, 1].legend(
        ncol=3,
        loc="upper left",
        fontsize=8.4,
        handlelength=2.5,
        columnspacing=1.0,
    )
    figure.legend(
        handles=handles,
        labels=[handle.get_label() for handle in handles],
        loc="outside lower center",
        ncol=6,
        handlelength=2.4,
        columnspacing=1.2,
        title="剖面时刻",
    )
    figure.suptitle(
        "问题1｜预热阶段的径向传热与水分迁移",
        x=0.01,
        ha="left",
        fontsize=16.5,
        fontweight="bold",
        color=COLORS["charcoal"],
    )
    return figure


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS_PATH)
    parser.add_argument("--environment", type=Path, default=DEFAULT_ENV_PATH)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_FIGURE_DIR)
    parser.add_argument(
        "--radii",
        type=float,
        nargs="+",
        default=DEFAULT_RADII_CM,
        metavar="CM",
        help="时间曲线使用的径向位置/cm",
    )
    parser.add_argument(
        "--profile-times",
        type=int,
        nargs="+",
        default=DEFAULT_PROFILE_TIMES_S,
        metavar="S",
        help="径向剖面使用的时刻/s",
    )
    parser.add_argument(
        "--formats",
        choices=("png", "pdf", "svg"),
        nargs="+",
        default=("png",),
        help="输出图片格式",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_argument_parser().parse_args(argv)
    configure_publication_style()
    time_s, radius_cm, temperature_c, moisture_kg_kg = load_results(args.results)
    ambient_temperature, _, _ = load_environment(args.environment, float(time_s[-1]))
    selected_radii_cm = tuple(float(value) for value in args.radii)
    if len(selected_radii_cm) > len(RADIAL_COLORS):
        raise ValueError(f"位置曲线最多绘制 {len(RADIAL_COLORS)} 条，以保证图例清晰")
    profile_times_s = tuple(int(value) for value in args.profile_times)
    formats = tuple(dict.fromkeys(args.formats))

    figures = {
        "q1_overview": plot_q1_overview(
            time_s,
            radius_cm,
            temperature_c,
            moisture_kg_kg,
            ambient_temperature,
            profile_times_s,
        ),
        "temperature_time_series": plot_temperature_history(
            time_s,
            radius_cm,
            temperature_c,
            ambient_temperature,
            selected_radii_cm,
        ),
        "moisture_time_series": plot_moisture_history(
            time_s, radius_cm, moisture_kg_kg, selected_radii_cm
        ),
        "spatiotemporal_heatmaps": plot_spatiotemporal_heatmaps(
            time_s, radius_cm, temperature_c, moisture_kg_kg
        ),
        "cross_section_evolution": plot_cross_section_evolution(
            time_s,
            radius_cm,
            temperature_c,
            moisture_kg_kg,
            profile_times_s,
        ),
        "state_trajectories": plot_state_trajectories(
            time_s,
            radius_cm,
            temperature_c,
            moisture_kg_kg,
            selected_radii_cm,
            profile_times_s,
        ),
        "field_landscapes": plot_field_landscapes(
            time_s,
            radius_cm,
            temperature_c,
            moisture_kg_kg,
        ),
        "radial_profiles": plot_radial_profiles(
            time_s,
            radius_cm,
            temperature_c,
            moisture_kg_kg,
            profile_times_s,
        ),
    }

    generated = []
    for stem, figure in figures.items():
        generated.extend(save_figure(figure, args.output_dir, stem, formats))
    print("已生成：")
    for path in generated:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
