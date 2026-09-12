"""绘制问题 2 变物性热–湿耦合模型的论文用图。

脚本读取 ``outputs/q2/data`` 中未四舍五入的 NPY 结果，复用全项目统一的
宋体/Times New Roman 论文视觉规范。当前数值结果来自尚未完成 Kirchhoff
界面修正的求解器；由于问题 2 的前三小时含水率仍高于 1 kg/kg，图可用于
展示耦合过程和排版，但正式定稿仍应在离散修正后重生并核对。
"""

from __future__ import annotations

import argparse
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path

matplotlib_cache = Path(tempfile.gettempdir()) / "cumcm2026_matplotlib"
matplotlib_cache.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(matplotlib_cache))

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

from q2_coupled import D_f, P, cp_f, k_f, load_env, rho_f
from visualization_style import (
    COLORS,
    RADIAL_COLORS,
    RADIAL_STYLES,
    add_panel_label,
    configure_publication_style,
    polish_axis,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA_DIR = PROJECT_ROOT / "outputs" / "q2" / "data"
DEFAULT_Q1_RESULTS = PROJECT_ROOT / "outputs" / "q1" / "data" / "q1_results_full.npz"
DEFAULT_FIGURE_DIR = PROJECT_ROOT / "outputs" / "q2" / "figures" / "fancy"
SELECTED_RADII_CM = (0.0, 0.5, 1.0, 1.5, 2.0)
PROFILE_TIMES_H = (0.5, 1.0, 1.5, 2.0, 2.5, 3.0)
CROSS_SECTION_TIMES_H = (0.5, 1.5, 2.5, 3.0)
TEMPERATURE_CMAP = "magma_r"
MOISTURE_CMAP = "YlGnBu"
TIME_CMAP = matplotlib.colors.LinearSegmentedColormap.from_list(
    "cumcm_q2_time",
    (COLORS["navy"], COLORS["blue"], COLORS["teal"], COLORS["gold"], COLORS["vermilion"]),
)


@dataclass(frozen=True)
class Q2Data:
    """Q2 温度、含水率和环境边界的完整时间序列。"""

    time_s: np.ndarray
    radius_cm: np.ndarray
    temperature_c: np.ndarray
    moisture: np.ndarray
    ambient_temperature_c: np.ndarray
    ambient_moisture: np.ndarray


@dataclass(frozen=True)
class Q1Data:
    """用于前 30 min 同尺度对比的 Q1 结果。"""

    time_s: np.ndarray
    radius_cm: np.ndarray
    temperature_c: np.ndarray
    moisture: np.ndarray


def load_q2_results(data_directory: Path, output_dt_s: float = 1.0) -> Q2Data:
    """读取、补入精确初态并校验 Q2 的 NPY 结果。"""

    if output_dt_s <= 0:
        raise ValueError("输出时间间隔必须为正")
    paths = {
        "temperature": data_directory / "Q2_T.npy",
        "moisture": data_directory / "Q2_C.npy",
        "radius": data_directory / "Q2_r.npy",
    }
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("缺少 Q2 数据文件：" + "、".join(missing))

    temperature_c = np.asarray(np.load(paths["temperature"]), dtype=float)
    moisture = np.asarray(np.load(paths["moisture"]), dtype=float)
    radius_cm = np.asarray(np.load(paths["radius"]), dtype=float)
    if temperature_c.shape != moisture.shape:
        raise ValueError("Q2 温度和含水率矩阵尺寸不一致")
    if temperature_c.ndim != 2 or temperature_c.shape[1] != radius_cm.size:
        raise ValueError("Q2 场矩阵尺寸与半径坐标不一致")
    if not np.all(np.diff(radius_cm) > 0):
        raise ValueError("Q2 半径坐标必须严格递增")
    if not np.isfinite(temperature_c).all() or not np.isfinite(moisture).all():
        raise ValueError("Q2 结果中存在 NaN 或无穷值")
    if np.any(moisture <= 0):
        raise ValueError("Q2 含水率必须为正")

    stored_time_s = np.arange(1, temperature_c.shape[0] + 1, dtype=float) * output_dt_s
    if not np.isclose(stored_time_s[-1], P.t_end):
        raise ValueError(
            f"Q2 数据终点应为 {P.t_end:g} s，实际为 {stored_time_s[-1]:g} s"
        )

    time_s = np.concatenate(([0.0], stored_time_s))
    temperature_c = np.vstack((np.full(radius_cm.size, P.T0), temperature_c))
    moisture = np.vstack((np.full(radius_cm.size, P.C0), moisture))
    ambient_temperature, ambient_moisture = load_env()
    ambient_temperature_c = np.asarray(ambient_temperature(time_s), dtype=float)
    ambient_moisture_values = np.asarray(ambient_moisture(time_s), dtype=float)

    if np.any(np.diff(moisture, axis=0) > 1e-8):
        raise ValueError("检测到 Q2 含水率随时间明显反弹，需先核查源数据")
    if np.any(np.diff(moisture, axis=1) > 1e-8):
        raise ValueError("检测到 Q2 含水率沿半径反向增大，需先核查源数据")

    return Q2Data(
        time_s=time_s,
        radius_cm=radius_cm,
        temperature_c=temperature_c,
        moisture=moisture,
        ambient_temperature_c=ambient_temperature_c,
        ambient_moisture=ambient_moisture_values,
    )


def load_q1_results(path: Path) -> Q1Data:
    """读取 Q1 结构化结果，用于与 Q2 做同时间、同半径对比。"""

    if not path.is_file():
        raise FileNotFoundError(f"未找到 Q1 结构化结果：{path}")
    with np.load(path) as result:
        required = {"time_s", "radius_cm", "temperature_c", "moisture_kg_kg"}
        missing = required.difference(result.files)
        if missing:
            raise ValueError(f"Q1 结果缺少字段：{', '.join(sorted(missing))}")
        data = Q1Data(
            time_s=np.asarray(result["time_s"], dtype=float),
            radius_cm=np.asarray(result["radius_cm"], dtype=float),
            temperature_c=np.asarray(result["temperature_c"], dtype=float),
            moisture=np.asarray(result["moisture_kg_kg"], dtype=float),
        )
    expected = (data.time_s.size, data.radius_cm.size)
    if data.temperature_c.shape != expected or data.moisture.shape != expected:
        raise ValueError("Q1 场矩阵尺寸与时间、半径坐标不一致")
    return data


def nearest_indices(values: np.ndarray, targets: tuple[float, ...]) -> list[int]:
    if min(targets) < values[0] or max(targets) > values[-1]:
        raise ValueError("指定的绘图位置或时刻超出结果范围")
    return [int(np.argmin(np.abs(values - target))) for target in targets]


def radial_label(radius_cm: float, maximum_radius_cm: float) -> str:
    if np.isclose(radius_cm, 0.0):
        return "中心"
    if np.isclose(radius_cm, maximum_radius_cm):
        return "表面"
    return f"r = {radius_cm:g} cm"


def radial_style(radius_cm: float, maximum_radius_cm: float) -> tuple[str, object]:
    slot = int(round(radius_cm / maximum_radius_cm * (len(RADIAL_COLORS) - 1)))
    slot = min(max(slot, 0), len(RADIAL_COLORS) - 1)
    return RADIAL_COLORS[slot], RADIAL_STYLES[slot]


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


def draw_temperature_histories(axis, data: Q2Data, *, legend: bool = True) -> None:
    time_h = data.time_s / 3600.0
    indices = nearest_indices(data.radius_cm, SELECTED_RADII_CM)
    for index in indices:
        color, style = radial_style(data.radius_cm[index], data.radius_cm[-1])
        axis.plot(
            time_h,
            data.temperature_c[:, index],
            color=color,
            linestyle=style,
            linewidth=2.05,
            label=radial_label(data.radius_cm[index], data.radius_cm[-1]),
        )
    axis.plot(
        time_h,
        data.ambient_temperature_c,
        color=COLORS["gray"],
        linestyle=(0, (3, 2)),
        linewidth=1.45,
        label="烘房空气",
    )
    axis.set(
        xlabel="干燥时间/h",
        ylabel="温度/°C",
        xlim=(0.0, time_h[-1]),
    )
    polish_axis(axis, grid_axis="y")
    if legend:
        axis.legend(ncol=3, loc="lower right", handlelength=2.7, columnspacing=1.2)


def draw_moisture_histories(axis, data: Q2Data, *, legend: bool = True) -> None:
    time_h = data.time_s / 3600.0
    indices = nearest_indices(data.radius_cm, SELECTED_RADII_CM)
    for index in indices:
        color, style = radial_style(data.radius_cm[index], data.radius_cm[-1])
        axis.plot(
            time_h,
            data.moisture[:, index],
            color=color,
            linestyle=style,
            linewidth=2.05,
            label=radial_label(data.radius_cm[index], data.radius_cm[-1]),
        )
    axis.plot(
        time_h,
        data.ambient_moisture,
        color=COLORS["gray"],
        linestyle=(0, (3, 2)),
        linewidth=1.35,
        label="等效环境含水率",
    )
    axis.set(
        xlabel="干燥时间/h",
        ylabel="干基含水率 C/(kg·kg⁻¹)",
        xlim=(0.0, time_h[-1]),
    )
    polish_axis(axis, grid_axis="y")
    if legend:
        axis.legend(ncol=3, loc="lower left", handlelength=2.7, columnspacing=1.2)


def draw_temperature_heatmap(axis, data: Q2Data, *, colorbar: bool = True):
    mesh = axis.pcolormesh(
        data.time_s / 3600.0,
        data.radius_cm,
        data.temperature_c.T,
        shading="auto",
        cmap=TEMPERATURE_CMAP,
        rasterized=True,
    )
    axis.set(
        xlabel="干燥时间/h",
        ylabel="径向位置 r/cm",
        xlim=(0.0, data.time_s[-1] / 3600.0),
        ylim=(data.radius_cm[0], data.radius_cm[-1]),
    )
    axis.set_yticks([0.0, 0.5, 1.0, 1.5, 2.0], ["中心", "0.5", "1.0", "1.5", "表面"])
    if colorbar:
        bar = axis.figure.colorbar(mesh, ax=axis, pad=0.025, fraction=0.045)
        bar.set_label("温度/°C")
    return mesh


def draw_moisture_heatmap(axis, data: Q2Data, *, colorbar: bool = True):
    mesh = axis.pcolormesh(
        data.time_s / 3600.0,
        data.radius_cm,
        data.moisture.T,
        shading="auto",
        cmap=MOISTURE_CMAP,
        rasterized=True,
    )
    axis.set(
        xlabel="干燥时间/h",
        ylabel="径向位置 r/cm",
        xlim=(0.0, data.time_s[-1] / 3600.0),
        ylim=(data.radius_cm[0], data.radius_cm[-1]),
    )
    axis.set_yticks([0.0, 0.5, 1.0, 1.5, 2.0], ["中心", "0.5", "1.0", "1.5", "表面"])
    if colorbar:
        bar = axis.figure.colorbar(mesh, ax=axis, pad=0.025, fraction=0.045)
        bar.set_label("干基含水率 C/(kg·kg⁻¹)")
    return mesh


def plot_q2_overview(
    data: Q2Data, output_directory: Path, formats: tuple[str, ...]
) -> list[Path]:
    """绘制 Q2 温度、含水率历史和连续场四面板总览。"""

    figure, axes = plt.subplots(2, 2, figsize=(13.0, 9.0), constrained_layout=True)
    draw_temperature_histories(axes[0, 0], data, legend=True)
    axes[0, 0].set(title="径向各位置的温度响应", ylim=(27.0, 51.5))
    add_panel_label(axes[0, 0], "a")

    draw_moisture_histories(axes[0, 1], data, legend=True)
    axes[0, 1].set(title="径向各位置的含水率演化", ylim=(0.0, 2.68))
    add_panel_label(axes[0, 1], "b")

    temperature_mesh = draw_temperature_heatmap(axes[1, 0], data, colorbar=False)
    axes[1, 0].set_title("温度场时空演化")
    temperature_bar = figure.colorbar(
        temperature_mesh, ax=axes[1, 0], pad=0.02, fraction=0.045
    )
    temperature_bar.set_label("温度/°C")
    add_panel_label(axes[1, 0], "c")

    moisture_mesh = draw_moisture_heatmap(axes[1, 1], data, colorbar=False)
    axes[1, 1].set_title("含水率场时空演化")
    moisture_bar = figure.colorbar(
        moisture_mesh, ax=axes[1, 1], pad=0.02, fraction=0.045
    )
    moisture_bar.set_label("C/(kg·kg⁻¹)")
    add_panel_label(axes[1, 1], "d")

    figure.suptitle("问题 2｜固定尺寸变物性热–湿耦合过程", fontsize=15.2)
    return save_figure(figure, output_directory, "q2_overview", formats)


def plot_radial_histories(
    data: Q2Data, output_directory: Path, formats: tuple[str, ...]
) -> list[Path]:
    """绘制温度与含水率的多半径时间曲线。"""

    figure, axes = plt.subplots(1, 2, figsize=(12.6, 5.1), constrained_layout=True)
    draw_temperature_histories(axes[0], data, legend=True)
    axes[0].set(title="升温响应：后期内外温度趋于一致", ylim=(27.0, 51.5))
    axes[0].annotate(
        f"3 h 表面–中心温差\n{data.temperature_c[-1, -1] - data.temperature_c[-1, 0]:.3f} °C",
        xy=(3.0, data.temperature_c[-1, -1]),
        xytext=(2.03, 42.2),
        arrowprops={"arrowstyle": "->", "color": COLORS["vermilion"], "lw": 1.1},
        color=COLORS["vermilion"],
        ha="center",
    )
    add_panel_label(axes[0], "a")

    draw_moisture_histories(axes[1], data, legend=True)
    axes[1].set(title="失水响应：径向含水率梯度持续扩大", ylim=(0.0, 2.68))
    axes[1].annotate(
        f"3 h 中心–表面差值\n{data.moisture[-1, 0] - data.moisture[-1, -1]:.3f} kg/kg",
        xy=(3.0, data.moisture[-1, -1]),
        xytext=(2.04, 0.58),
        arrowprops={"arrowstyle": "->", "color": COLORS["vermilion"], "lw": 1.1},
        color=COLORS["vermilion"],
        ha="center",
    )
    add_panel_label(axes[1], "b")
    return save_figure(figure, output_directory, "q2_radial_histories", formats)


def plot_spatiotemporal_fields(
    data: Q2Data, output_directory: Path, formats: tuple[str, ...]
) -> list[Path]:
    """绘制温度与含水率双场时空图。"""

    figure, axes = plt.subplots(2, 1, figsize=(11.2, 8.0), constrained_layout=True)
    draw_temperature_heatmap(axes[0], data, colorbar=True)
    axes[0].set_title("温度场：低温浅色，高温深色")
    add_panel_label(axes[0], "a")

    draw_moisture_heatmap(axes[1], data, colorbar=True)
    axes[1].set_title("含水率场：低含水率浅色，高含水率深色")
    add_panel_label(axes[1], "b")
    figure.suptitle("Q2 温度—含水率双场时空演化", fontsize=14.5)
    return save_figure(figure, output_directory, "q2_spatiotemporal_fields", formats)


def plot_radial_profiles(
    data: Q2Data, output_directory: Path, formats: tuple[str, ...]
) -> list[Path]:
    """绘制六个代表时刻的温度和含水率径向剖面。"""

    indices = nearest_indices(
        data.time_s, tuple(hour * 3600.0 for hour in PROFILE_TIMES_H)
    )
    colors = TIME_CMAP(np.linspace(0.03, 0.97, len(indices)))
    figure, axes = plt.subplots(1, 2, figsize=(12.2, 5.2), constrained_layout=True)
    for index, color in zip(indices, colors, strict=True):
        label = f"{data.time_s[index] / 3600:g} h"
        axes[0].plot(
            data.radius_cm,
            data.temperature_c[index],
            color=color,
            linewidth=2.05,
            label=label,
        )
        axes[1].plot(
            data.radius_cm,
            data.moisture[index],
            color=color,
            linewidth=2.05,
            label=label,
        )

    axes[0].set(
        title="温度径向剖面",
        xlabel="径向位置 r/cm",
        ylabel="温度/°C",
        xlim=(0.0, 2.0),
    )
    axes[1].set(
        title="含水率径向剖面",
        xlabel="径向位置 r/cm",
        ylabel="干基含水率 C/(kg·kg⁻¹)",
        xlim=(0.0, 2.0),
        ylim=(0.9, 2.58),
    )
    for panel, axis in zip(("a", "b"), axes, strict=True):
        polish_axis(axis, grid_axis="y")
        add_panel_label(axis, panel)
    axes[0].legend(title="干燥时间", ncol=2, loc="lower right")
    axes[1].legend(title="干燥时间", ncol=2, loc="lower left")
    return save_figure(figure, output_directory, "q2_radial_profiles", formats)


def _circular_field(
    radial_distance: np.ndarray,
    radius_cm: np.ndarray,
    values: np.ndarray,
) -> np.ma.MaskedArray:
    field = np.interp(radial_distance, radius_cm, values)
    return np.ma.masked_where(radial_distance > radius_cm[-1], field)


def plot_cross_section_evolution(
    data: Q2Data, output_directory: Path, formats: tuple[str, ...]
) -> list[Path]:
    """绘制温度、含水率两行圆截面小多图。"""

    indices = nearest_indices(
        data.time_s, tuple(hour * 3600.0 for hour in CROSS_SECTION_TIMES_H)
    )
    radius = data.radius_cm[-1]
    coordinates = np.linspace(-radius, radius, 281)
    xx, yy = np.meshgrid(coordinates, coordinates)
    radial_distance = np.sqrt(xx**2 + yy**2)
    temperature_norm = Normalize(
        vmin=float(data.temperature_c[indices].min()),
        vmax=float(data.temperature_c[indices].max()),
    )
    moisture_norm = Normalize(
        vmin=float(data.moisture[indices].min()),
        vmax=float(data.moisture[indices].max()),
    )
    temperature_cmap = plt.get_cmap(TEMPERATURE_CMAP).copy()
    moisture_cmap = plt.get_cmap(MOISTURE_CMAP).copy()
    temperature_cmap.set_bad("white", alpha=0.0)
    moisture_cmap.set_bad("white", alpha=0.0)

    figure, axes = plt.subplots(2, 4, figsize=(13.0, 6.65), constrained_layout=True)
    for column, index in enumerate(indices):
        temperature_image = axes[0, column].imshow(
            _circular_field(radial_distance, data.radius_cm, data.temperature_c[index]),
            extent=(-radius, radius, -radius, radius),
            origin="lower",
            cmap=temperature_cmap,
            norm=temperature_norm,
            interpolation="bilinear",
        )
        moisture_image = axes[1, column].imshow(
            _circular_field(radial_distance, data.radius_cm, data.moisture[index]),
            extent=(-radius, radius, -radius, radius),
            origin="lower",
            cmap=moisture_cmap,
            norm=moisture_norm,
            interpolation="bilinear",
        )
        for row in range(2):
            axes[row, column].add_patch(
                plt.Circle(
                    (0.0, 0.0), radius, fill=False, color=COLORS["charcoal"], linewidth=0.75
                )
            )
            axes[row, column].scatter([0.0], [0.0], s=5, color="white", alpha=0.85)
            axes[row, column].set(aspect="equal")
            axes[row, column].axis("off")
        axes[0, column].set_title(
            f"{data.time_s[index] / 3600:g} h\n"
            f"中心 {data.temperature_c[index, 0]:.1f} °C · 表面 {data.temperature_c[index, -1]:.1f} °C",
            fontsize=10.3,
        )
        axes[1, column].set_title(
            f"中心 {data.moisture[index, 0]:.3f} · 表面 {data.moisture[index, -1]:.3f}",
            fontsize=10.3,
        )

    temperature_bar = figure.colorbar(
        temperature_image, ax=axes[0, :].tolist(), pad=0.018, fraction=0.025
    )
    temperature_bar.set_label("温度/°C")
    moisture_bar = figure.colorbar(
        moisture_image, ax=axes[1, :].tolist(), pad=0.018, fraction=0.025
    )
    moisture_bar.set_label("干基含水率 C/(kg·kg⁻¹)")
    axes[0, 0].text(
        -0.18,
        0.5,
        "温度场",
        transform=axes[0, 0].transAxes,
        rotation=90,
        ha="center",
        va="center",
        fontsize=11.5,
        fontweight="semibold",
        color=COLORS["charcoal"],
    )
    axes[1, 0].text(
        -0.18,
        0.5,
        "含水率场",
        transform=axes[1, 0].transAxes,
        rotation=90,
        ha="center",
        va="center",
        fontsize=11.5,
        fontweight="semibold",
        color=COLORS["charcoal"],
    )
    figure.suptitle("Q2 圆截面双场演化", fontsize=14.5)
    return save_figure(figure, output_directory, "q2_cross_section_evolution", formats)


def plot_coupling_trajectories(
    data: Q2Data, output_directory: Path, formats: tuple[str, ...]
) -> list[Path]:
    """绘制中心与表面的温度—含水率耦合轨迹。"""

    time_h = data.time_s / 3600.0
    sample = np.unique(np.linspace(0, time_h.size - 1, 550).astype(int))
    figure, axes = plt.subplots(1, 2, figsize=(11.4, 5.25), constrained_layout=True)
    scatter = None
    for panel, index, title, color in (
        ("a", 0, "中心状态轨迹", COLORS["navy"]),
        ("b", -1, "表面状态轨迹", COLORS["vermilion"]),
    ):
        axis = axes[0] if index == 0 else axes[1]
        axis.plot(
            data.temperature_c[:, index],
            data.moisture[:, index],
            color=COLORS["light_gray"],
            linewidth=3.0,
            zorder=1,
        )
        scatter = axis.scatter(
            data.temperature_c[sample, index],
            data.moisture[sample, index],
            c=time_h[sample],
            cmap=TIME_CMAP,
            s=24,
            linewidth=0.0,
            zorder=3,
        )
        axis.scatter(
            [data.temperature_c[0, index], data.temperature_c[-1, index]],
            [data.moisture[0, index], data.moisture[-1, index]],
            s=(42, 68),
            color=(COLORS["gold"], color),
            edgecolor="white",
            linewidth=0.9,
            zorder=4,
        )
        axis.annotate(
            "初态",
            xy=(data.temperature_c[0, index], data.moisture[0, index]),
            xytext=(30.2, 2.42),
            arrowprops={"arrowstyle": "->", "color": COLORS["gold"], "lw": 1.0},
            color=COLORS["gold"],
        )
        axis.annotate(
            "3 h",
            xy=(data.temperature_c[-1, index], data.moisture[-1, index]),
            xytext=(45.4, data.moisture[-1, index] + 0.24),
            arrowprops={"arrowstyle": "->", "color": color, "lw": 1.0},
            color=color,
        )
        axis.set(
            title=title,
            xlabel="温度/°C",
            ylabel="干基含水率 C/(kg·kg⁻¹)",
            xlim=(27.3, 50.8),
            ylim=(0.9, 2.62),
        )
        polish_axis(axis, grid_axis="both")
        add_panel_label(axis, panel)

    colorbar = figure.colorbar(scatter, ax=axes.tolist(), pad=0.02, fraction=0.035)
    colorbar.set_label("干燥时间/h")
    figure.suptitle("温度—含水率耦合轨迹", fontsize=14.2)
    return save_figure(figure, output_directory, "q2_coupling_trajectories", formats)


def plot_property_evolution(
    data: Q2Data, output_directory: Path, formats: tuple[str, ...]
) -> list[Path]:
    """绘制 Q2 四种状态相关物性的中心/表面演化。"""

    time_h = data.time_s / 3600.0
    density = rho_f(data.moisture)
    heat_capacity = cp_f(data.moisture)
    conductivity = k_f(data.moisture)
    diffusivity_scaled = D_f(data.moisture, data.temperature_c) * 1e9
    figure, axes = plt.subplots(2, 2, figsize=(11.4, 8.0), constrained_layout=True)
    panels = (
        (density, "有效密度", "ρ/(kg·m⁻³)"),
        (heat_capacity, "比热容", "cp/(J·kg⁻¹·K⁻¹)"),
        (conductivity, "导热系数", "k/(W·m⁻¹·K⁻¹)"),
        (diffusivity_scaled, "有效水分扩散系数", "D/(10⁻⁹ m²·s⁻¹)"),
    )
    for panel, axis, (values, title, ylabel) in zip(
        ("a", "b", "c", "d"), axes.ravel(), panels, strict=True
    ):
        axis.plot(time_h, values[:, 0], color=COLORS["navy"], linewidth=2.2, label="中心")
        axis.plot(
            time_h,
            values[:, -1],
            color=COLORS["vermilion"],
            linewidth=2.2,
            linestyle=(0, (6, 2)),
            label="表面",
        )
        axis.set(title=title, xlabel="干燥时间/h", ylabel=ylabel, xlim=(0.0, 3.0))
        polish_axis(axis, grid_axis="y")
        add_panel_label(axis, panel)
        axis.legend(loc="best")
    figure.suptitle("Q2 变物性指纹：含水率变化反向作用于传热传质", fontsize=14.5)
    return save_figure(figure, output_directory, "q2_property_evolution", formats)


def plot_q1_q2_comparison(
    q2: Q2Data,
    q1: Q1Data,
    output_directory: Path,
    formats: tuple[str, ...],
) -> list[Path]:
    """在前 30 min 同尺度比较 Q1 常物性与 Q2 变物性结果。"""

    q1_indices = nearest_indices(q1.radius_cm, (0.0, q1.radius_cm[-1]))
    q2_indices = nearest_indices(q2.radius_cm, (0.0, q2.radius_cm[-1]))
    if not np.allclose(q1.radius_cm, q2.radius_cm):
        raise ValueError("Q1 与 Q2 半径网格不一致，不能直接进行同尺度比较")
    common_end_s = min(float(q1.time_s[-1]), 1800.0)
    q1_mask = q1.time_s <= common_end_s
    q2_mask = q2.time_s <= common_end_s
    q1_time_min = q1.time_s[q1_mask] / 60.0
    q2_time_min = q2.time_s[q2_mask] / 60.0

    figure, axes = plt.subplots(1, 2, figsize=(12.2, 5.15), constrained_layout=True)
    for q1_index, q2_index, color, label in zip(
        q1_indices,
        q2_indices,
        (COLORS["navy"], COLORS["vermilion"]),
        ("中心", "表面"),
        strict=True,
    ):
        axes[0].plot(
            q1_time_min,
            q1.temperature_c[q1_mask, q1_index],
            color=color,
            linewidth=1.65,
            linestyle=(0, (3, 2)),
            alpha=0.88,
            label=f"Q1 {label}",
        )
        axes[0].plot(
            q2_time_min,
            q2.temperature_c[q2_mask, q2_index],
            color=color,
            linewidth=2.25,
            label=f"Q2 {label}",
        )
        axes[1].plot(
            q1_time_min,
            q1.moisture[q1_mask, q1_index],
            color=color,
            linewidth=1.65,
            linestyle=(0, (3, 2)),
            alpha=0.88,
            label=f"Q1 {label}",
        )
        axes[1].plot(
            q2_time_min,
            q2.moisture[q2_mask, q2_index],
            color=color,
            linewidth=2.25,
            label=f"Q2 {label}",
        )

    q1_last = int(np.flatnonzero(q1_mask)[-1])
    q2_last = int(np.flatnonzero(q2_mask)[-1])
    axes[0].annotate(
        f"30 min 中心温差\nQ2−Q1 = {q2.temperature_c[q2_last, 0] - q1.temperature_c[q1_last, 0]:+.2f} °C",
        xy=(30.0, q2.temperature_c[q2_last, 0]),
        xytext=(18.2, 29.0),
        arrowprops={"arrowstyle": "->", "color": COLORS["navy"], "lw": 1.0},
        color=COLORS["navy"],
        ha="center",
    )
    axes[1].annotate(
        f"30 min 表面含水率差\nQ2−Q1 = {q2.moisture[q2_last, -1] - q1.moisture[q1_last, -1]:+.3f}",
        xy=(30.0, q2.moisture[q2_last, -1]),
        xytext=(18.1, 1.76),
        arrowprops={"arrowstyle": "->", "color": COLORS["vermilion"], "lw": 1.0},
        color=COLORS["vermilion"],
        ha="center",
    )
    axes[0].set(
        title="温度响应对比",
        xlabel="干燥时间/min",
        ylabel="温度/°C",
        xlim=(0.0, 30.0),
        ylim=(27.0, 38.5),
    )
    axes[1].set(
        title="含水率响应对比",
        xlabel="干燥时间/min",
        ylabel="干基含水率 C/(kg·kg⁻¹)",
        xlim=(0.0, 30.0),
        ylim=(1.35, 2.62),
    )
    for panel, axis in zip(("a", "b"), axes, strict=True):
        polish_axis(axis, grid_axis="y")
        add_panel_label(axis, panel)
        axis.legend(ncol=2, loc="best", handlelength=2.6)
    figure.suptitle("Q1 常物性模型与 Q2 变物性模型的前 30 min 对比", fontsize=14.3)
    return save_figure(figure, output_directory, "q1_q2_30min_comparison", formats)


def create_all_figures(
    q2: Q2Data,
    q1: Q1Data,
    output_directory: Path,
    formats: tuple[str, ...],
) -> list[Path]:
    paths = []
    for plotter in (
        plot_q2_overview,
        plot_radial_histories,
        plot_spatiotemporal_fields,
        plot_radial_profiles,
        plot_cross_section_evolution,
        plot_coupling_trajectories,
        plot_property_evolution,
    ):
        paths.extend(plotter(q2, output_directory, formats))
    paths.extend(plot_q1_q2_comparison(q2, q1, output_directory, formats))
    return paths


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="绘制问题 2 的论文用变物性耦合图")
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--q1-results", type=Path, default=DEFAULT_Q1_RESULTS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_FIGURE_DIR)
    parser.add_argument("--output-dt-s", type=float, default=1.0)
    parser.add_argument(
        "--formats",
        nargs="+",
        default=("png", "pdf", "svg"),
        choices=("png", "pdf", "svg"),
    )
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    configure_publication_style()
    q2 = load_q2_results(arguments.data_dir.resolve(), arguments.output_dt_s)
    q1 = load_q1_results(arguments.q1_results.resolve())
    paths = create_all_figures(
        q2,
        q1,
        arguments.output_dir.resolve(),
        tuple(arguments.formats),
    )
    print(f"Q2 图已生成：{len(paths)} 个文件")
    for path in paths:
        print(path)


if __name__ == "__main__":
    main()
