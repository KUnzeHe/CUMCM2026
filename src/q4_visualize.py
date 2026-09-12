"""绘制问题 4 移动边界热湿耦合模型的论文用可视化图。

本脚本只从未四舍五入的结构化结果读取 Q4 正式解，并在物理坐标图中显式
留白当前药材表面以外的区域。为单独识别收缩效应，首次运行会用相同的
Kirchhoff 离散补算“附录 4 物性 + 固定初始半径”对照，并缓存对照结果。
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
import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm, Normalize

from visualization_style import (
    COLORS,
    RADIAL_COLORS,
    RADIAL_STYLES,
    add_panel_label,
    configure_publication_style,
    polish_axis,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA_DIR = PROJECT_ROOT / "outputs" / "q4" / "data"
DEFAULT_FIGURE_DIR = PROJECT_ROOT / "outputs" / "q4" / "figures" / "fancy"
CONTROL_CACHE_NAME = "q4_fixed_control.npz"
TARGET_MOISTURE = 0.15
INITIAL_RADIUS_CM = 2.0
SELECTED_XI = (0.0, 0.25, 0.50, 0.75, 1.0)
MOISTURE_CMAP = "YlGnBu"
TEMPERATURE_CMAP = "magma_r"
TIME_CMAP = LinearSegmentedColormap.from_list(
    "cumcm_q4_time",
    (COLORS["navy"], COLORS["blue"], COLORS["teal"], COLORS["gold"], COLORS["vermilion"]),
)


@dataclass(frozen=True)
class Q4Data:
    """Q4 正式移动半径解。full_* 字段位于材料坐标 xi 上。"""

    time_s: np.ndarray
    radius_cm: np.ndarray
    physical_radius_cm: np.ndarray
    physical_moisture: np.ndarray
    surface_moisture: np.ndarray
    full_time_s: np.ndarray
    xi: np.ndarray
    full_moisture: np.ndarray
    full_temperature_c: np.ndarray
    full_radius_cm: np.ndarray


@dataclass(frozen=True)
class ControlData:
    """同一数值口径下的固定半径对照。"""

    fixed_time_s: np.ndarray
    fixed_cmax: np.ndarray
    fixed_cbar: np.ndarray
    fixed_full_time_s: np.ndarray
    fixed_full_moisture: np.ndarray
    fixed_xi: np.ndarray
    fixed_tf_s: float
    app3_tf_s: float


def _load_array(data_directory: Path, name: str) -> np.ndarray:
    path = data_directory / name
    if not path.is_file():
        raise FileNotFoundError(
            f"缺少 Q4 结构化结果 {path}。请先运行 "
            "python src/q4_moving_boundary.py --quick；新版文件名用于避免 macOS "
            "大小写不敏感导致的数据覆盖。"
        )
    return np.asarray(np.load(path), dtype=float)


def _prepend_initial(data: Q4Data) -> Q4Data:
    return Q4Data(
        time_s=np.r_[0.0, data.time_s],
        radius_cm=np.r_[INITIAL_RADIUS_CM, data.radius_cm],
        physical_radius_cm=data.physical_radius_cm,
        physical_moisture=np.vstack(
            (np.full(data.physical_radius_cm.size, 2.55), data.physical_moisture)
        ),
        surface_moisture=np.r_[2.55, data.surface_moisture],
        full_time_s=np.r_[0.0, data.full_time_s],
        xi=data.xi,
        full_moisture=np.vstack((np.full(data.xi.size, 2.55), data.full_moisture)),
        full_temperature_c=np.vstack(
            (np.full(data.xi.size, 28.0), data.full_temperature_c)
        ),
        full_radius_cm=np.r_[INITIAL_RADIUS_CM, data.full_radius_cm],
    )


def load_q4_results(data_directory: Path) -> Q4Data:
    """读取并交叉校验新版、无文件名冲突的 Q4 结果。"""

    data = Q4Data(
        time_s=_load_array(data_directory, "Q4_time.npy"),
        radius_cm=100.0 * _load_array(data_directory, "Q4_Rhist.npy"),
        physical_radius_cm=_load_array(data_directory, "Q4_rout.npy"),
        physical_moisture=_load_array(data_directory, "Q4_C.npy"),
        surface_moisture=_load_array(data_directory, "Q4_Csurf.npy"),
        full_time_s=_load_array(data_directory, "Q4_full_time.npy"),
        xi=_load_array(data_directory, "Q4_xi.npy"),
        full_moisture=_load_array(data_directory, "Q4_full_C.npy"),
        full_temperature_c=_load_array(data_directory, "Q4_full_Temp.npy"),
        full_radius_cm=100.0 * _load_array(data_directory, "Q4_full_Rhist.npy"),
    )
    if data.physical_moisture.shape != (
        data.time_s.size,
        data.physical_radius_cm.size,
    ):
        raise ValueError("Q4 固定物理位置含水率矩阵尺寸错误")
    expected_full = (data.full_time_s.size, data.xi.size)
    if data.full_moisture.shape != expected_full or data.full_temperature_c.shape != expected_full:
        raise ValueError("Q4 材料坐标双场矩阵尺寸错误")
    if data.radius_cm.shape != data.time_s.shape or data.surface_moisture.shape != data.time_s.shape:
        raise ValueError("Q4 半径或表面含水率历史长度错误")
    if data.full_radius_cm.shape != data.full_time_s.shape:
        raise ValueError("Q4 完整剖面与半径历史长度不一致")
    if not np.all(np.diff(data.time_s) > 0) or not np.all(np.diff(data.full_time_s) > 0):
        raise ValueError("Q4 时间坐标必须严格递增")
    if not np.all(np.diff(data.xi) > 0) or not np.allclose((data.xi[0], data.xi[-1]), (0.0, 1.0)):
        raise ValueError("Q4 材料坐标必须严格递增并覆盖 [0,1]")
    if np.any(np.diff(data.radius_cm) > 2e-12) or np.any(np.diff(data.full_radius_cm) > 2e-12):
        raise ValueError("Q4 半径历史出现回升")
    if not np.isfinite(data.full_moisture).all() or not np.isfinite(data.full_temperature_c).all():
        raise ValueError("Q4 材料坐标双场含有非有限值")
    if np.nanmin(data.physical_moisture) <= 0 or np.min(data.full_moisture) <= 0:
        raise ValueError("Q4 含水率中出现非正值")
    outside = data.physical_radius_cm[None, :] > data.radius_cm[:, None] * (1.0 + 1e-12)
    if not np.array_equal(np.isnan(data.physical_moisture), outside):
        raise ValueError("Q4 物理位置留白区域与移动半径不一致")
    if np.any(np.diff(np.max(data.full_moisture, axis=1)) > 1e-9):
        raise ValueError("Q4 全域最大含水率随时间反弹")
    return _prepend_initial(data)


def _control_from_npz(path: Path) -> ControlData:
    with np.load(path) as raw:
        return ControlData(
            fixed_time_s=np.asarray(raw["fixed_time_s"], dtype=float),
            fixed_cmax=np.asarray(raw["fixed_cmax"], dtype=float),
            fixed_cbar=np.asarray(raw["fixed_cbar"], dtype=float),
            fixed_full_time_s=np.asarray(raw["fixed_full_time_s"], dtype=float),
            fixed_full_moisture=np.asarray(raw["fixed_full_moisture"], dtype=float),
            fixed_xi=np.asarray(raw["fixed_xi"], dtype=float),
            fixed_tf_s=float(raw["fixed_tf_s"]),
            app3_tf_s=float(raw["app3_tf_s"]),
        )


def load_or_compute_control(data_directory: Path, recompute: bool = False) -> ControlData:
    """补算并缓存固定半径对照；不使用调和平均旧结果。"""

    cache = data_directory / CONTROL_CACHE_NAME
    if cache.is_file() and not recompute:
        return _control_from_npz(cache)

    from q4_moving_boundary import P4, constant_radius, solve_q4

    fixed = solve_q4(
        constant_radius(P4.R0),
        props="app4",
        keep_full_dt=600.0,
        verbose=False,
    )
    app3 = solve_q4(
        constant_radius(P4.R0),
        props="app3",
        keep_full_dt=None,
        verbose=False,
    )
    np.savez_compressed(
        cache,
        fixed_time_s=fixed["t"],
        fixed_cmax=fixed["Cmax"],
        fixed_cbar=fixed["Cbar"],
        fixed_full_time_s=fixed["full_t"],
        fixed_full_moisture=fixed["full_C"],
        fixed_xi=fixed["xi"],
        fixed_tf_s=fixed["t_f"],
        app3_tf_s=app3["t_f"],
    )
    return _control_from_npz(cache)


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


def nearest_indices(values: np.ndarray, targets: tuple[float, ...]) -> list[int]:
    return [int(np.argmin(np.abs(values - target))) for target in targets]


def material_label(xi: float) -> str:
    if np.isclose(xi, 0.0):
        return "中心 ξ = 0"
    if np.isclose(xi, 1.0):
        return "表面 ξ = 1"
    return f"ξ = {xi:g}"


def material_style(xi: float) -> tuple[str, object]:
    slot = int(round(xi * (len(RADIAL_COLORS) - 1)))
    slot = min(max(slot, 0), len(RADIAL_COLORS) - 1)
    return RADIAL_COLORS[slot], RADIAL_STYLES[slot]


def material_weights(xi: np.ndarray) -> np.ndarray:
    dxi = float(xi[1] - xi[0])
    weights = xi * dxi
    weights[0] = dxi**2 / 8.0
    weights[-1] = dxi / 2.0 * (1.0 - dxi / 4.0)
    if not np.isclose(weights.sum(), 0.5, atol=1e-12):
        raise ValueError("材料面积权重不守恒")
    return weights


def material_mean(moisture: np.ndarray, xi: np.ndarray) -> np.ndarray:
    return 2.0 * (moisture * material_weights(xi)[None, :]).sum(axis=1)


def moving_event_time_s(data: Q4Data) -> float:
    """在 60 s 未舍入结果间线性定位 Cmax=0.15。"""

    maximum = np.nanmax(
        np.column_stack((data.physical_moisture, data.surface_moisture)), axis=1
    )
    hit = np.flatnonzero(maximum < TARGET_MOISTURE)
    if hit.size == 0:
        raise ValueError("Q4 正式结果未达到终止阈值")
    index = int(hit[0])
    if index == 0:
        return float(data.time_s[0])
    t0, t1 = data.time_s[index - 1 : index + 1]
    c0, c1 = maximum[index - 1 : index + 1]
    return float(t0 + (TARGET_MOISTURE - c0) / (c1 - c0) * (t1 - t0))


def physical_field(
    time_s: np.ndarray,
    xi: np.ndarray,
    radius_cm: np.ndarray,
    material_values: np.ndarray,
    radial_points: int = 241,
) -> tuple[np.ndarray, np.ndarray, np.ma.MaskedArray]:
    """将材料坐标剖面映射到公共物理半径网格，域外保持掩膜。"""

    radial_grid = np.linspace(0.0, INITIAL_RADIUS_CM, radial_points)
    values = np.full((time_s.size, radial_points), np.nan)
    for index, radius in enumerate(radius_cm):
        inside = radial_grid <= radius * (1.0 + 1e-12)
        query_xi = np.minimum(radial_grid[inside] / radius, 1.0)
        values[index, inside] = np.interp(query_xi, xi, material_values[index])
    return time_s / 3600.0, radial_grid, np.ma.masked_invalid(values.T)


def _boundary_line(axis, time_h: np.ndarray, radius_cm: np.ndarray) -> None:
    line = axis.plot(
        time_h,
        radius_cm,
        color=COLORS["charcoal"],
        linewidth=1.7,
        label="移动表面 R(t)",
        zorder=5,
    )[0]
    line.set_path_effects([path_effects.Stroke(linewidth=3.5, foreground="white"), path_effects.Normal()])


def draw_moving_field(
    axis,
    time_h: np.ndarray,
    radial_grid: np.ndarray,
    field: np.ma.MaskedArray,
    radius_cm: np.ndarray,
    *,
    cmap: str,
    norm,
    colorbar_label: str,
):
    palette = plt.get_cmap(cmap).copy()
    palette.set_bad("white")
    mesh = axis.pcolormesh(
        time_h,
        radial_grid,
        field,
        shading="auto",
        cmap=palette,
        norm=norm,
        rasterized=True,
    )
    _boundary_line(axis, time_h, radius_cm)
    axis.fill_between(
        time_h,
        radius_cm,
        INITIAL_RADIUS_CM,
        color="white",
        alpha=0.18,
        hatch="///",
        edgecolor=COLORS["light_gray"],
        linewidth=0.0,
        zorder=4,
    )
    axis.set(
        xlabel="干燥时间/h",
        ylabel="物理径向位置 r/cm",
        xlim=(time_h[0], time_h[-1]),
        ylim=(0.0, INITIAL_RADIUS_CM),
    )
    axis.set_yticks([0.0, 0.5, 1.0, 1.5, 2.0], ["中心", "0.5", "1.0", "1.5", "初始表面"])
    colorbar = axis.figure.colorbar(mesh, ax=axis, pad=0.02, fraction=0.035)
    colorbar.set_label(colorbar_label)
    return mesh


def _duration_axis(axis, data: Q4Data, control: ControlData, *, title: str) -> None:
    moving_h = moving_event_time_s(data) / 3600.0
    values = np.array([control.app3_tf_s, control.fixed_tf_s, moving_event_time_s(data)]) / 3600.0
    labels = ("附录 3\n固定半径", "附录 4\n固定半径", "附录 4\n移动半径")
    colors = (COLORS["gray"], COLORS["gold"], COLORS["teal"])
    positions = np.arange(3)
    axis.vlines(positions, 0.0, values, color=colors, linewidth=8.0, alpha=0.28)
    axis.scatter(positions, values, s=105, color=colors, edgecolor="white", linewidth=1.2, zorder=4)
    for x, value, color in zip(positions, values, colors, strict=True):
        axis.text(x, value + 5.0, f"{value:.2f} h", ha="center", color=color, fontweight="semibold")
    reduction = (control.fixed_tf_s / 3600.0 - moving_h) / (control.fixed_tf_s / 3600.0) * 100.0
    axis.annotate(
        f"同物性下缩短 {reduction:.1f}%",
        xy=(2.0, moving_h),
        xytext=(1.55, 93.0),
        arrowprops={"arrowstyle": "->", "color": COLORS["teal"], "lw": 1.2},
        ha="center",
        color=COLORS["teal"],
    )
    axis.set(title=title, ylabel="连续终止时间/h", xticks=positions, xticklabels=labels, ylim=(0.0, 150.0))
    polish_axis(axis, grid_axis="y")


def plot_q4_overview(
    data: Q4Data,
    control: ControlData,
    output_directory: Path,
    formats: tuple[str, ...],
) -> list[Path]:
    """绘制适合论文和答辩首页的四面板总览。"""

    figure, axes = plt.subplots(2, 2, figsize=(13.2, 9.1), constrained_layout=True)
    time_h = data.time_s / 3600.0
    radius_ratio = data.radius_cm / INITIAL_RADIUS_CM

    axes[0, 0].plot(time_h, data.radius_cm, color=COLORS["navy"], linewidth=2.5)
    axes[0, 0].fill_between(time_h, data.radius_cm, INITIAL_RADIUS_CM, color=COLORS["teal"], alpha=0.12)
    axes[0, 0].axvspan(0.0, 18.0, color=COLORS["gold"], alpha=0.09)
    axes[0, 0].scatter([0.0, time_h[-1]], [INITIAL_RADIUS_CM, data.radius_cm[-1]],
                       color=(COLORS["navy"], COLORS["vermilion"]), s=58, zorder=4)
    axes[0, 0].annotate("前 18 h 为主要收缩阶段", xy=(12.0, 1.25), xytext=(21.0, 1.57),
                        arrowprops={"arrowstyle": "->", "color": COLORS["gold"], "lw": 1.15},
                        color=COLORS["charcoal"])
    axes[0, 0].set(title="给定半径历史", xlabel="干燥时间/h", ylabel="当前半径 R(t)/cm",
                   xlim=(0.0, time_h[-1]), ylim=(1.12, 2.06))
    polish_axis(axes[0, 0], grid_axis="y")
    axes[0, 0].text(0.98, 0.93, f"截面积降至 {radius_ratio[-1] ** 2 * 100:.1f}%",
                    transform=axes[0, 0].transAxes, ha="right", va="top", color=COLORS["vermilion"])
    add_panel_label(axes[0, 0], "a")

    full_time_h = data.full_time_s / 3600.0
    maximum = np.max(data.full_moisture, axis=1)
    axes[0, 1].plot(full_time_h, maximum, color=COLORS["navy"], linewidth=2.4, label="全域最大值")
    axes[0, 1].plot(full_time_h, data.full_moisture[:, -1], color=COLORS["vermilion"],
                    linewidth=1.9, label="表面")
    axes[0, 1].axhline(TARGET_MOISTURE, color=COLORS["vermilion"], linestyle="--", linewidth=1.3,
                       label="终止阈值 0.15")
    tf_h = moving_event_time_s(data) / 3600.0
    axes[0, 1].axvline(tf_h, color=COLORS["gold"], linestyle=(0, (3, 2)), linewidth=1.4)
    axes[0, 1].annotate(f"连续终止 {tf_h:.3f} h", xy=(tf_h, TARGET_MOISTURE), xytext=(31.0, 0.31),
                        arrowprops={"arrowstyle": "->", "color": COLORS["gold"], "lw": 1.15},
                        color=COLORS["charcoal"])
    axes[0, 1].set(title="中心控制最终达标", xlabel="干燥时间/h",
                   ylabel="干基含水率 C/(kg·kg⁻¹)", xlim=(0.0, time_h[-1]),
                   yscale="log", ylim=(0.045, 3.1))
    axes[0, 1].set_yticks([0.05, 0.10, 0.15, 0.30, 0.60, 1.20, 2.40])
    axes[0, 1].set_yticklabels(["0.05", "0.10", "0.15", "0.30", "0.60", "1.20", "2.40"])
    polish_axis(axes[0, 1], grid_axis="y")
    axes[0, 1].legend(loc="upper right")
    add_panel_label(axes[0, 1], "b")

    map_time_h, radial_grid, moisture_field = physical_field(
        data.full_time_s, data.xi, data.full_radius_cm, data.full_moisture
    )
    mesh = draw_moving_field(
        axes[1, 0], map_time_h, radial_grid, moisture_field, data.full_radius_cm,
        cmap=MOISTURE_CMAP,
        norm=LogNorm(vmin=0.05, vmax=2.55),
        colorbar_label="C/(kg·kg⁻¹)",
    )
    axes[1, 0].contour(map_time_h, radial_grid, moisture_field, levels=[TARGET_MOISTURE],
                       colors=[COLORS["vermilion"]], linewidths=1.35)
    axes[1, 0].set_title("物理域含水率：边界以内才有材料")
    axes[1, 0].text(0.98, 0.91, "斜线区：收缩后已在药材外部", transform=axes[1, 0].transAxes,
                    ha="right", color=COLORS["gray"], fontsize=8.9,
                    bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "none", "alpha": 0.86})
    add_panel_label(axes[1, 0], "c")
    _ = mesh

    _duration_axis(axes[1, 1], data, control, title="三组同离散模型拆分物性与收缩效应")
    add_panel_label(axes[1, 1], "d")

    figure.suptitle("问题 4｜移动边界下的收缩—热湿耦合", fontsize=15.3)
    return save_figure(figure, output_directory, "q4_overview", formats)


def plot_moving_domain_fields(
    data: Q4Data,
    control: ControlData,
    output_directory: Path,
    formats: tuple[str, ...],
) -> list[Path]:
    del control
    time_h, radial_grid, temperature = physical_field(
        data.full_time_s, data.xi, data.full_radius_cm, data.full_temperature_c
    )
    _, _, moisture = physical_field(
        data.full_time_s, data.xi, data.full_radius_cm, data.full_moisture
    )
    figure, axes = plt.subplots(2, 1, figsize=(11.5, 8.2), constrained_layout=True)
    draw_moving_field(
        axes[0], time_h, radial_grid, temperature, data.full_radius_cm,
        cmap=TEMPERATURE_CMAP,
        norm=Normalize(vmin=28.0, vmax=50.05),
        colorbar_label="温度 T/°C",
    )
    axes[0].set_title("温度场：前几小时完成预热，几何边界继续收缩")
    add_panel_label(axes[0], "a")
    draw_moving_field(
        axes[1], time_h, radial_grid, moisture, data.full_radius_cm,
        cmap=MOISTURE_CMAP,
        norm=LogNorm(vmin=0.05, vmax=2.55),
        colorbar_label="干基含水率 C/(kg·kg⁻¹)，对数色标",
    )
    contour = axes[1].contour(time_h, radial_grid, moisture, levels=[TARGET_MOISTURE],
                              colors=[COLORS["vermilion"]], linewidths=1.55)
    if contour.allsegs and contour.allsegs[0]:
        axes[1].clabel(contour, fmt={TARGET_MOISTURE: "C = 0.15"}, fontsize=8.8)
    axes[1].set_title("含水率场：低含水层由移动表面向中心推进")
    add_panel_label(axes[1], "b")
    figure.suptitle("Q4 移动物理域中的温度—含水率双场", fontsize=14.6)
    return save_figure(figure, output_directory, "q4_moving_domain_fields", formats)


def plot_radius_transport(
    data: Q4Data,
    control: ControlData,
    output_directory: Path,
    formats: tuple[str, ...],
) -> list[Path]:
    del control
    time_h = data.time_s / 3600.0
    ratio = data.radius_cm / INITIAL_RADIUS_CM
    figure, axes = plt.subplots(1, 2, figsize=(12.4, 5.15), constrained_layout=True)

    axes[0].plot(time_h, ratio, color=COLORS["navy"], linewidth=2.4, label="半径比 R/R₀")
    axes[0].plot(time_h, ratio**2, color=COLORS["teal"], linewidth=2.4, label="截面积比 (R/R₀)²")
    axes[0].fill_between(time_h, ratio**2, ratio, color=COLORS["teal"], alpha=0.08)
    axes[0].axvspan(0.0, 18.0, color=COLORS["gold"], alpha=0.09)
    i18 = int(np.argmin(np.abs(time_h - 18.0)))
    axes[0].scatter([time_h[i18]], [ratio[i18]], color=COLORS["gold"], s=58, zorder=4)
    axes[0].annotate(f"18 h：R = {data.radius_cm[i18]:.3f} cm", xy=(time_h[i18], ratio[i18]),
                     xytext=(24.0, 0.74), arrowprops={"arrowstyle": "->", "color": COLORS["gold"], "lw": 1.1},
                     color=COLORS["charcoal"])
    axes[0].set(title="几何尺度压缩", xlabel="干燥时间/h", ylabel="相对初始值", xlim=(0.0, time_h[-1]), ylim=(0.30, 1.04))
    polish_axis(axes[0], grid_axis="y")
    axes[0].legend(loc="upper right")
    add_panel_label(axes[0], "a")

    internal = ratio**-2
    boundary = ratio**-1
    axes[1].plot(time_h, internal, color=COLORS["vermilion"], linewidth=2.5,
                 label="内部传输尺度 (R₀/R)²")
    axes[1].plot(time_h, boundary, color=COLORS["blue"], linewidth=2.3,
                 label="表面项尺度 R₀/R")
    axes[1].fill_between(time_h, 1.0, internal, color=COLORS["vermilion"], alpha=0.09)
    axes[1].scatter([time_h[-1], time_h[-1]], [internal[-1], boundary[-1]],
                    color=(COLORS["vermilion"], COLORS["blue"]), s=55, zorder=4)
    axes[1].annotate(f"内部项放大 {internal[-1]:.2f} 倍", xy=(time_h[-1], internal[-1]),
                     xytext=(29.0, 2.42), arrowprops={"arrowstyle": "->", "color": COLORS["vermilion"], "lw": 1.1},
                     color=COLORS["vermilion"])
    axes[1].set(title="收缩进入方程的两个通道", xlabel="干燥时间/h", ylabel="相对初始半径的尺度因子",
                xlim=(0.0, time_h[-1]), ylim=(0.92, 3.02))
    polish_axis(axes[1], grid_axis="y")
    axes[1].legend(loc="center right")
    add_panel_label(axes[1], "b")
    figure.suptitle("Q4 收缩机制：更短路径与更大的表面积—体积比共同加速传输", fontsize=14.4)
    return save_figure(figure, output_directory, "q4_radius_transport", formats)


def plot_material_histories(
    data: Q4Data,
    control: ControlData,
    output_directory: Path,
    formats: tuple[str, ...],
) -> list[Path]:
    del control
    indices = nearest_indices(data.xi, SELECTED_XI)
    time_h = data.full_time_s / 3600.0
    figure, axes = plt.subplots(1, 2, figsize=(12.6, 5.2), constrained_layout=True)
    for index in indices:
        color, style = material_style(data.xi[index])
        axes[0].plot(time_h, data.full_moisture[:, index], color=color, linestyle=style,
                     linewidth=2.05, label=material_label(data.xi[index]))
        axes[1].plot(time_h, data.full_temperature_c[:, index], color=color, linestyle=style,
                     linewidth=2.05, label=material_label(data.xi[index]))
    axes[0].axhline(TARGET_MOISTURE, color=COLORS["vermilion"], linestyle="--", linewidth=1.3)
    axes[0].set(title="材料层含水率响应", xlabel="干燥时间/h", ylabel="干基含水率 C/(kg·kg⁻¹)",
                xlim=(0.0, time_h[-1]), yscale="log", ylim=(0.045, 3.1))
    axes[0].set_yticks([0.05, 0.10, 0.15, 0.30, 0.60, 1.20, 2.40])
    axes[0].set_yticklabels(["0.05", "0.10", "0.15", "0.30", "0.60", "1.20", "2.40"])
    polish_axis(axes[0], grid_axis="y")
    axes[0].legend(ncol=2, loc="upper right", handlelength=2.6)
    add_panel_label(axes[0], "a")

    axes[1].axhline(50.0049, color=COLORS["gray"], linestyle="--", linewidth=1.2, label="平台环境温度")
    axes[1].set(title="材料层温度响应（前 6 h 放大）", xlabel="干燥时间/h", ylabel="温度 T/°C",
                xlim=(0.0, 6.0), ylim=(27.3, 50.8))
    polish_axis(axes[1], grid_axis="y")
    axes[1].legend(ncol=2, loc="lower right", handlelength=2.5)
    add_panel_label(axes[1], "b")
    figure.suptitle("Q4 固定材料坐标 ξ 上的双场历史", fontsize=14.4)
    return save_figure(figure, output_directory, "q4_material_histories", formats)


def plot_physical_profiles(
    data: Q4Data,
    control: ControlData,
    output_directory: Path,
    formats: tuple[str, ...],
) -> list[Path]:
    del control
    moisture_times = (0.0, 6.0, 12.0, 24.0, 36.0, 48.0, 51.0)
    temperature_times = (0.0, 0.5, 1.0, 2.0, 4.0, 6.0)
    m_indices = nearest_indices(data.full_time_s / 3600.0, moisture_times)
    t_indices = nearest_indices(data.full_time_s / 3600.0, temperature_times)
    figure, axes = plt.subplots(1, 2, figsize=(12.7, 5.25), constrained_layout=True)

    for index, color in zip(m_indices, TIME_CMAP(np.linspace(0.02, 0.98, len(m_indices))), strict=True):
        r = data.xi * data.full_radius_cm[index]
        axes[0].plot(r, data.full_moisture[index], color=color, linewidth=2.0,
                     label=f"{data.full_time_s[index] / 3600:g} h")
        axes[0].scatter([r[-1]], [data.full_moisture[index, -1]], s=18, color=[color], zorder=4)
    axes[0].axhline(TARGET_MOISTURE, color=COLORS["vermilion"], linestyle="--", linewidth=1.25)
    axes[0].axvline(INITIAL_RADIUS_CM, color=COLORS["light_gray"], linestyle=":", linewidth=1.2)
    axes[0].set(title="含水率剖面：曲线右端就是当前表面", xlabel="物理径向位置 r/cm",
                ylabel="干基含水率 C/(kg·kg⁻¹)", xlim=(0.0, 2.05), yscale="log", ylim=(0.045, 3.1))
    axes[0].set_yticks([0.05, 0.10, 0.15, 0.30, 0.60, 1.20, 2.40])
    axes[0].set_yticklabels(["0.05", "0.10", "0.15", "0.30", "0.60", "1.20", "2.40"])
    polish_axis(axes[0], grid_axis="both")
    axes[0].legend(title="干燥时间", ncol=2, loc="upper right")
    add_panel_label(axes[0], "a")

    for index, color in zip(t_indices, TIME_CMAP(np.linspace(0.02, 0.88, len(t_indices))), strict=True):
        r = data.xi * data.full_radius_cm[index]
        axes[1].plot(r, data.full_temperature_c[index], color=color, linewidth=2.05,
                     label=f"{data.full_time_s[index] / 3600:g} h")
        axes[1].scatter([r[-1]], [data.full_temperature_c[index, -1]], s=18, color=[color], zorder=4)
    axes[1].axvline(INITIAL_RADIUS_CM, color=COLORS["light_gray"], linestyle=":", linewidth=1.2)
    axes[1].set(title="温度剖面：预热与收缩同步发生", xlabel="物理径向位置 r/cm", ylabel="温度 T/°C",
                xlim=(0.0, 2.05), ylim=(27.0, 51.0))
    polish_axis(axes[1], grid_axis="both")
    axes[1].legend(title="干燥时间", ncol=2, loc="lower left")
    add_panel_label(axes[1], "b")
    figure.suptitle("Q4 映射回物理半径后的径向剖面", fontsize=14.4)
    return save_figure(figure, output_directory, "q4_physical_profiles", formats)


def _circular_material_field(xi_grid: np.ndarray, values: np.ndarray, radius_cm: float):
    coordinates = np.linspace(-INITIAL_RADIUS_CM, INITIAL_RADIUS_CM, 301)
    xx, yy = np.meshgrid(coordinates, coordinates)
    rr = np.sqrt(xx**2 + yy**2)
    query = np.minimum(rr / radius_cm, 1.0)
    field = np.interp(query, xi_grid, values)
    return np.ma.masked_where(rr > radius_cm, field)


def plot_cross_section_evolution(
    data: Q4Data,
    control: ControlData,
    output_directory: Path,
    formats: tuple[str, ...],
) -> list[Path]:
    del control
    times_h = (0.0, 3.0, 6.0, 12.0, 24.0, 51.0)
    indices = nearest_indices(data.full_time_s / 3600.0, times_h)
    moisture_norm = LogNorm(vmin=0.05, vmax=2.55)
    temperature_norm = Normalize(vmin=28.0, vmax=50.05)
    moisture_cmap = plt.get_cmap(MOISTURE_CMAP).copy()
    temperature_cmap = plt.get_cmap(TEMPERATURE_CMAP).copy()
    moisture_cmap.set_bad("white", alpha=0.0)
    temperature_cmap.set_bad("white", alpha=0.0)
    figure, axes = plt.subplots(2, len(indices), figsize=(15.2, 5.65), constrained_layout=True)

    for column, index in enumerate(indices):
        radius = data.full_radius_cm[index]
        temp_image = axes[0, column].imshow(
            _circular_material_field(data.xi, data.full_temperature_c[index], radius),
            extent=(-INITIAL_RADIUS_CM, INITIAL_RADIUS_CM, -INITIAL_RADIUS_CM, INITIAL_RADIUS_CM),
            origin="lower", cmap=temperature_cmap, norm=temperature_norm, interpolation="bilinear",
        )
        moist_image = axes[1, column].imshow(
            _circular_material_field(data.xi, data.full_moisture[index], radius),
            extent=(-INITIAL_RADIUS_CM, INITIAL_RADIUS_CM, -INITIAL_RADIUS_CM, INITIAL_RADIUS_CM),
            origin="lower", cmap=moisture_cmap, norm=moisture_norm, interpolation="bilinear",
        )
        for row in range(2):
            axes[row, column].add_patch(plt.Circle((0, 0), INITIAL_RADIUS_CM, fill=False,
                                                   color=COLORS["light_gray"], linestyle=":", linewidth=0.8))
            axes[row, column].add_patch(plt.Circle((0, 0), radius, fill=False,
                                                   color=COLORS["charcoal"], linewidth=0.85))
            axes[row, column].set(xlim=(-2.08, 2.08), ylim=(-2.08, 2.08), aspect="equal")
            axes[row, column].axis("off")
        axes[0, column].set_title(f"{data.full_time_s[index] / 3600:g} h\nR = {radius:.3f} cm", fontsize=10.2)

    tbar = figure.colorbar(temp_image, ax=axes[0, :].tolist(), pad=0.012, fraction=0.025)
    tbar.set_label("温度 T/°C")
    cbar = figure.colorbar(moist_image, ax=axes[1, :].tolist(), pad=0.012, fraction=0.025)
    cbar.set_label("干基含水率 C/(kg·kg⁻¹)")
    cbar.set_ticks([0.05, 0.10, 0.15, 0.30, 0.60, 1.20, 2.40])
    cbar.set_ticklabels(["0.05", "0.10", "0.15", "0.30", "0.60", "1.20", "2.40"])
    figure.suptitle("Q4 真实尺寸圆截面演化：颜色变化与几何收缩同步可见", fontsize=14.5)
    figure.text(0.5, 0.005, "灰色虚线为初始半径 2 cm；黑色实线为各时刻实际表面",
                ha="center", fontsize=9.0, color=COLORS["gray"])
    return save_figure(figure, output_directory, "q4_cross_section_evolution", formats)


def plot_shrinkage_comparison(
    data: Q4Data,
    control: ControlData,
    output_directory: Path,
    formats: tuple[str, ...],
) -> list[Path]:
    moving_time_h = data.full_time_s / 3600.0
    moving_cmax = np.max(data.full_moisture, axis=1)
    moving_cbar = material_mean(data.full_moisture, data.xi)
    fixed_time_h = control.fixed_time_s / 3600.0
    figure, axes = plt.subplots(2, 2, figsize=(13.0, 8.7), constrained_layout=True)

    axes[0, 0].plot(fixed_time_h, control.fixed_cmax, color=COLORS["gold"], linewidth=2.3,
                    label="固定半径 R = 2 cm")
    axes[0, 0].plot(moving_time_h, moving_cmax, color=COLORS["teal"], linewidth=2.5,
                    label="移动半径 R(t)")
    axes[0, 0].axhline(TARGET_MOISTURE, color=COLORS["vermilion"], linestyle="--", linewidth=1.3)
    axes[0, 0].set(title="全域最大含水率", xlabel="干燥时间/h", ylabel="Cmax/(kg·kg⁻¹)",
                   xlim=(0.0, fixed_time_h[-1]), yscale="log", ylim=(0.045, 3.1))
    axes[0, 0].set_yticks([0.05, 0.10, 0.15, 0.30, 0.60, 1.20, 2.40])
    axes[0, 0].set_yticklabels(["0.05", "0.10", "0.15", "0.30", "0.60", "1.20", "2.40"])
    polish_axis(axes[0, 0], grid_axis="y")
    axes[0, 0].legend(loc="upper right")
    add_panel_label(axes[0, 0], "a")

    axes[0, 1].plot(control.fixed_time_s / 3600.0, control.fixed_cbar, color=COLORS["gold"], linewidth=2.3,
                    label="固定半径")
    axes[0, 1].plot(moving_time_h, moving_cbar, color=COLORS["teal"], linewidth=2.5,
                    label="移动半径")
    axes[0, 1].set(title="材料面积加权平均含水率", xlabel="干燥时间/h", ylabel="材料平均含水率 C̄/(kg·kg⁻¹)",
                   xlim=(0.0, fixed_time_h[-1]), yscale="log", ylim=(0.045, 3.1))
    axes[0, 1].set_yticks([0.05, 0.10, 0.15, 0.30, 0.60, 1.20, 2.40])
    axes[0, 1].set_yticklabels(["0.05", "0.10", "0.15", "0.30", "0.60", "1.20", "2.40"])
    polish_axis(axes[0, 1], grid_axis="y")
    axes[0, 1].legend(loc="upper right")
    add_panel_label(axes[0, 1], "b")

    target_h = 48.0
    moving_index = nearest_indices(moving_time_h, (target_h,))[0]
    fixed_full_time_h = control.fixed_full_time_s / 3600.0
    fixed_index = nearest_indices(fixed_full_time_h, (target_h,))[0]
    axes[1, 0].plot(control.fixed_xi * INITIAL_RADIUS_CM, control.fixed_full_moisture[fixed_index],
                    color=COLORS["gold"], linewidth=2.4, label="固定半径")
    axes[1, 0].plot(data.xi * data.full_radius_cm[moving_index], data.full_moisture[moving_index],
                    color=COLORS["teal"], linewidth=2.5, label="移动半径")
    axes[1, 0].scatter([INITIAL_RADIUS_CM, data.full_radius_cm[moving_index]],
                       [control.fixed_full_moisture[fixed_index, -1], data.full_moisture[moving_index, -1]],
                       color=(COLORS["gold"], COLORS["teal"]), s=48, zorder=4)
    axes[1, 0].set(title="同一时刻的物理剖面对照（48 h）", xlabel="物理径向位置 r/cm",
                   ylabel="干基含水率 C/(kg·kg⁻¹)", xlim=(0.0, 2.05), ylim=(0.035, 0.44))
    polish_axis(axes[1, 0], grid_axis="both")
    axes[1, 0].legend(loc="upper right")
    add_panel_label(axes[1, 0], "c")

    _duration_axis(axes[1, 1], data, control, title="烘干时长分解")
    add_panel_label(axes[1, 1], "d")
    figure.suptitle("Q4 收缩效应：必须在相同附录 4 物性下比较", fontsize=14.6)
    return save_figure(figure, output_directory, "q4_shrinkage_comparison", formats)


def plot_geometry_state_portrait(
    data: Q4Data,
    control: ControlData,
    output_directory: Path,
    formats: tuple[str, ...],
) -> list[Path]:
    del control
    time_h = data.full_time_s / 3600.0
    ratio = data.full_radius_cm / INITIAL_RADIUS_CM
    shrinkage_pct = 100.0 * (1.0 - ratio)
    cbar = material_mean(data.full_moisture, data.xi)
    gradient = data.full_moisture[:, 0] - data.full_moisture[:, -1]
    acceleration = ratio**-2
    sample = np.unique(np.linspace(0, time_h.size - 1, 260).astype(int))
    norm = Normalize(vmin=0.0, vmax=float(time_h[-1]))
    figure, axes = plt.subplots(1, 2, figsize=(12.2, 5.55), constrained_layout=True)

    axes[0].plot(shrinkage_pct, cbar, color=COLORS["light_gray"], linewidth=3.3, zorder=1)
    scatter = axes[0].scatter(shrinkage_pct[sample], cbar[sample], c=time_h[sample], cmap=TIME_CMAP,
                              norm=norm, s=30, linewidth=0.0, zorder=3)
    for hour in (6.0, 12.0, 24.0, 48.0):
        index = nearest_indices(time_h, (hour,))[0]
        axes[0].annotate(f"{hour:g} h", xy=(shrinkage_pct[index], cbar[index]), xytext=(5, 5),
                         textcoords="offset points", fontsize=8.7, color=COLORS["charcoal"])
    axes[0].set(title="几何—含水状态轨迹", xlabel="半径收缩率/%", ylabel="材料平均含水率 C̄/(kg·kg⁻¹)",
                xlim=(-0.8, max(shrinkage_pct) + 1.6), yscale="log", ylim=(0.06, 2.85))
    polish_axis(axes[0], grid_axis="both")
    add_panel_label(axes[0], "a")

    axes[1].plot(acceleration, gradient, color=COLORS["light_gray"], linewidth=3.3, zorder=1)
    axes[1].scatter(acceleration[sample], gradient[sample], c=time_h[sample], cmap=TIME_CMAP,
                    norm=norm, s=30, linewidth=0.0, zorder=3)
    peak = int(np.argmax(gradient))
    axes[1].scatter([acceleration[peak]], [gradient[peak]], color=COLORS["vermilion"], s=70,
                    edgecolor="white", linewidth=1.0, zorder=5)
    axes[1].annotate(f"径向差最大\n{time_h[peak]:.1f} h", xy=(acceleration[peak], gradient[peak]),
                     xytext=(1.58, 1.72), arrowprops={"arrowstyle": "->", "color": COLORS["vermilion"], "lw": 1.1},
                     color=COLORS["vermilion"])
    axes[1].set(title="传输加速—径向滞后轨迹", xlabel="内部传输尺度因子 (R₀/R)²",
                ylabel="中心—表面含水率差 ΔC/(kg·kg⁻¹)", xlim=(0.96, max(acceleration) + 0.08), ylim=(-0.03, 2.08))
    polish_axis(axes[1], grid_axis="both")
    add_panel_label(axes[1], "b")
    colorbar = figure.colorbar(scatter, ax=axes.ravel().tolist(), pad=0.025, fraction=0.035)
    colorbar.set_label("干燥时间/h")
    figure.suptitle("Q4 状态轨迹：把收缩、传输加速和水分梯度放在同一张图上", fontsize=14.3)
    return save_figure(figure, output_directory, "q4_geometry_state_portrait", formats)


def create_all_figures(
    data: Q4Data,
    control: ControlData,
    output_directory: Path,
    formats: tuple[str, ...],
) -> list[Path]:
    paths = []
    for plotter in (
        plot_q4_overview,
        plot_moving_domain_fields,
        plot_radius_transport,
        plot_material_histories,
        plot_physical_profiles,
        plot_cross_section_evolution,
        plot_shrinkage_comparison,
        plot_geometry_state_portrait,
    ):
        paths.extend(plotter(data, control, output_directory, formats))
    return paths


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="绘制问题 4 的移动边界论文图")
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_FIGURE_DIR)
    parser.add_argument("--recompute-control", action="store_true", help="忽略缓存，重新计算固定半径对照")
    parser.add_argument(
        "--formats", nargs="+", default=("png", "pdf", "svg"), choices=("png", "pdf", "svg")
    )
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    configure_publication_style()
    data_directory = arguments.data_dir.resolve()
    data = load_q4_results(data_directory)
    control = load_or_compute_control(data_directory, recompute=arguments.recompute_control)
    paths = create_all_figures(data, control, arguments.output_dir.resolve(), tuple(arguments.formats))
    print(
        f"Q4 图已生成：{len(paths)} 个文件；移动半径连续终止时间 "
        f"{moving_event_time_s(data) / 3600.0:.3f} h，固定半径对照 {control.fixed_tf_s / 3600.0:.3f} h"
    )
    for path in paths:
        print(path)


if __name__ == "__main__":
    main()
