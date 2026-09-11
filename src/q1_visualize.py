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
from matplotlib import font_manager

from q1_fvm import DEFAULT_ENV_PATH, PROJECT_ROOT, load_environment


DEFAULT_RESULTS_PATH = PROJECT_ROOT / "outputs" / "q1" / "data" / "q1_results_full.npz"
DEFAULT_FIGURE_DIR = PROJECT_ROOT / "outputs" / "q1" / "figures"
DEFAULT_RADII_CM = (0.0, 0.5, 1.0, 1.5, 2.0)
DEFAULT_PROFILE_TIMES_S = (300, 600, 900, 1200, 1500, 1800)
SERIES_COLORS = ("#173F5F", "#20639B", "#3CAEA3", "#F6A623", "#C44536")
SERIES_STYLES = ("-", "--", "-.", ":", (0, (5, 1)))


def configure_style() -> None:
    """设置适合中文论文图的字体、颜色和尺寸。"""

    installed = {font.name for font in font_manager.fontManager.ttflist}
    for family in (
        "Noto Sans CJK SC",
        "Source Han Sans SC",
        "PingFang SC",
        "Microsoft YaHei",
        "SimHei",
        "Arial Unicode MS",
    ):
        if family in installed:
            plt.rcParams["font.sans-serif"] = [family, "DejaVu Sans"]
            break
    plt.rcParams.update(
        {
            "axes.unicode_minus": False,
            "axes.edgecolor": "#333333",
            "axes.labelcolor": "#222222",
            "axes.titleweight": "semibold",
            "axes.titlesize": 13,
            "axes.labelsize": 11,
            "xtick.color": "#444444",
            "ytick.color": "#444444",
            "font.size": 10,
            "legend.frameon": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


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
        figure.savefig(path, dpi=300, bbox_inches="tight")
        paths.append(path.resolve())
    plt.close(figure)
    return paths


def plot_temperature_history(
    time_s: np.ndarray,
    radius_cm: np.ndarray,
    temperature_c: np.ndarray,
    ambient_temperature,
    selected_radii_cm: tuple[float, ...],
) -> plt.Figure:
    figure, axis = plt.subplots(figsize=(8.4, 5.2), constrained_layout=True)
    time_min = time_s / 60.0
    indices = nearest_indices(radius_cm, selected_radii_cm)
    for index, color, style in zip(
        indices,
        SERIES_COLORS[: len(indices)],
        SERIES_STYLES[: len(indices)],
        strict=True,
    ):
        actual_radius = radius_cm[index]
        axis.plot(
            time_min,
            temperature_c[:, index],
            color=color,
            linestyle=style,
            linewidth=2.0,
            label=f"r = {actual_radius:g} cm",
        )
    axis.plot(
        time_min,
        np.asarray(ambient_temperature(time_s), dtype=float),
        color="#222222",
        linestyle=(0, (3, 2)),
        linewidth=1.6,
        alpha=0.8,
        label="烘房空气",
    )
    axis.set(
        title="不同径向位置的温度随时间变化",
        xlabel="时间（min）",
        ylabel="温度（°C）",
        xlim=(0.0, time_min[-1]),
    )
    axis.grid(True, color="#D9DEE5", linewidth=0.7, alpha=0.8)
    axis.legend(ncol=2, loc="upper left")
    return figure


def plot_moisture_history(
    time_s: np.ndarray,
    radius_cm: np.ndarray,
    moisture_kg_kg: np.ndarray,
    selected_radii_cm: tuple[float, ...],
) -> plt.Figure:
    figure, axis = plt.subplots(figsize=(8.4, 5.2), constrained_layout=True)
    time_min = time_s / 60.0
    indices = nearest_indices(radius_cm, selected_radii_cm)
    for index, color, style in zip(
        indices,
        SERIES_COLORS[: len(indices)],
        SERIES_STYLES[: len(indices)],
        strict=True,
    ):
        actual_radius = radius_cm[index]
        axis.plot(
            time_min,
            moisture_kg_kg[:, index],
            color=color,
            linestyle=style,
            linewidth=2.0,
            label=f"r = {actual_radius:g} cm",
        )
    axis.set(
        title="不同径向位置的含水率随时间变化",
        xlabel="时间（min）",
        ylabel="干基含水率（kg/kg）",
        xlim=(0.0, time_min[-1]),
    )
    axis.grid(True, color="#D9DEE5", linewidth=0.7, alpha=0.8)
    axis.legend(ncol=2, loc="lower left")
    return figure


def plot_spatiotemporal_heatmaps(
    time_s: np.ndarray,
    radius_cm: np.ndarray,
    temperature_c: np.ndarray,
    moisture_kg_kg: np.ndarray,
) -> plt.Figure:
    figure, axes = plt.subplots(2, 1, figsize=(9.0, 7.6), constrained_layout=True)
    time_min = time_s / 60.0
    temperature_mesh = axes[0].pcolormesh(
        time_min,
        radius_cm,
        temperature_c.T,
        shading="auto",
        cmap="magma",
    )
    axes[0].set(
        title="温度时空分布",
        xlabel="时间（min）",
        ylabel="径向位置（cm）",
    )
    temperature_bar = figure.colorbar(temperature_mesh, ax=axes[0], pad=0.02)
    temperature_bar.set_label("温度（°C）")

    moisture_mesh = axes[1].pcolormesh(
        time_min,
        radius_cm,
        moisture_kg_kg.T,
        shading="auto",
        cmap="viridis",
    )
    axes[1].set(
        title="含水率时空分布",
        xlabel="时间（min）",
        ylabel="径向位置（cm）",
    )
    moisture_bar = figure.colorbar(moisture_mesh, ax=axes[1], pad=0.02)
    moisture_bar.set_label("干基含水率（kg/kg）")
    ticks = np.arange(0.0, time_min[-1] + 0.1, 5.0)
    for axis in axes:
        axis.set_xlim(0.0, time_min[-1])
        axis.set_xticks(ticks)
    return figure


def plot_radial_profiles(
    time_s: np.ndarray,
    radius_cm: np.ndarray,
    temperature_c: np.ndarray,
    moisture_kg_kg: np.ndarray,
    profile_times_s: tuple[int, ...],
) -> plt.Figure:
    figure, axes = plt.subplots(1, 2, figsize=(11.0, 4.6), constrained_layout=True)
    indices = nearest_indices(time_s, tuple(float(value) for value in profile_times_s))
    colors = plt.get_cmap("cividis")(np.linspace(0.12, 0.88, len(indices)))
    for index, color in zip(indices, colors, strict=True):
        label = f"{time_s[index] / 60:g} min"
        axes[0].plot(radius_cm, temperature_c[index], color=color, linewidth=2, label=label)
        axes[1].plot(radius_cm, moisture_kg_kg[index], color=color, linewidth=2, label=label)

    axes[0].set(
        title="不同时刻的径向温度剖面",
        xlabel="径向位置（cm）",
        ylabel="温度（°C）",
        xlim=(radius_cm[0], radius_cm[-1]),
    )
    axes[1].set(
        title="不同时刻的径向含水率剖面",
        xlabel="径向位置（cm）",
        ylabel="干基含水率（kg/kg）",
        xlim=(radius_cm[0], radius_cm[-1]),
    )
    for axis in axes:
        axis.grid(True, color="#D9DEE5", linewidth=0.7, alpha=0.8)
    axes[1].legend(title="时刻", loc="lower left")
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
    configure_style()
    time_s, radius_cm, temperature_c, moisture_kg_kg = load_results(args.results)
    ambient_temperature, _, _ = load_environment(args.environment, float(time_s[-1]))
    selected_radii_cm = tuple(float(value) for value in args.radii)
    if len(selected_radii_cm) > len(SERIES_COLORS):
        raise ValueError(f"位置曲线最多绘制 {len(SERIES_COLORS)} 条，以保证图例清晰")
    profile_times_s = tuple(int(value) for value in args.profile_times)
    formats = tuple(dict.fromkeys(args.formats))

    figures = {
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
