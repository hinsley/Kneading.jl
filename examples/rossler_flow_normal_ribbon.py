#!/usr/bin/env python3
"""Render a Rössler attractor with an accepted-step flow-normal tangent ribbon."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection


@dataclass(frozen=True)
class RosslerParameters:
    a: float = 0.2
    b: float = 0.2
    c: float = 5.7


@dataclass(frozen=True)
class IntegrationSettings:
    attractor_dt: float = 0.01
    ribbon_dt: float = 0.0025
    attractor_time: float = 1000.0
    ribbon_state_burn_in: float = 200.0
    ribbon_tangent_alignment: float = 80.0
    ribbon_time: float = 30.0
    attractor_plot_stride: int = 2
    ribbon_plot_stride: int = 4
    ribbon_vector_scale: float = 1.0


def rossler_flow(state: np.ndarray, parameters: RosslerParameters) -> np.ndarray:
    x, y, z = state
    return np.array(
        [-y - z, x + parameters.a * y, parameters.b + z * (x - parameters.c)],
        dtype=float,
    )


def rossler_jacobian(
    state: np.ndarray, parameters: RosslerParameters
) -> np.ndarray:
    x, _, z = state
    return np.array(
        [
            [0.0, -1.0, -1.0],
            [1.0, parameters.a, 0.0],
            [z, 0.0, x - parameters.c],
        ],
        dtype=float,
    )


def rk4_state_step(
    state: np.ndarray, dt: float, parameters: RosslerParameters
) -> np.ndarray:
    k1 = rossler_flow(state, parameters)
    k2 = rossler_flow(state + 0.5 * dt * k1, parameters)
    k3 = rossler_flow(state + 0.5 * dt * k2, parameters)
    k4 = rossler_flow(state + dt * k3, parameters)
    return state + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def rk4_augmented_step(
    state: np.ndarray,
    tangent: np.ndarray,
    dt: float,
    parameters: RosslerParameters,
) -> tuple[np.ndarray, np.ndarray]:
    k1_state = rossler_flow(state, parameters)
    k1_tangent = rossler_jacobian(state, parameters) @ tangent

    state_2 = state + 0.5 * dt * k1_state
    tangent_2 = tangent + 0.5 * dt * k1_tangent
    k2_state = rossler_flow(state_2, parameters)
    k2_tangent = rossler_jacobian(state_2, parameters) @ tangent_2

    state_3 = state + 0.5 * dt * k2_state
    tangent_3 = tangent + 0.5 * dt * k2_tangent
    k3_state = rossler_flow(state_3, parameters)
    k3_tangent = rossler_jacobian(state_3, parameters) @ tangent_3

    state_4 = state + dt * k3_state
    tangent_4 = tangent + dt * k3_tangent
    k4_state = rossler_flow(state_4, parameters)
    k4_tangent = rossler_jacobian(state_4, parameters) @ tangent_4

    next_state = state + (dt / 6.0) * (
        k1_state + 2.0 * k2_state + 2.0 * k3_state + k4_state
    )
    next_tangent = tangent + (dt / 6.0) * (
        k1_tangent + 2.0 * k2_tangent + 2.0 * k3_tangent + k4_tangent
    )
    return next_state, next_tangent


def project_flow_normal(
    tangent: np.ndarray,
    state: np.ndarray,
    parameters: RosslerParameters,
) -> tuple[np.ndarray, float]:
    flow = rossler_flow(state, parameters)
    flow_squared = float(flow @ flow)
    if flow_squared <= np.finfo(float).eps:
        raise FloatingPointError("The flow is too small for tangent projection.")

    projected = tangent - (float(tangent @ flow) / flow_squared) * flow
    projected_norm = float(np.linalg.norm(projected))
    if projected_norm <= np.finfo(float).eps:
        raise FloatingPointError("The projected tangent is too small to normalize.")
    return projected / projected_norm, projected_norm


def integrate_state(
    initial_state: np.ndarray,
    duration: float,
    dt: float,
    parameters: RosslerParameters,
    sample_stride: int,
) -> np.ndarray:
    steps = int(round(duration / dt))
    sample_count = steps // sample_stride + 1
    samples = np.empty((sample_count, 3), dtype=float)
    samples[0] = initial_state

    state = initial_state.copy()
    sample_index = 1
    for step in range(1, steps + 1):
        state = rk4_state_step(state, dt, parameters)
        if step % sample_stride == 0:
            samples[sample_index] = state
            sample_index += 1
    return samples[:sample_index]


def advance_state(
    initial_state: np.ndarray,
    duration: float,
    dt: float,
    parameters: RosslerParameters,
) -> np.ndarray:
    state = initial_state.copy()
    for _ in range(int(round(duration / dt))):
        state = rk4_state_step(state, dt, parameters)
    return state


def integrate_flow_normal_ribbon(
    initial_state: np.ndarray,
    parameters: RosslerParameters,
    settings: IntegrationSettings,
) -> tuple[np.ndarray, np.ndarray, dict[str, float]]:
    state = advance_state(
        initial_state,
        settings.ribbon_state_burn_in,
        settings.ribbon_dt,
        parameters,
    )
    tangent, _ = project_flow_normal(
        np.array([1.0, 0.0, 0.0], dtype=float), state, parameters
    )

    alignment_steps = int(
        round(settings.ribbon_tangent_alignment / settings.ribbon_dt)
    )
    for _ in range(alignment_steps):
        state, tangent = rk4_augmented_step(
            state, tangent, settings.ribbon_dt, parameters
        )
        tangent, _ = project_flow_normal(tangent, state, parameters)

    record_steps = int(round(settings.ribbon_time / settings.ribbon_dt))
    sample_count = record_steps // settings.ribbon_plot_stride + 1
    states = np.empty((sample_count, 3), dtype=float)
    tangents = np.empty((sample_count, 3), dtype=float)
    states[0] = state
    tangents[0] = tangent

    log_growth = 0.0
    max_unit_error = 0.0
    max_flow_dot = 0.0
    sample_index = 1
    for step in range(1, record_steps + 1):
        state, tangent = rk4_augmented_step(
            state, tangent, settings.ribbon_dt, parameters
        )
        tangent, growth = project_flow_normal(tangent, state, parameters)
        log_growth += np.log(growth)
        max_unit_error = max(max_unit_error, abs(np.linalg.norm(tangent) - 1.0))
        max_flow_dot = max(
            max_flow_dot, abs(float(tangent @ rossler_flow(state, parameters)))
        )

        if step % settings.ribbon_plot_stride == 0:
            states[sample_index] = state
            tangents[sample_index] = tangent
            sample_index += 1

    diagnostics = {
        "max_unit_error": max_unit_error,
        "max_flow_dot": max_flow_dot,
        "transverse_growth_rate": log_growth / settings.ribbon_time,
    }
    return states[:sample_index], tangents[:sample_index], diagnostics


def y_minimum_half_nullcline(
    parameters: RosslerParameters,
) -> np.ndarray:
    yz_vertices = np.array(
        [
            [-12.0, -3.0],
            [0.0, -3.0],
            [0.0, 5.0],
            [-12.0, 5.0],
        ],
        dtype=float,
    )
    return np.column_stack(
        (-parameters.a * yz_vertices[:, 0], yz_vertices[:, 0], yz_vertices[:, 1])
    )


def add_half_nullcline(
    axis,
    parameters: RosslerParameters,
) -> None:
    vertices = y_minimum_half_nullcline(parameters)
    face = Poly3DCollection(
        [vertices],
        facecolor="#171923",
        edgecolor="#20232d",
        linewidth=1.0,
        alpha=0.58,
        zorder=1,
    )
    axis.add_collection3d(face)

    mesh_segments = []
    for y_value in np.linspace(-12.0, 0.0, 7):
        mesh_segments.append(
            [
                (-parameters.a * y_value, y_value, -3.0),
                (-parameters.a * y_value, y_value, 5.0),
            ]
        )
    for z_value in np.linspace(-3.0, 5.0, 5):
        mesh_segments.append(
            [
                (-parameters.a * -12.0, -12.0, z_value),
                (0.0, 0.0, z_value),
            ]
        )
    axis.add_collection3d(
        Line3DCollection(mesh_segments, colors="#3d414c", linewidths=0.45, alpha=0.52)
    )


def render_figure(
    attractor: np.ndarray,
    ribbon_states: np.ndarray,
    ribbon_tangents: np.ndarray,
    parameters: RosslerParameters,
    settings: IntegrationSettings,
    output_path: Path,
    pdf_path: Path | None,
) -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 11,
            "axes.labelsize": 14,
            "axes.titlesize": 18,
        }
    )
    figure = plt.figure(figsize=(10.8, 9.6), constrained_layout=True)
    axis = figure.add_subplot(111, projection="3d")
    figure.patch.set_facecolor("white")
    axis.set_facecolor("white")

    attractor_line = axis.plot(
        attractor[:, 0],
        attractor[:, 1],
        attractor[:, 2],
        color="#4169e1",
        linewidth=0.38,
        alpha=0.30,
        zorder=2,
    )[0]
    attractor_line.set_rasterized(True)

    add_half_nullcline(axis, parameters)

    tangent_edge = ribbon_states + settings.ribbon_vector_scale * ribbon_tangents
    ribbon_x = np.column_stack((ribbon_states[:, 0], tangent_edge[:, 0]))
    ribbon_y = np.column_stack((ribbon_states[:, 1], tangent_edge[:, 1]))
    ribbon_z = np.column_stack((ribbon_states[:, 2], tangent_edge[:, 2]))
    axis.plot_surface(
        ribbon_x,
        ribbon_y,
        ribbon_z,
        color="#ff3b30",
        alpha=0.38,
        shade=False,
        linewidth=0,
        antialiased=True,
        zorder=4,
    )
    axis.plot(
        ribbon_states[:, 0],
        ribbon_states[:, 1],
        ribbon_states[:, 2],
        color="#6e001f",
        linewidth=1.0,
        alpha=0.95,
        zorder=6,
    )
    axis.plot(
        tangent_edge[:, 0],
        tangent_edge[:, 1],
        tangent_edge[:, 2],
        color="#c1123f",
        linewidth=0.72,
        alpha=0.94,
        zorder=5,
    )

    crossbar_indices = np.arange(0, len(ribbon_states), 40)
    crossbars = [
        [ribbon_states[index], tangent_edge[index]] for index in crossbar_indices
    ]
    axis.add_collection3d(
        Line3DCollection(
            crossbars,
            colors="#ff5a36",
            linewidths=0.72,
            alpha=0.94,
            zorder=7,
        )
    )

    axis.set_xlim(-10.5, 12.0)
    axis.set_ylim(-12.5, 9.5)
    axis.set_zlim(-3.2, 24.0)
    axis.set_box_aspect((22.5, 22.0, 27.2), zoom=0.88)
    axis.view_init(elev=22.5, azim=-102.0, roll=0.0)
    axis.set_xlabel("x", labelpad=10)
    axis.set_ylabel("y", labelpad=10)
    axis.set_zlabel("z", labelpad=8)
    axis.set_title(
        "Unstable flow-normal tangent ribbon on the Rössler attractor",
        pad=18,
        weight="semibold",
    )

    for pane in (axis.xaxis.pane, axis.yaxis.pane, axis.zaxis.pane):
        pane.set_facecolor((1.0, 1.0, 1.0, 0.0))
        pane.set_edgecolor((0.80, 0.81, 0.84, 0.55))
    for coordinate_axis in (axis.xaxis, axis.yaxis, axis.zaxis):
        coordinate_axis._axinfo["grid"]["color"] = (0.68, 0.70, 0.74, 0.32)
        coordinate_axis._axinfo["grid"]["linewidth"] = 0.55
        coordinate_axis._axinfo["axisline"]["color"] = (0.25, 0.27, 0.31, 0.68)
    axis.tick_params(colors="#343741", labelsize=10, pad=1)

    legend_handles = [
        Line2D([0], [0], color="#4169e1", linewidth=1.2, alpha=0.70),
        Patch(facecolor="#ff3b30", edgecolor="#c1123f", alpha=0.48),
        Patch(facecolor="#171923", edgecolor="#20232d", alpha=0.64),
    ]
    legend_labels = [
        r"long trajectory ($t=1000$)",
        r"transported unit flow-normal ribbon ($t=30$)",
        r"$y$-minimum half-nullcline: $\dot y=0$, $y\leq 0$",
    ]
    axis.legend(
        legend_handles,
        legend_labels,
        loc="upper left",
        bbox_to_anchor=(0.015, 0.965),
        frameon=True,
        framealpha=0.94,
        facecolor="white",
        edgecolor="#d5d7dd",
        borderpad=0.75,
        labelspacing=0.65,
    )
    axis.text2D(
        0.025,
        0.025,
        "State burn-in: 200; tangent relaxation: 80; "
        r"$\Delta t_{\mathrm{ribbon}}=0.0025$" "\n"
        r"Accepted-step gauge: $\|v\|=1$ and $v\!\cdot\!f=0$",
        transform=axis.transAxes,
        color="#4a1c2b",
        bbox={
            "boxstyle": "round,pad=0.45",
            "facecolor": "white",
            "edgecolor": "#e0cbd2",
            "alpha": 0.94,
        },
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=300, facecolor=figure.get_facecolor())
    if pdf_path is not None:
        figure.savefig(pdf_path, dpi=300, facecolor=figure.get_facecolor())
    plt.close(figure)


def parse_arguments() -> argparse.Namespace:
    script_directory = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=script_directory / "rossler_flow_normal_ribbon.png",
        help="PNG output path",
    )
    parser.add_argument(
        "--pdf",
        type=Path,
        default=script_directory / "rossler_flow_normal_ribbon.pdf",
        help="PDF output path; pass an empty string to disable it",
    )
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    parameters = RosslerParameters()
    settings = IntegrationSettings()

    attractor = integrate_state(
        np.array([-1.0, 0.0, 0.0], dtype=float),
        settings.attractor_time,
        settings.attractor_dt,
        parameters,
        settings.attractor_plot_stride,
    )
    ribbon_states, ribbon_tangents, diagnostics = integrate_flow_normal_ribbon(
        np.array([0.0, -5.0, 0.0], dtype=float),
        parameters,
        settings,
    )
    pdf_path = None if str(arguments.pdf) == "" else arguments.pdf
    render_figure(
        attractor,
        ribbon_states,
        ribbon_tangents,
        parameters,
        settings,
        arguments.output,
        pdf_path,
    )

    print(f"saved PNG: {arguments.output}")
    if pdf_path is not None:
        print(f"saved PDF: {pdf_path}")
    print(f"attractor samples: {len(attractor)}")
    print(f"ribbon samples: {len(ribbon_states)}")
    print(f"maximum unit-norm error: {diagnostics['max_unit_error']:.3e}")
    print(f"maximum |v dot f|: {diagnostics['max_flow_dot']:.3e}")
    print(
        "finite-time transverse growth rate: "
        f"{diagnostics['transverse_growth_rate']:.6f}"
    )


if __name__ == "__main__":
    main()
