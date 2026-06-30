"""Plotting utilities for strain simulation results."""

from pathlib import Path

import matplotlib.pyplot as plt

from strain_wave.models.result import StrainSimulationResult


def plot_strain_results(
    result: StrainSimulationResult,
    save_path: str | Path | None = None,
    show: bool = True,
) -> plt.Figure:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    axes[0].plot(result.z * 1e9, result.displacement * 1e9, "-k")
    axes[0].set_xlabel("distance (nm)", size=15)
    axes[0].set_ylabel("displacement (nm)", size=15)
    axes[0].tick_params(labelsize=15)

    axes[1].plot(result.z * 1e9, result.strain, "-k")
    axes[1].set_xlabel("distance (nm)", size=15)
    axes[1].set_ylabel("strain", size=15)
    axes[1].tick_params(labelsize=15)

    fig.suptitle(f"{result.model}: {result.film_material}/{result.substrate_material}")
    fig.tight_layout()

    if save_path is not None:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=150, bbox_inches="tight")

    if show:
        plt.show()

    return fig
