"""CLI entry point."""

import argparse
from pathlib import Path

from strain_wave.config import SimulationConfig
from strain_wave.io import save_strain_profile
from strain_wave.models.base import list_models
from strain_wave.pipeline import run_simulation
from strain_wave.plotting import plot_strain_results


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run strain and elastic wave simulation."
    )
    parser.add_argument(
        "--model",
        default="ttm_cr_gaas",
        choices=list_models(),
        help="Strain simulation model",
    )
    parser.add_argument(
        "--t-max",
        type=float,
        default=300e-12,
        help="Simulation end time in seconds",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/strain_profile.npz"),
        help="Path to save strain profile (.npz)",
    )
    parser.add_argument(
        "--figure",
        type=Path,
        default=Path("results/strain_figure.png"),
        help="Path to save displacement/strain figure",
    )
    parser.add_argument("--no-show", action="store_true")
    args = parser.parse_args()

    config = SimulationConfig(model=args.model, t_max=args.t_max)
    result = run_simulation(config=config)
    save_strain_profile(result.to_profile(), args.output)
    plot_strain_results(result, save_path=args.figure, show=not args.no_show)

    print(f"Saved strain profile to {args.output}")
    print(f"Saved figure to {args.figure}")


if __name__ == "__main__":
    main()
