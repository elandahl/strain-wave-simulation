"""Save and load portable strain profiles."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from strain_wave.models.result import FORMAT_VERSION, StrainProfile


def save_strain_profile(profile: StrainProfile, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        path,
        format_version=profile.format_version,
        model=profile.model,
        film_material=profile.film_material,
        substrate_material=profile.substrate_material,
        z=profile.z,
        displacement=profile.displacement,
        strain=profile.strain,
        substrate_strain=profile.substrate_strain,
        dz=profile.dz,
        n_bin_film=profile.n_bin_film,
        L_film=profile.L_film,
        L_sub=profile.L_sub,
        t_max=profile.t_max,
        metadata_json=json.dumps(profile.metadata),
    )


def load_strain_profile(path: str | Path) -> StrainProfile:
    path = Path(path)
    with np.load(path, allow_pickle=False) as data:
        metadata = json.loads(str(data["metadata_json"]))
        return StrainProfile(
            format_version=str(data["format_version"]),
            model=str(data["model"]),
            film_material=str(data["film_material"]),
            substrate_material=str(data["substrate_material"]),
            z=data["z"],
            displacement=data["displacement"],
            strain=data["strain"],
            substrate_strain=data["substrate_strain"],
            dz=float(data["dz"]),
            n_bin_film=int(data["n_bin_film"]),
            L_film=float(data["L_film"]),
            L_sub=float(data["L_sub"]),
            t_max=float(data["t_max"]),
            metadata=metadata,
        )


def validate_profile(profile: StrainProfile) -> None:
    if profile.format_version != FORMAT_VERSION:
        raise ValueError(
            f"Unsupported strain profile format {profile.format_version!r}; "
            f"expected {FORMAT_VERSION!r}."
        )
