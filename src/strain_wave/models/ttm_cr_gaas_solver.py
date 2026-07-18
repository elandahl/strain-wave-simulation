"""Two-temperature model and elastic wave solver for Cr/GaAs."""

import numpy as np
from numba import njit

V_GAAS = 4726.5


@njit
def source(R, knu, J, tp, z_film, t):
    beta = 4 * np.log(2)

    I = (beta / np.pi) ** 0.5 * (1 - R) * J * knu / tp * np.exp(
        -knu * z_film - beta * ((t - 2 * tp) / tp) ** 2
    )

    return I


@njit
def solver(
    G, R_ps, k_e_factor, k_s_factor, J, d, d_v, t_Cr, t_sub, dz, t_max,
    monitor_idx=-1,
):
    # monitor_idx: optional grid index (in the substrate) at which the strain
    # and substrate temperature time-histories are recorded every step. Used
    # by the d'Alembert reconstruction model; recording does not affect the
    # solution. Pass -1 to disable (empty arrays are returned).
    R_Cr = 0.6
    knu_Cr = (4 * np.pi / (800e-9) * 3.135)
    tp_laser = 50e-15

    k_Cr = 93.7
    Cp_Cr = 460
    rho_Cr = 7190 * d
    C_Cr = rho_Cr * Cp_Cr

    alpha_Cr = k_Cr / C_Cr
    Mol_M_Cr = 51.9961 * 1e-3
    gamma_Cr = rho_Cr / Mol_M_Cr * 1.4 * 1e-3

    v_Cr = 6477 * d_v
    beta_Cr = 4.9e-6
    B_Cr = rho_Cr * (v_Cr**2)

    factor = k_e_factor
    k_p = k_Cr
    C_p = C_Cr
    alpha_p = k_p / C_p

    v_GaAs = V_GAAS
    rho_GaAs = 5.32 * 1e3
    beta_GaAs = 5.73e-6
    B_GaAs = rho_GaAs * (v_GaAs**2)

    k_s = 55 * k_s_factor
    C_s = 350 * rho_GaAs
    alpha_s = k_s / C_s

    L_film = t_Cr
    sub_film = t_sub
    L_tot = L_film + sub_film

    n_total = int(L_tot / dz)
    z = np.linspace(0, L_tot, n_total)

    n_bin_film = int(L_film / dz)

    T0 = 300

    T_p = np.zeros(n_total) + T0
    T_e = np.zeros(n_total) + T0
    T_s = np.zeros(n_total) + T0

    Tp1_p = np.zeros(n_bin_film)
    Tp1_e = np.zeros(n_bin_film)
    Tp1_s = np.zeros(n_total)

    uf = np.zeros(n_total)
    ufp1 = np.zeros(n_total)
    ufm1 = np.zeros(n_total)

    us = np.zeros(n_total)
    usp1 = np.zeros(n_total)
    usm1 = np.zeros(n_total)

    dt = dz**2 / 2 / 0.002
    n_iter = int(t_max / dt)

    n_rec = n_iter if monitor_idx >= 0 else 0
    strain_rec = np.zeros(n_rec)
    Ts_rec = np.zeros(n_rec)

    for i in range(0, n_iter):
        C_e = gamma_Cr * T_e
        k_e = k_Cr * T_e / T_p / factor
        alpha_e = k_e / C_e

        Tp1_e[0] = (
            T_e[0]
            + 2 * alpha_e[0] * dt / (dz**2) * (T_e[1] - T_e[0])
            - dt * G / C_e[0] * (T_e[0] - T_p[0])
            + source(R_Cr, knu_Cr, J, tp_laser, 0, i * dt) * dt / C_e[0]
        )

        Tp1_p[0] = (
            T_p[0]
            + 2 * alpha_p * dt / (dz**2) * (T_p[1] - T_p[0])
            + dt * G / C_p * (T_e[0] - T_p[0])
        )

        if i == 0:
            ufp1[0] = uf[0] + 0.5 * v_Cr**2 * dt**2 / (dz**2) * (uf[1] - uf[0])
        else:
            ufp1[0] = (
                2 * uf[0]
                - ufm1[0]
                + v_Cr**2 * dt**2 / (dz**2) * (uf[1] - uf[0])
                - 3 * beta_Cr * B_Cr * (dt**2) / dz * (T_p[0] - T0) / rho_Cr
            )

        Tp1_e[1 : n_bin_film - 1] = (
            T_e[1 : n_bin_film - 1]
            + alpha_e[1 : n_bin_film - 1]
            * dt
            / (dz**2)
            * (
                T_e[2:n_bin_film]
                - 2 * T_e[1 : n_bin_film - 1]
                + T_e[0 : n_bin_film - 2]
            )
            - G / C_e[1 : n_bin_film - 1]
            * dt
            * (T_e[1 : n_bin_film - 1] - T_p[1 : n_bin_film - 1])
            + source(
                R_Cr,
                knu_Cr,
                J,
                tp_laser,
                z[1 : n_bin_film - 1],
                i * dt,
            )
            * dt
            / C_e[1 : n_bin_film - 1]
        )

        Tp1_p[1 : n_bin_film - 1] = (
            T_p[1 : n_bin_film - 1]
            + alpha_p
            * dt
            / (dz**2)
            * (
                T_p[2:n_bin_film]
                - 2 * T_p[1 : n_bin_film - 1]
                + T_p[0 : n_bin_film - 2]
            )
            + G / C_p * dt * (T_e[1 : n_bin_film - 1] - T_p[1 : n_bin_film - 1])
        )

        if i == 0:
            ufp1[1 : n_bin_film - 1] = (
                uf[1 : n_bin_film - 1]
                + 0.5
                * v_Cr**2
                * (dt**2)
                / (dz**2)
                * (
                    uf[2:n_bin_film]
                    - 2 * uf[1 : n_bin_film - 1]
                    + uf[0 : n_bin_film - 2]
                )
            )
        else:
            ufp1[1 : n_bin_film - 1] = (
                2 * uf[1 : n_bin_film - 1]
                - ufm1[1 : n_bin_film - 1]
                + v_Cr**2 * (dt**2) / (dz**2)
                * (
                    uf[2:n_bin_film]
                    - 2 * uf[1 : n_bin_film - 1]
                    + uf[0 : n_bin_film - 2]
                )
                - 3
                * beta_Cr
                * B_Cr
                * (T_p[2:n_bin_film] - T_p[0 : n_bin_film - 2])
                * dt**2
                / (2 * dz)
                / rho_Cr
            )

        Tp1_e[n_bin_film - 1] = (
            T_e[n_bin_film - 1]
            + 2
            * alpha_e[n_bin_film - 1]
            * dt
            / (dz**2)
            * (T_e[n_bin_film - 2] - T_e[n_bin_film - 1])
            - G / C_e[n_bin_film - 1]
            * dt
            * (T_e[n_bin_film - 1] - T_p[n_bin_film - 1])
            + source(
                R_Cr,
                knu_Cr,
                J,
                tp_laser,
                z[n_bin_film - 1],
                i * dt,
            )
            * dt
            / C_e[n_bin_film - 1]
        )

        Tp1_p[n_bin_film - 1] = (
            T_p[n_bin_film - 1]
            + 2
            * alpha_p
            * dt
            / (dz**2)
            * (
                T_p[n_bin_film - 2]
                - T_p[n_bin_film - 1]
                - (T_p[n_bin_film - 1] - T_s[n_bin_film - 1]) * dz / R_ps / k_p
            )
            + G / C_p * dt * (T_e[n_bin_film - 1] - T_p[n_bin_film - 1])
        )

        Tp1_s[n_bin_film - 1] = (
            T_s[n_bin_film - 1]
            + 2
            * alpha_s
            * dt
            / (dz**2)
            * (
                T_s[n_bin_film]
                - T_s[n_bin_film - 1]
                + dz / k_s * (1 / R_ps * (T_p[n_bin_film - 1] - T_s[n_bin_film - 1]))
            )
        )

        if i == 0:
            ufIp1 = (
                uf[n_bin_film - 1]
                + rho_GaAs
                * (v_GaAs**2)
                / (rho_Cr * (v_Cr**2))
                * (us[n_bin_film] - us[n_bin_film - 1])
                + 3
                * beta_Cr
                * B_Cr
                / (rho_Cr * (v_Cr**2))
                * (T_p[n_bin_film - 1] - T0)
                * dz
                - 3
                * beta_GaAs
                * B_GaAs
                / (rho_Cr * (v_Cr**2))
                * (T_s[n_bin_film - 1] - T0)
                * dz
            )

            ufp1[n_bin_film - 1] = (
                uf[n_bin_film - 1]
                + 0.5 * v_Cr**2 * (dt / dz) ** 2
                * (ufIp1 - 2 * uf[n_bin_film - 1] + uf[n_bin_film - 2])
                + 1.5
                * dt**2
                * beta_Cr
                * B_Cr
                / rho_Cr
                * (1 / R_ps / k_p)
                * (T_p[n_bin_film - 1] - T_s[n_bin_film - 1])
            )

            usp1[n_bin_film - 1] = ufp1[n_bin_film - 1]

        else:
            ufIp1 = (
                uf[n_bin_film - 1]
                + rho_GaAs
                * (v_GaAs**2)
                / (rho_Cr * (v_Cr**2))
                * (us[n_bin_film] - us[n_bin_film - 1])
                + 3
                * beta_Cr
                * B_Cr
                / (rho_Cr * (v_Cr**2))
                * (T_p[n_bin_film - 1] - T0)
                * dz
                - 3
                * beta_GaAs
                * B_GaAs
                / (rho_Cr * (v_Cr**2))
                * (T_s[n_bin_film - 1] - T0)
                * dz
            )

            ufp1[n_bin_film - 1] = (
                2 * uf[n_bin_film - 1]
                - ufm1[n_bin_film - 1]
                + v_Cr**2 * (dt / dz) ** 2
                * (ufIp1 - 2 * uf[n_bin_film - 1] + uf[n_bin_film - 2])
                + 3
                * dt**2
                * beta_Cr
                * B_Cr
                / rho_Cr
                * (1 / (R_ps * k_p))
                * (T_p[n_bin_film - 1] - T_s[n_bin_film - 1])
            )

            usp1[n_bin_film - 1] = ufp1[n_bin_film - 1]

        Tp1_s[n_bin_film : n_total - 1] = (
            T_s[n_bin_film : n_total - 1]
            + alpha_s
            * dt
            / (dz**2)
            * (
                T_s[n_bin_film + 1 : n_total]
                - 2 * T_s[n_bin_film : n_total - 1]
                + T_s[n_bin_film - 1 : n_total - 2]
            )
        )

        if i == 0:
            usp1[n_bin_film : n_total - 1] = (
                us[n_bin_film : n_total - 1]
                + 0.5
                * v_GaAs**2
                * (dt / dz) ** 2
                * (
                    us[n_bin_film + 1 : n_total]
                    - 2 * us[n_bin_film : n_total - 1]
                    + us[n_bin_film - 1 : n_total - 2]
                )
            )
        else:
            usp1[n_bin_film : n_total - 1] = (
                2 * us[n_bin_film : n_total - 1]
                - usm1[n_bin_film : n_total - 1]
                + v_GaAs**2 * (dt / dz) ** 2
                * (
                    us[n_bin_film + 1 : n_total]
                    - 2 * us[n_bin_film : n_total - 1]
                    + us[n_bin_film - 1 : n_total - 2]
                )
                - 3
                / 2
                * beta_GaAs
                * B_GaAs
                * dt**2
                / dz
                * (T_s[n_bin_film + 1 : n_total] - T_s[n_bin_film - 1 : n_total - 2])
                / rho_GaAs
            )

        Tp1_s[n_total - 1] = (
            T_s[n_total - 1]
            + 2 * alpha_s * dt / (dz**2) * (-T_s[n_total - 1] + T_s[n_total - 2])
        )

        usp1[n_total - 1] = us[n_total - 2] + (v_GaAs * dt / dz - 1) / (
            v_GaAs * dt / dz + 1
        ) * (usp1[n_total - 2] - us[n_total - 1])

        T_p = Tp1_p.copy()
        T_e = Tp1_e.copy()
        T_s = Tp1_s.copy()

        ufm1 = uf.copy()
        uf = ufp1.copy()

        usm1 = us.copy()
        us = usp1.copy()

        if monitor_idx >= 0:
            strain_rec[i] = (us[monitor_idx + 1] - us[monitor_idx - 1]) / (2 * dz)
            Ts_rec[i] = T_s[monitor_idx]

    displacement = np.concatenate((uf[: n_bin_film - 1], us[n_bin_film - 1 :]))
    return displacement, T_e, T_p, T_s, n_iter, dt, strain_rec, Ts_rec


def simulation_diagnostics(dz: float, t_max: float, n_iter: int, dt: float) -> float:
    """Return how far the elastic wave reaches in nm."""
    return round(n_iter * dt * V_GAAS * 1e9, 3)
