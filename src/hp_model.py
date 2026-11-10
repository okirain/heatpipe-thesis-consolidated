"""
hp_model.py -- Reduced-order two-phase model for alkali-metal heat pipes
==========================================================================
Konsolidasi fisika dari enam berkas sesi (yuisa, jeremiah, gamma, beta,
yukita, reduced-order) menjadi satu implementasi tunggal yang dapat
diverifikasi.

Referensi formulasi:
  * Fink & Leibowitz (ANL/RE-95/2, 1995)  -- sifat natrium
  * Cotter (1965); Busse (1973); Levy (1968) -- batas operasi
  * Faghri (2015) -- jaringan hambatan termal, aliran uap 1D
  * Chi (1976) -- korelasi sumbu kasa / Blake-Kozeny

Model ini 1D aksial, tunak, dengan asumsi:
  1. Uap mengikuti kesetimbangan saturasi lokal T = Tsat(P).
  2. Aliran uap laminer/turbulen bergantung Re, dengan koreksi injeksi
     (evaporator) dan suction (kondensor) melalui koefisien beta.
  3. Aliran cairan balik mengikuti Darcy-Ergun pada sumbu berpori seragam.
  4. Gravitasi diperhitungkan melalui head hidrostatik.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable

import numpy as np

# ----------------------------------------------------------------------------
# 0. KONSTANTA
# ----------------------------------------------------------------------------
T_BOIL = 1156.2        # K, titik didih natrium pada 101325 Pa
T_CRIT = 2503.7        # K, temperatur kritis natrium
P_ATM = 101325.0       # Pa
R_U = 8.314462618      # J/(mol K)
M_NA = 0.02298977      # kg/mol
R_G = R_U / M_NA       # 361.66 J/(kg K)
GAMMA = 5.0 / 3.0      # rasio panas spesifik uap monatomik

# Toleransi numerik
_TOL_P = 1e-3


# ----------------------------------------------------------------------------
# 1. SIFAT TERMOFISIKA NATRIUM
# ----------------------------------------------------------------------------
def h_fg(T: float) -> float:
    """Panas laten penguapan natrium [J/kg], Fink & Leibowitz (1995)."""
    x = max(1.0e-6, 1.0 - T / T_CRIT)
    return 4.76e6 * x ** 0.38


def _build_psat_table() -> tuple[np.ndarray, np.ndarray]:
    """Integrasi Clausius-Clapeyron dari titik didih normal.

    d(ln P)/dT = h_fg(T) / (R_g T^2)

    Integrasi numeris (aturan trapesium) pada grid halus, dua arah dari
    T_BOIL. Hasilnya konsisten secara termodinamika dengan h_fg(T) yang
    dipakai model, sehingga tidak ada inkonsistensi P_sat vs h_fg.
    """
    Ts = np.arange(400.0, 1500.0 + 0.25, 0.25)
    lnP = np.empty_like(Ts)
    i0 = int(np.argmin(np.abs(Ts - T_BOIL)))
    lnP[i0] = math.log(P_ATM)
    # ke atas
    for i in range(i0 + 1, len(Ts)):
        T1, T2 = Ts[i - 1], Ts[i]
        f = lambda T: h_fg(T) / (R_G * T * T)
        lnP[i] = lnP[i - 1] + 0.5 * (f(T1) + f(T2)) * (T2 - T1)
    # ke bawah
    for i in range(i0 - 1, -1, -1):
        T1, T2 = Ts[i], Ts[i + 1]
        f = lambda T: h_fg(T) / (R_G * T * T)
        lnP[i] = lnP[i + 1] - 0.5 * (f(T1) + f(T2)) * (T2 - T1)
    return Ts, lnP


_TS, _LNP = _build_psat_table()


def psat(T: float) -> float:
    """Tekanan uap jenuh natrium [Pa]."""
    return float(np.exp(np.interp(T, _TS, _LNP)))


def tsat(P: float) -> float:
    """Temperatur saturasi natrium [K] dari tekanan [Pa]."""
    return float(np.interp(math.log(max(P, 1e-6)), _LNP, _TS))


def rho_l(T: float) -> float:
    """Densitas natrium cair [kg/m3] (ANL/RE-95/2, +-1%)."""
    return 1014.0 - 0.235 * T


def rho_v(T: float) -> float:
    """Densitas uap natrium [kg/m3], gas ideal monomer."""
    return psat(T) / (R_G * T)


def mu_l(T: float) -> float:
    """Viskositas natrium cair [Pa s]."""
    return math.exp(-6.4406 - 0.3958 * math.log(T) + 556.835 / T)


def mu_v(T: float) -> float:
    """Viskositas uap natrium [Pa s], Chapman-Enskog (sigma=3.72 A, Om=1)."""
    return 2.6693e-6 * math.sqrt(M_NA * 1e3 * T) / (3.72 ** 2)


def k_l(T: float) -> float:
    """Konduktivitas termal natrium cair [W/(m K)]."""
    return 124.5 - 0.0475 * T


def cp_l(T: float) -> float:
    """Kalor jenis natrium cair [J/(kg K)]."""
    return 1436.0 - 0.5806 * T + 4.627e-4 * T * T


def sigma(T: float) -> float:
    """Tegangan permukaan natrium [N/m]."""
    return 0.2405 * max(1e-6, 1.0 - T / T_CRIT) ** 1.126


def sound_speed(T: float) -> float:
    """Kecepatan suara uap natrium [m/s]."""
    return math.sqrt(GAMMA * R_G * T)


def properties(T: float) -> dict:
    return dict(
        T=T, P_sat=psat(T), rho_l=rho_l(T), rho_v=rho_v(T),
        mu_l=mu_l(T), mu_v=mu_v(T), k_l=k_l(T), cp_l=cp_l(T),
        sigma=sigma(T), h_fg=h_fg(T), c_s=sound_speed(T),
    )


# ----------------------------------------------------------------------------
# 2. GEOMETRI
# ----------------------------------------------------------------------------
@dataclass
class HeatPipe:
    """Geometri pipa kalor (semua panjang dalam meter)."""
    name: str
    L_e: float          # panjang evaporator
    L_a: float          # panjang adiabatik
    L_c: float          # panjang kondensor
    r_o: float          # jari-jari luar selubung
    t_wall: float       # tebal dinding selubung
    t_wick: float       # tebal sumbu (annular)
    r_eff: float        # jari-jari pori efektif meniskus
    K: float            # permeabilitas sumbu [m2]
    eps: float          # porositas sumbu (-)
    k_wall: float       # konduktivitas selubung [W/(m K)]
    rho_wall: float     # densitas selubung [kg/m3]
    k_solid: float = 15.0   # konduktivitas material kasa [W/(m K)]
    r_n: float = 2.5e-6     # jari-jari rongga nukleasi [m]
    fill_frac: float = 0.35  # fraksi area annulus yang terbuka untuk aliran

    # --- properti turunan -------------------------------------------------
    @property
    def L(self) -> float:
        return self.L_e + self.L_a + self.L_c

    @property
    def r_i(self) -> float:
        return self.r_o - self.t_wall

    @property
    def r_v(self) -> float:
        return self.r_i - self.t_wick

    @property
    def D_v(self) -> float:
        return 2.0 * self.r_v

    @property
    def A_v(self) -> float:
        return math.pi * self.r_v ** 2

    @property
    def A_w(self) -> float:
        """Luas penampang aliran cairan efektif pada annulus sumbu."""
        return math.pi * (self.r_i ** 2 - self.r_v ** 2) * self.fill_frac

    @property
    def L_eff(self) -> float:
        return 0.5 * self.L_e + self.L_a + 0.5 * self.L_c

    @property
    def A_wall(self) -> float:
        return math.pi * (self.r_o ** 2 - self.r_i ** 2)

    def mdot_profile(self, z: np.ndarray) -> np.ndarray:
        """Laju alir massa uap ternormalisasi m(z)/m_max."""
        z = np.asarray(z, dtype=float)
        m = np.empty_like(z)
        ev = z <= self.L_e
        ad = (z > self.L_e) & (z <= self.L_e + self.L_a)
        cd = z > self.L_e + self.L_a
        m[ev] = z[ev] / self.L_e
        m[ad] = 1.0
        m[cd] = 1.0 - (z[cd] - self.L_e - self.L_a) / self.L_c
        return np.clip(m, 0.0, 1.0)

    def beta(self, z: np.ndarray) -> np.ndarray:
        """Koefisien momentum: 1.0 (injeksi) -> 0.65 (suction)."""
        z = np.asarray(z, dtype=float)
        b = np.full_like(z, 1.0)
        b[z > self.L_e + self.L_a] = 0.65
        return b

    # --- massa ------------------------------------------------------------
    def wall_mass(self) -> float:
        return self.A_wall * self.L * self.rho_wall

    def wick_mass(self) -> float:
        """Massa struktur kasa (padatan) pada annulus sumbu."""
        V = math.pi * (self.r_i ** 2 - self.r_v ** 2) * self.L
        return V * (1.0 - self.eps) * self.rho_wall * 0.6

    def fluid_mass(self) -> float:
        """Massa natrium terisi (sumbu + sedikit uap)."""
        V_wick = math.pi * (self.r_i ** 2 - self.r_v ** 2) * self.L * self.eps
        V_vap = self.A_v * self.L
        T = 1000.0
        return V_wick * rho_l(T) + V_vap * rho_v(T)

    def total_mass(self) -> float:
        return self.wall_mass() + self.wick_mass() + self.fluid_mass()


# ----------------------------------------------------------------------------
# 3. BATAS OPERASI
# ----------------------------------------------------------------------------
def q_sonic(hp: HeatPipe, T: float) -> float:
    """Batas sonik Levy-Busse [W]."""
    rho = rho_v(T)
    return hp.A_v * rho * h_fg(T) * math.sqrt(
        GAMMA * R_G * T / (2.0 * (GAMMA + 1.0))
    )


def q_entrainment(hp: HeatPipe, T: float) -> float:
    """Batas entrainment Cotter (kriteria Weber = 1) [W]."""
    rho = rho_v(T)
    return hp.A_v * h_fg(T) * math.sqrt(sigma(T) * rho / (2.0 * hp.r_eff))


def q_viscous(hp: HeatPipe, T: float) -> float:
    """Batas viskos Busse [W]."""
    return (hp.r_v ** 2 * h_fg(T) * rho_v(T) * psat(T) * hp.A_v
            / (16.0 * mu_v(T) * hp.L_eff))


def q_boiling(hp: HeatPipe, T: float) -> float:
    """Batas pendidihan (insipiensi nukleat) [W]."""
    dT_crit = 2.0 * sigma(T) * T / (rho_v(T) * h_fg(T) * hp.r_n)
    k_eff = k_wick_eff(hp, T)
    return 2.0 * math.pi * hp.L_e * k_eff * dT_crit / math.log(hp.r_i / hp.r_v)


def k_wick_eff(hp: HeatPipe, T: float) -> float:
    """Konduktivitas efektif sumbu jenuh cairan (Maxwell-Eucken)."""
    kl, ks, e = k_l(T), hp.k_solid, hp.eps
    num = 2 * kl + ks - 2 * (1 - e) * (kl - ks)
    den = 2 * kl + ks + (1 - e) * (kl - ks)
    return kl * num / den


# ----------------------------------------------------------------------------
# 4. PROFIL AKSIAL UAP DAN CAIRAN
# ----------------------------------------------------------------------------
@dataclass
class AxialSolution:
    z: np.ndarray
    P_v: np.ndarray
    T_v: np.ndarray
    u_v: np.ndarray
    Ma: np.ndarray
    Re_v: np.ndarray
    P_l: np.ndarray
    r_c: np.ndarray
    Q: float
    T_evap: float
    T_cond: float
    dP_v: float
    dP_l: float
    dP_g: float

    @property
    def dP_req(self) -> float:
        return self.dP_v + self.dP_l + self.dP_g


def _friction_factor(Re: float) -> float:
    if Re < 1e-9:
        return 0.0
    if Re < 2300.0:
        return 64.0 / Re
    return 0.3164 * Re ** -0.25


def solve_axial(hp: HeatPipe, Q: float, T_cond: float, g: float = 0.0,
                angle_deg: float = 90.0, N: int = 240) -> AxialSolution:
    """Integrasi persamaan momentum 1D dari kondensor (z = L) ke evaporator.

    Keadaan batas: pada z = L, P = P_sat(T_cond) dan u = 0. Integrasi
    mundur (RK4) menghasilkan P_v(z). Temperatur mengikuti saturasi lokal.
    """
    z = np.linspace(0.0, hp.L, N + 1)
    dz = z[1] - z[0]
    mdot_max = Q / h_fg(T_cond)
    m_norm = hp.mdot_profile(z)
    beta = hp.beta(z)
    A_v = hp.A_v

    def rho_of(P: float) -> float:
        T = tsat(P)
        return P / (R_G * T)

    def rhs(P: float, zz: float, m: float, bb: float) -> float:
        if m <= 1e-12:
            return 0.0
        T = tsat(P)
        rho = P / (R_G * T)
        u = m * mdot_max / (rho * A_v)
        Re = rho * u * hp.D_v / mu_v(T)
        f = _friction_factor(Re)
        # d(rho u^2)/dz secara numerik dari profil m(z)
        eps = 1e-4 * hp.L
        z1, z2 = max(0.0, zz - eps), min(hp.L, zz + eps)
        m1 = float(np.interp(z1, z, m_norm))
        m2 = float(np.interp(z2, z, m_norm))
        rho_u2_1 = (m1 * mdot_max) ** 2 / (rho * A_v)
        rho_u2_2 = (m2 * mdot_max) ** 2 / (rho * A_v)
        dinert = (rho_u2_2 - rho_u2_1) / (z2 - z1) if z2 > z1 else 0.0
        return -(f * rho * u * u / (2.0 * hp.D_v)) - bb * dinert

    # --- integrasi mundur: P[i] dari P[i+1]
    P = np.empty(N + 1)
    P[N] = psat(T_cond)
    for i in range(N - 1, -1, -1):
        Pn = P[i + 1]
        k1 = rhs(Pn, z[i + 1], m_norm[i + 1], beta[i + 1])
        k2 = rhs(Pn + 0.5 * dz * k1, z[i] + 0.5 * dz,
                 float(np.interp(z[i] + 0.5 * dz, z, m_norm)),
                 float(np.interp(z[i] + 0.5 * dz, z, beta)))
        k3 = rhs(Pn + 0.5 * dz * k2, z[i] + 0.5 * dz,
                 float(np.interp(z[i] + 0.5 * dz, z, m_norm)),
                 float(np.interp(z[i] + 0.5 * dz, z, beta)))
        k4 = rhs(Pn + dz * k3, z[i], m_norm[i], beta[i])
        # integrasi mundur: P(z) = P(z+dz) - dz * (dP/dz)
        P[i] = Pn - (dz / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        P[i] = max(P[i], 1.0)

    T_v = np.array([tsat(p) for p in P])
    rho = P / (R_G * T_v)
    u = np.where(rho > 0, m_norm * mdot_max / (rho * A_v), 0.0)
    c_s = np.sqrt(GAMMA * R_G * T_v)
    Ma = u / c_s
    Re = rho * u * hp.D_v / np.array([mu_v(t) for t in T_v])

    # --- cairan: Darcy-Ergun, integrasi mundur dari z = L
    T_ref = 0.5 * (T_v[0] + T_v[-1])
    m_l = m_norm  # kekekalan massa: m_l(z) = m_v(z)
    mdot_l = m_l * mdot_max
    mu_l_T = np.array([mu_l(t) for t in T_v])
    rho_l_T = np.array([rho_l(t) for t in T_v])
    dPdz_l = (mu_l_T * mdot_l / (rho_l_T * hp.A_w * hp.K)
              * (1.0 + 0.0 * m_l))
    dP_l_tot = float(np.trapezoid(dPdz_l, z))

    dP_g = rho_l(T_ref) * g * hp.L * math.sin(math.radians(angle_deg))

    # --- jari-jari meniskus lokal r_c(z) = 2 sigma / (P_v - P_l)
    P_l = np.empty_like(P)
    P_l[-1] = P[-1]
    for i in range(N - 1, -1, -1):
        # cairan mengalir dari kondensor ke evaporator: P_l turun searah aliran
        P_l[i] = P_l[i + 1] - 0.5 * dz * (dPdz_l[i] + dPdz_l[i + 1])
    dPc = np.maximum(P - P_l, 1.0)
    r_c = 2.0 * np.array([sigma(t) for t in T_v]) / dPc

    return AxialSolution(
        z=z, P_v=P, T_v=T_v, u_v=u, Ma=Ma, Re_v=Re, P_l=P_l, r_c=r_c,
        Q=Q, T_evap=float(T_v[0]), T_cond=float(T_v[-1]),
        dP_v=float(np.max(P) - np.min(P)), dP_l=dP_l_tot, dP_g=dP_g,
    )


def dP_cap_max(hp: HeatPipe, T: float) -> float:
    return 2.0 * sigma(T) / hp.r_eff


def q_capillary(hp: HeatPipe, T: float, g: float = 0.0,
                angle_deg: float = 90.0) -> float:
    """Batas kapiler [W] dengan bisection pada Q.

    Delta P_req(Q) monoton naik; cari Q sehingga dP_req = dP_c,max.
    """
    cap = dP_cap_max(hp, T) - rho_l(T) * g * hp.L * math.sin(math.radians(angle_deg))
    if cap <= 0.0:
        return 0.0

    def req(Q: float) -> float:
        s = solve_axial(hp, Q, T, g=g, angle_deg=angle_deg, N=80)
        return s.dP_req

    if req(1.0) > cap:
        return 0.0
    lo, hi = 1.0, 1.0
    while req(hi) < cap and hi < 5.0e5:
        hi *= 1.6
    if hi >= 5.0e5:
        return 5.0e5
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if req(mid) < cap:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-3 * hi:
            break
    return 0.5 * (lo + hi)


@dataclass
class Envelope:
    T: float
    Q_cap: float
    Q_sonic: float
    Q_ent: float
    Q_boil: float
    Q_visc: float

    @property
    def Q_min(self) -> float:
        return min(self.Q_cap, self.Q_sonic, self.Q_ent, self.Q_boil, self.Q_visc)

    @property
    def limiting(self) -> str:
        pairs = dict(cap=self.Q_cap, sonic=self.Q_sonic, entrainment=self.Q_ent,
                     boiling=self.Q_boil, viscous=self.Q_visc)
        return min(pairs, key=pairs.get)


def envelope(hp: HeatPipe, T: float, g: float = 0.0,
             angle_deg: float = 90.0) -> Envelope:
    return Envelope(
        T=T,
        Q_cap=q_capillary(hp, T, g=g, angle_deg=angle_deg),
        Q_sonic=q_sonic(hp, T),
        Q_ent=q_entrainment(hp, T),
        Q_boil=q_boiling(hp, T),
        Q_visc=q_viscous(hp, T),
    )


# ----------------------------------------------------------------------------
# 5. JARINGAN HAMBATAN TERMAL DAN KONDUKSI BAHAN BAKAR
# ----------------------------------------------------------------------------
def thermal_resistance(hp: HeatPipe, T: float) -> dict:
    """Hambatan termal radial + aksial pipa kalor [K/W]."""
    kw = hp.k_wall
    kwe = k_wick_eff(hp, T)
    R_wall = math.log(hp.r_o / hp.r_i) / (2 * math.pi * kw * hp.L_e)
    R_wick = math.log(hp.r_i / hp.r_v) / (2 * math.pi * kwe * hp.L_e)
    R_cond_wall = math.log(hp.r_o / hp.r_i) / (2 * math.pi * kw * hp.L_c)
    R_vap = 8.0 * mu_v(T) * T / (math.pi * rho_v(T) ** 2 * h_fg(T) ** 2 * hp.r_v ** 4)
    return dict(R_wall=R_wall, R_wick=R_wick, R_cond_wall=R_cond_wall,
                R_vap=R_vap, R_total=R_wall + R_wick + R_cond_wall + R_vap)


def annular_fuel_dT(Q: float, r_i: float, r_o: float, H: float,
                    k_fuel: float = 16.0) -> dict:
    """Kenaikan temperatur puncak pelet UN anular (adiabatik di OD).

    Solusi tunak 1D radial dengan pembangkitan volumetrik seragam q''':
        dT/dr = q'''/(2k) (r_o^2/r - r)
        dT_puncak = q'''/(2k) [r_o^2 ln(r_o/r_i) - (r_o^2 - r_i^2)/2]
    """
    V = math.pi * (r_o ** 2 - r_i ** 2) * H
    q3 = Q / V
    dT = q3 / (2 * k_fuel) * (r_o ** 2 * math.log(r_o / r_i)
                              - 0.5 * (r_o ** 2 - r_i ** 2))
    return dict(dT=dT, q3=q3, V=V)


# ----------------------------------------------------------------------------
# 6. DESAIN ACUAN
# ----------------------------------------------------------------------------
def krusty_pipe() -> HeatPipe:
    """Pipa kalor uji NASA KRUSTY (Haynes 230, kasa 400 mesh)."""
    return HeatPipe(
        name="KRUSTY-HP",
        L_e=0.25, L_a=0.30, L_c=0.45,
        r_o=0.00794, t_wall=0.00089, t_wick=0.00060,
        r_eff=3.0e-5, K=1.2e-10, eps=0.70,
        k_wall=24.0, rho_wall=8970.0, k_solid=15.0,
        fill_frac=1.0,
    )


def isr_pipe(OD: float = 0.018, t_wall: float = 0.0009,
             t_wick: float = 0.0008, K: float = 5.0e-9,
             r_eff: float = 2.0e-5) -> HeatPipe:
    """Pipa kalor natrium in-pin untuk teras ISR (selubung Mo-30Re)."""
    return HeatPipe(
        name=f"ISR-HP-OD{OD*1e3:.1f}",
        L_e=0.25, L_a=0.05, L_c=0.16,
        r_o=OD / 2.0, t_wall=t_wall, t_wick=t_wick,
        r_eff=r_eff, K=K, eps=0.55,
        k_wall=95.0, rho_wall=10800.0, k_solid=30.0,
        fill_frac=0.40,
    )
