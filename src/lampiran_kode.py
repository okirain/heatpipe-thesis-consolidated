"""Kutipan inti hp_model.py; berkas lengkap tersedia pada folder src."""

# --- Sifat natrium: panas laten, tekanan uap jenuh, temperatur saturasi ---
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

# ...

# --- Geometri pipa kalor dan properti turunannya ---
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

# ...

# --- Penyelesaian aksial uap dan cairan (integrasi mundur Runge-Kutta) ---
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

# ...

# --- Batas kapiler dengan metode bagi dua ---
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

# ...

# --- Konduksi radial pelet bahan bakar anular ---
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

# ...
