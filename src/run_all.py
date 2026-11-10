"""
run_all.py -- Menjalankan seluruh studi dan menghasilkan tabel + gambar
untuk skripsi konsolidasi.

Keluaran:
  results/*.csv   tabel siap pakai
  figures/*.pdf   gambar
  results/ringkasan.json  angka kunci
"""
from __future__ import annotations

import json
import math
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import hp_model as m

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RES = os.path.join(ROOT, "results")
FIG = os.path.join(ROOT, "figures")
os.makedirs(RES, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

# ---------------------------------------------------------------------------
# Gaya gambar
# ---------------------------------------------------------------------------
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["DejaVu Serif"],
    "font.size": 9,
    "axes.labelsize": 9,
    "axes.titlesize": 9,
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linestyle": ":",
    "figure.dpi": 160,
    "savefig.bbox": "tight",
})
C = dict(blue="#1f4e79", red="#c00000", green="#2e7d32", orange="#e07b00",
         purple="#6a3d9a", gray="#555555")

# ---------------------------------------------------------------------------
# Parameter desain acuan (hasil konsolidasi)
# ---------------------------------------------------------------------------
T_DESIGN = 1075.0          # K, temperatur uap di evaporator (titik desain)
T_COND = 1073.0            # K, temperatur uap di kondensor (selisih ~2 K)
N_PIN = 61
P_CORE = 500e3             # W
FR = np.array([1.28, 1.21, 1.09, 0.93, 0.81])
N_RING = np.array([1, 6, 12, 18, 24])
RING_NAME = ["Cincin 1 (pusat)", "Cincin 2", "Cincin 3", "Cincin 4", "Cincin 5 (tepi)"]
Q_AVG = P_CORE / N_PIN
_norm = float((FR * N_RING).sum() / N_PIN)
Q_PEAK = Q_AVG * FR[0] / _norm
PIN = dict(r_o=16.3e-3, H=0.250, k_UN=16.0, rho_UN=13.55e3, rho_frac=0.95)
LIMIT_FUEL = 1650.0        # K, batas desain konservatif integritas UN
HP = m.isr_pipe(OD=0.018, t_wall=0.0009, t_wick=0.0008, K=5.0e-9, r_eff=1.5e-5)

summary = {}
summary["Q_avg"] = Q_AVG
summary["Q_peak"] = Q_PEAK
summary["f_peak"] = FR[0] / _norm
summary["D_v_mm"] = HP.D_v * 1e3


# ===========================================================================
# 1. TABEL SIFAT TERMOFISIKA
# ===========================================================================
def study_properties():
    rows = []
    for T in range(800, 1226, 25):
        p = m.properties(T)
        rows.append(dict(T=T, P_sat_kPa=p["P_sat"] / 1e3, rho_l=p["rho_l"],
                         rho_v=p["rho_v"], mu_l=p["mu_l"] * 1e6,
                         mu_v=p["mu_v"] * 1e6, k_l=p["k_l"],
                         cp_l=p["cp_l"], sigma=p["sigma"] * 1e3,
                         h_fg=p["h_fg"] / 1e3, c_s=p["c_s"]))
    df = pd.DataFrame(rows)
    df.to_csv(f"{RES}/tabel_sifat.csv", index=False)
    return df


# ===========================================================================
# 2. VERIFIKASI NUMERIK
# ===========================================================================
def study_verification():
    rows = []
    # (a) konvergensi grid
    ref = m.solve_axial(HP, Q_PEAK, T_COND, N=1920)
    for N in (30, 60, 120, 240, 480):
        s = m.solve_axial(HP, Q_PEAK, T_COND, N=N)
        rows.append(dict(uji="konvergensi_grid", parameter=f"N={N}",
                         nilai=s.dP_v, referensi=ref.dP_v,
                         galat_rel=abs(s.dP_v - ref.dP_v) / ref.dP_v * 100))
    # (b) neraca energi: m_dot * h_fg = Q
    s = m.solve_axial(HP, Q_PEAK, T_COND, N=240)
    mdot = s.Q / m.h_fg(s.T_cond)
    rows.append(dict(uji="neraca_energi", parameter="mdot*hfg,Q",
                     nilai=mdot * m.h_fg(s.T_cond), referensi=Q_PEAK,
                     galat_rel=abs(mdot * m.h_fg(s.T_cond) - Q_PEAK) / Q_PEAK * 100))
    # (c) Mach maksimum pada titik desain (informasi, kriteria < 1)
    s_des = m.solve_axial(HP, Q_PEAK, T_COND, N=240)
    rows.append(dict(uji="Mach_desain", parameter="Ma_maks pada Q_puncak",
                     nilai=float(s_des.Ma.max()), referensi=1.0,
                     galat_rel=float("nan")))
    # (d) aliran cairan: solver vs Darcy analitik pada titik desain
    T = T_DESIGN
    F_l = m.mu_l(T) / (m.rho_l(T) * HP.A_w * HP.K * m.h_fg(T))
    dP_l_an = Q_PEAK * F_l * HP.L_eff
    rows.append(dict(uji="tekanan_cairan", parameter="solver vs Darcy analitik",
                     nilai=s_des.dP_l, referensi=dP_l_an,
                     galat_rel=abs(s_des.dP_l - dP_l_an) / dP_l_an * 100))
    # (e) aliran uap bermassa konstan vs Hagen-Poiseuille
    import dataclasses
    hp_c = dataclasses.replace(HP, L_e=2e-4, L_a=0.29, L_c=2e-4)
    Qc = 50.0
    s_c = m.solve_axial(hp_c, Qc, T_DESIGN, N=400)
    rho_c = m.rho_v(T_DESIGN)
    u_c = Qc / (rho_c * hp_c.A_v * m.h_fg(T_DESIGN))
    dP_hp = 32.0 * m.mu_v(T_DESIGN) * u_c * hp_c.L / hp_c.D_v ** 2
    rows.append(dict(uji="tekanan_uap", parameter="solver vs Hagen-Poiseuille",
                     nilai=s_c.dP_v, referensi=dP_hp,
                     galat_rel=abs(s_c.dP_v - dP_hp) / dP_hp * 100))
    # (e) jangkar sifat termofisika
    anchors = [("T_didih pada 101,3 kPa [K]", m.tsat(101325.0), 1156.2),
               ("rho_l(1000 K) [kg/m3]", m.rho_l(1000.0), 780.0),
               ("sigma(1000 K) [N/m]", m.sigma(1000.0), 0.135),
               ("h_fg(1000 K) [MJ/kg]", m.h_fg(1000.0) / 1e6, 3.9),
               ("P_sat(1000 K) [kPa]", m.psat(1000.0) / 1e3, 24.0)]
    for name, val, refv in anchors:
        rows.append(dict(uji="jangkar_sifat", parameter=name, nilai=val,
                         referensi=refv, galat_rel=abs(val - refv) / refv * 100))
    df = pd.DataFrame(rows)
    df.to_csv(f"{RES}/verifikasi.csv", index=False)
    return df


# ===========================================================================
# 3. VALIDASI KRUSTY
# ===========================================================================
KRUSTY_DATA = [
    # Q [W], T_evap_exp [C], T_cond_exp [C]
    (312.5, 718.0, 692.0),
    (450.0, 752.0, 720.0),
    (625.0, 795.0, 758.0),
    (820.0, 834.0, 788.0),
    (1000.0, 868.0, 812.0),
]


def study_krusty():
    hp = m.krusty_pipe()
    rows = []
    for i, (Q, Te, Tc) in enumerate(KRUSTY_DATA):
        TeK, TcK = Te + 273.15, Tc + 273.15
        s = m.solve_axial(hp, Q, TcK, N=240)
        Rth = m.thermal_resistance(hp, 0.5 * (TeK + TcK))
        dT_hp = (s.T_evap - s.T_cond) + Q * Rth["R_wall"] + Q * Rth["R_wick"]
        rows.append(dict(run=f"KRUSTY-SS-{i+1}", Q=Q, T_evap_exp=TeK,
                         T_cond_exp=TcK, dT_exp=TeK - TcK, dT_hp=dT_hp,
                         Rth_exp=(TeK - TcK) / Q))
    df = pd.DataFrame(rows)

    # --- kalibrasi dua konstanta (R_int, b) dengan kuadrat terkecil
    #     DeltaT_pred = b + a*Q + dT_hp(Q)
    def fit(Q, resid, w=None):
        M = np.column_stack([Q, np.ones_like(Q)])
        coef, *_ = np.linalg.lstsq(M, resid, rcond=None)
        return coef[0], coef[1]

    a_all, b_all = fit(df["Q"].values, (df["dT_exp"] - df["dT_hp"]).values)
    # --- leave-one-out: latih pada 4 titik, uji pada titik kelima
    pred = []
    for i in range(len(df)):
        tr = df.drop(index=i)
        a, b = fit(tr["Q"].values, (tr["dT_exp"] - tr["dT_hp"]).values)
        pred.append(b + a * df.loc[i, "Q"] + df.loc[i, "dT_hp"])
    df["dT_pred"] = b_all + a_all * df["Q"] + df["dT_hp"]
    df["dT_loo"] = pred
    df["galat_rel"] = (df["dT_loo"] - df["dT_exp"]).abs() / df["dT_exp"] * 100
    df["Rth_meas"] = df["dT_exp"] / df["Q"]
    df["Rth_hp"] = df["dT_hp"] / df["Q"]
    df.to_csv(f"{RES}/krusty_validasi.csv", index=False)

    met = dict(R_int=float(a_all), b_offset=float(b_all),
               MAPE=float(df["galat_rel"].mean()),
               RMSE=float(np.sqrt(((df["dT_loo"] - df["dT_exp"]) ** 2).mean())),
               MAE=float((df["dT_loo"] - df["dT_exp"]).abs().mean()),
               max_err=float(df["galat_rel"].max()))
    y, yh = df["dT_exp"].values, df["dT_loo"].values
    met["R2"] = float(1 - ((y - yh) ** 2).sum() / ((y - y.mean()) ** 2).sum())
    met["Rth_hp_mean"] = float(df["Rth_hp"].mean())
    summary["krusty"] = met
    return df, met


# ===========================================================================
# 4. AMPLOP BATAS OPERASI
# ===========================================================================
def study_envelope():
    rows = []
    for g, tag in ((0.0, "0g"), (9.81, "1g")):
        for T in range(900, 1176, 25):
            e = m.envelope(HP, T, g=g, angle_deg=90.0)
            rows.append(dict(kondisi=tag, T=T, Q_cap=e.Q_cap, Q_sonic=e.Q_sonic,
                             Q_ent=e.Q_ent, Q_boil=e.Q_boil, Q_visc=e.Q_visc,
                             Q_env=e.Q_min, pembatas=e.limiting))
    df = pd.DataFrame(rows)
    df.to_csv(f"{RES}/amplop.csv", index=False)
    return df


# ===========================================================================
# 5. PROFIL AKSIAL TITIK DESAIN
# ===========================================================================
def study_axial():
    s = m.solve_axial(HP, Q_PEAK, T_COND, N=240)
    Rth = m.thermal_resistance(HP, 0.5 * (s.T_evap + s.T_cond))
    T_wall = s.T_v + Q_PEAK * (Rth["R_wick"] + Rth["R_wall"]) * (
        np.where(s.z <= HP.L_e + HP.L_a, 1.0, 0.5))
    df = pd.DataFrame(dict(z=s.z, P_v_kPa=s.P_v / 1e3, P_l_kPa=s.P_l / 1e3,
                           T_v=s.T_v, T_wall=T_wall, u_v=s.u_v, Ma=s.Ma,
                           r_c_um=s.r_c * 1e6))
    df.to_csv(f"{RES}/profil_aksial.csv", index=False)
    summary["axial"] = dict(T_evap=s.T_evap, T_cond=s.T_cond,
                            dP_v=s.dP_v, dP_l=s.dP_l, dP_g=s.dP_g,
                            dP_req=s.dP_req, dP_cap=m.dP_cap_max(HP, s.T_evap),
                            Ma_max=float(s.Ma.max()), u_max=float(s.u_v.max()),
                            r_c_min_um=float(s.r_c.min() * 1e6),
                            R_wick=Rth["R_wick"], R_wall=Rth["R_wall"],
                            R_vap=Rth["R_vap"])
    return s


# ===========================================================================
# 6. STUDI GEOMETRI (TRADE-OFF)
# ===========================================================================
def study_geometry():
    rows = []
    r_fuel_o = PIN["r_o"]
    hole0 = 0.00765 / 2
    for OD in np.arange(0.010, 0.0241, 0.002):
        hp = m.isr_pipe(OD=float(OD), t_wall=0.0009, t_wick=0.0008,
                        K=5.0e-9, r_eff=1.5e-5)
        e = m.envelope(hp, T_DESIGN)
        s = m.solve_axial(hp, Q_PEAK, T_COND, N=120)
        fuel_loss = 1 - (r_fuel_o ** 2 - (OD / 2) ** 2) / (r_fuel_o ** 2 - hole0 ** 2)
        rows.append(dict(OD_mm=OD * 1e3, D_v_mm=hp.D_v * 1e3,
                         Q_cap=e.Q_cap / 1e3, Q_sonic=e.Q_sonic / 1e3,
                         Q_ent=e.Q_ent / 1e3, Q_env=e.Q_min / 1e3,
                         pembatas=e.limiting,
                         SF=e.Q_min / Q_PEAK,
                         dP_req=s.dP_req, Ma=s.Ma.max(),
                         fuel_loss=fuel_loss * 100,
                         merit=(e.Q_min / Q_PEAK) / max(fuel_loss, 1e-6)))
    df = pd.DataFrame(rows)
    df.to_csv(f"{RES}/trade_geometri.csv", index=False)
    return df


# ===========================================================================
# 7. INTEGRASI TERAS 61 PIN
# ===========================================================================
def study_core():
    rows = []
    T_op = T_DESIGN
    for i, (fr, n) in enumerate(zip(FR, N_RING)):
        Q = Q_AVG * fr / _norm
        e = m.envelope(HP, T_op)
        sf = e.Q_min / Q
        s = m.solve_axial(HP, Q, T_COND, N=160)
        Rth = m.thermal_resistance(HP, T_op)
        T_clad = s.T_evap + Q * (Rth["R_wick"] + Rth["R_wall"])
        r_i = HP.r_o + 0.0001
        fd = m.annular_fuel_dT(Q, r_i, PIN["r_o"], PIN["H"], PIN["k_UN"])
        rows.append(dict(cincin=RING_NAME[i], n=n, fr=fr, Q=Q, Q_tot=Q * n,
                         T_evap=s.T_evap, T_clad=T_clad,
                         T_fuel=T_clad + fd["dT"], dT_fuel=fd["dT"],
                         q3=fd["q3"], dP_req=s.dP_req, Ma=s.Ma.max(),
                         margin_cap=(m.dP_cap_max(HP, T_op) - s.dP_req)
                         / m.dP_cap_max(HP, T_op) * 100,
                         Q_env=e.Q_min, SF=sf,
                         SF_cap=m.q_capillary(HP, T_op) / Q,
                         SF_sonic=e.Q_sonic / Q, SF_ent=e.Q_ent / Q))
    df = pd.DataFrame(rows)
    df.to_csv(f"{RES}/teras.csv", index=False)
    summary["core"] = dict(
        T_clad_max=float(df["T_clad"].max()),
        T_fuel_max=float(df["T_fuel"].max()),
        SF_min=float(df["SF"].min()),
        margin_cap_min=float(df["margin_cap"].min()),
        Q_tot=float(df["Q_tot"].sum()))
    return df


# ===========================================================================
# 8. TOLERANSI KEGAGALAN SATU PIPA
# ===========================================================================
def study_failure():
    rows = []
    h_gap = 2000.0          # W/(m2 K) konduktansi celah efektif
    A_lat = 0.020           # m2 luas kontak lateral efektif per pin
    R_lat = 1.0 / (h_gap * A_lat)
    for n_fail in (0, 1, 2):
        fr_eff = FR.copy()
        # beban pin gagal dibagi ke enam tetangga pada cincin yang sama
        q_fail = Q_AVG * FR[0] / _norm * n_fail
        q_neigh = q_fail * n_fail / 6.0 / 6.0 if n_fail else 0.0
        rows.append(dict(
            n_gagal=n_fail,
            Q_tetangga_tambahan=q_neigh,
            dT_fuel_gagal=R_lat * q_fail if n_fail else 0.0,
            catatan=("nominal" if n_fail == 0 else
                     f"{n_fail} pipa gagal; beban dibagi ke 6 pin tetangga")))
    df = pd.DataFrame(rows)
    df.to_csv(f"{RES}/kegagalan.csv", index=False)
    # hitung eksplisit kasus 1 pipa
    q1 = Q_PEAK
    dT_rise = R_lat * q1
    summary["failure"] = dict(R_lat=R_lat, dT_rise=dT_rise,
                              T_fuel_failed=summary["core"]["T_fuel_max"] + dT_rise,
                              sf_neighbor=summary["core"]["SF_min"] / 1.167)
    return df


# ===========================================================================
# 9. NERACA MASSA DAN DAYA PARASITIK
# ===========================================================================
def study_mass():
    # --- baseline loop Na-78K
    base = [
        ("Pompa elektromagnetik (2 unit, 1 cadangan)", 19.0),
        ("Pengondisian daya pompa & kabel", 5.5),
        ("Tangki ekspansi & gas penutup He", 12.0),
        ("Perpipaan & manifold tekanan", 15.0),
        ("Inventaris Na-78K (14 L)", 10.5),
        ("Instrumentasi aliran & katup", 3.0),
    ]
    retr = [
        ("61 pipa kalor Mo-30Re (OD 18 mm, t 0,9 mm)", 61 * 0.253),
        ("Struktur kasa & arteri komposit", 61 * 0.032),
        ("Inventaris natrium (61 pipa)", 61 * 0.011),
        ("Blok sadel pipa-kondensor & header", 6.0),
    ]
    mb = sum(v for _, v in base)
    mr = sum(v for _, v in retr)
    rows = []
    for k, v in base:
        rows.append(dict(kategori="Baseline Na-78K", item=k, massa=v))
    for k, v in retr:
        rows.append(dict(kategori="Retrofit pipa kalor", item=k, massa=v))
    rows.append(dict(kategori="Baseline Na-78K", item="SUBTOTAL", massa=mb))
    rows.append(dict(kategori="Retrofit pipa kalor", item="SUBTOTAL", massa=mr))
    df = pd.DataFrame(rows)
    df.to_csv(f"{RES}/neraca_massa.csv", index=False)

    # --- daya pompa
    cp_nak = 1200.0
    dT_loop = 200.0
    mdot = P_CORE / (cp_nak * dT_loop)
    dP_loop = 30000.0
    P_hyd = mdot * dP_loop / 780.0
    eta = 0.15
    rows2 = [dict(besaran="Laju alir massa Na-78K", nilai=mdot, satuan="kg/s"),
             dict(besaran="Tekanan pompa total", nilai=dP_loop / 1e3, satuan="kPa"),
             dict(besaran="Daya hidraulik", nilai=P_hyd, satuan="W"),
             dict(besaran="Efisiensi pompa EM", nilai=eta * 100, satuan="%"),
             dict(besaran="Daya listrik pompa", nilai=P_hyd / eta, satuan="W"),
             dict(besaran="Fraksi daya termal", nilai=P_hyd / eta / P_CORE * 100,
                  satuan="%")]
    pd.DataFrame(rows2).to_csv(f"{RES}/daya_pompa.csv", index=False)

    # --- massa bahan bakar & umur
    hole0 = 0.00765 / 2
    r_i = HP.r_o + 0.0001
    V0 = math.pi * (PIN["r_o"] ** 2 - hole0 ** 2) * PIN["H"]
    V1 = math.pi * (PIN["r_o"] ** 2 - r_i ** 2) * PIN["H"]
    rho_eff = PIN["rho_UN"] * PIN["rho_frac"]
    m0, m1 = V0 * rho_eff * N_PIN, V1 * rho_eff * N_PIN
    life0 = 17.0
    rows3 = [dict(besaran="Volume UN per pin baseline", nilai=V0 * 1e6, satuan="cm3"),
             dict(besaran="Volume UN per pin retrofit", nilai=V1 * 1e6, satuan="cm3"),
             dict(besaran="Massa UN baseline (61 pin)", nilai=m0, satuan="kg"),
             dict(besaran="Massa UN retrofit (61 pin)", nilai=m1, satuan="kg"),
             dict(besaran="Kehilangan UN", nilai=(m0 - m1) / m0 * 100, satuan="%"),
             dict(besaran="Umur teras ekuivalen", nilai=life0 * m1 / m0, satuan="tahun"),
             dict(besaran="Tinggi aktif untuk paritas volume",
                  nilai=PIN["H"] * (V0 / V1) * 1e3, satuan="mm")]
    pd.DataFrame(rows3).to_csv(f"{RES}/bahan_bakar.csv", index=False)

    summary["mass"] = dict(baseline=mb, retrofit=mr, saving=mb - mr,
                           saving_pct=(mb - mr) / mb * 100,
                           fuel0=m0, fuel1=m1, fuel_loss=(m0 - m1) / m0 * 100,
                           life=life0 * m1 / m0,
                           P_pump=P_hyd / eta)
    return df


# ===========================================================================
# 10. SENSITIVITAS
# ===========================================================================
def study_sensitivity():
    rows = []
    for T in range(1000, 1151, 25):
        e = m.envelope(HP, T)
        s = m.solve_axial(HP, Q_PEAK, T - 25, N=120)
        Rth = m.thermal_resistance(HP, T)
        T_clad = s.T_evap + Q_PEAK * (Rth["R_wick"] + Rth["R_wall"])
        fd = m.annular_fuel_dT(Q_PEAK, HP.r_o + 0.0001, PIN["r_o"], PIN["H"],
                               PIN["k_UN"])
        rows.append(dict(T=T, Q_env=e.Q_min, SF=e.Q_min / Q_PEAK,
                         Ma=s.Ma.max(), dP_req=s.dP_req,
                         T_fuel=T_clad + fd["dT"]))
    df = pd.DataFrame(rows)
    df.to_csv(f"{RES}/sensitivitas.csv", index=False)
    return df


# ===========================================================================
# GAMBAR
# ===========================================================================
def figures(df_prop, df_krusty, df_env, s_ax, df_geo, df_core, df_mass, df_sens):
    # --- Gambar 3.1: sifat termofisika
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.7))
    T = df_prop["T"]
    ax[0].semilogy(T, df_prop["P_sat_kPa"], color=C["blue"], label=r"$P_{sat}$ (kPa)")
    ax[0].semilogy(T, df_prop["rho_v"] * 1e3, color=C["red"], ls="--",
                   label=r"$\rho_v\times10^{3}$ (kg/m$^3$)")
    ax[0].set_xlabel("Temperatur (K)"); ax[0].set_ylabel("Nilai")
    ax[0].legend(loc="lower right"); ax[0].set_title("(a) Tekanan dan densitas uap")
    ax[1].plot(T, df_prop["sigma"], color=C["green"], label=r"$\sigma$ (mN/m)")
    ax[1].plot(T, df_prop["k_l"], color=C["purple"], ls="-.",
               label=r"$k_l$ (W/m K)")
    ax[1].set_xlabel("Temperatur (K)"); ax[1].set_ylim(40, 180)
    axb = ax[1].twinx()
    axb.plot(T, df_prop["h_fg"] / 1e3, color=C["orange"], ls="--",
             label=r"$h_{fg}$ (MJ/kg)")
    axb.set_ylabel("Panas laten (MJ/kg)"); axb.set_ylim(3.0, 4.4)
    h1, l1 = ax[1].get_legend_handles_labels()
    h2, l2 = axb.get_legend_handles_labels()
    ax[1].legend(h1 + h2, l1 + l2, loc="lower left", fontsize=7)
    ax[1].set_title("(b) Tegangan permukaan, konduktivitas, panas laten")
    fig.tight_layout(); fig.savefig(f"{FIG}/sifat_natrium.pdf"); plt.close(fig)

    # --- Gambar 5.1: validasi KRUSTY + profil aksial
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.9))
    ax[0].plot(df_krusty["Q"], df_krusty["dT_exp"], "o", color=C["blue"],
               ms=6, label="Eksperimen KRUSTY")
    ax[0].plot(df_krusty["Q"], df_krusty["dT_loo"], "s--", color=C["red"],
               ms=4, label="Model ROM (leave-one-out)")
    lo = min(df_krusty["dT_exp"].min(), df_krusty["dT_loo"].min()) - 3
    hi = max(df_krusty["dT_exp"].max(), df_krusty["dT_loo"].max()) + 3
    ax[0].plot([lo, hi], [lo, hi], ":", color=C["gray"], lw=1)
    ax[0].set_xlabel("Daya termal per pipa (W)")
    ax[0].set_ylabel(r"$\Delta T$ evaporator-kondensor (K)")
    ax[0].legend(); ax[0].set_title("(a) Karakteristik $\\Delta T$ vs daya")

    ax[1].plot([20, 60], [20, 60], ":", color=C["gray"], lw=1)
    ax[1].plot(df_krusty["dT_exp"], df_krusty["dT_loo"], "o", color=C["blue"],
               ms=6, label="Titik uji (LOO)")
    ax[1].fill_between([20, 60], [16, 56], [24, 64], color=C["red"], alpha=0.10,
                       label="pita $\\pm10\\%$")
    ax[1].set_xlim(20, 60); ax[1].set_ylim(20, 60)
    ax[1].set_xlabel(r"$\Delta T$ terukur (K)"); ax[1].set_ylabel(r"$\Delta T$ model (K)")
    ax[1].set_title("(b) Diagram paritas validasi")
    ax[1].legend(loc="lower right")
    fig.tight_layout(); fig.savefig(f"{FIG}/validasi_profil.pdf"); plt.close(fig)

    # --- Gambar 5.2: amplop operasi
    fig, ax = plt.subplots(figsize=(5.6, 3.4))
    g0 = df_env[df_env["kondisi"] == "0g"]
    ax.plot(g0["T"], g0["Q_cap"] / 1e3, color=C["blue"], label=r"$Q_{cap}$ (0g)")
    ax.plot(g0["T"], g0["Q_sonic"] / 1e3, color=C["red"], ls="--", label=r"$Q_{sonic}$")
    ax.plot(g0["T"], g0["Q_ent"] / 1e3, color=C["green"], ls="-.", label=r"$Q_{ent}$")
    ax.plot(g0["T"], g0["Q_boil"] / 1e3, color=C["orange"], ls=":", label=r"$Q_{boil}$")
    g1 = df_env[df_env["kondisi"] == "1g"]
    ax.plot(g1["T"], g1["Q_cap"] / 1e3, color=C["gray"], ls="--", lw=1.0,
            label=r"$Q_{cap}$ (1g, lawan gravitasi)")
    ax.axhline(Q_PEAK / 1e3, color=C["purple"], lw=1.2, ls="-")
    ax.annotate(f"Beban puncak pin = {Q_PEAK/1e3:.2f} kW", (1040, Q_PEAK / 1e3 * 0.55),
                color=C["purple"], fontsize=8)
    ax.set_yscale("log"); ax.set_xlabel("Temperatur uap (K)")
    ax.set_ylabel("Daya transpor (kW)"); ax.legend(ncol=2, loc="upper left")
    fig.tight_layout(); fig.savefig(f"{FIG}/amplop_operasi.pdf"); plt.close(fig)

    # --- Gambar 5.3: profil tekanan, Mach, meniskus
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.5))
    ax[0].plot(s_ax.z, s_ax.P_v / 1e3, color=C["blue"], label=r"$P_v$")
    ax[0].plot(s_ax.z, s_ax.P_l / 1e3, color=C["red"], ls="--", label=r"$P_l$")
    ax[0].set_xlabel("z (m)"); ax[0].set_ylabel("Tekanan (kPa)"); ax[0].legend()
    ax[0].set_title("(a) Tekanan uap dan cairan")
    ax[1].plot(s_ax.z, s_ax.Ma, color=C["green"])
    ax[1].set_ylim(0, 0.25); ax[1].set_xlabel("z (m)")
    ax[1].set_ylabel("Bilangan Mach"); ax[1].set_title("(b) Mach aksial")
    ax[1].annotate(f"$Ma_{{maks}}$ = {s_ax.Ma.max():.3f}", xy=(0.26, 0.17),
                   fontsize=7, color=C["green"])
    dPc = (s_ax.P_v - s_ax.P_l) / 1e3
    ax[2].plot(s_ax.z, dPc, color=C["purple"], label=r"$\Delta P_c(z)$")
    ax[2].axhline(m.dP_cap_max(HP, T_DESIGN) / 1e3, color=C["red"], ls=":", lw=1,
                  label=r"$\Delta P_{c,maks}$")
    ax[2].set_ylim(0, m.dP_cap_max(HP, T_DESIGN) / 1e3 * 1.15)
    ax[2].set_xlabel("z (m)"); ax[2].set_ylabel(r"$\Delta P_c$ (kPa)")
    ax[2].set_title("(c) Tekanan kapiler"); ax[2].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(f"{FIG}/profil_detail.pdf"); plt.close(fig)

    # --- Gambar 5.4: trade-off geometri
    fig, ax = plt.subplots(2, 1, figsize=(5.4, 4.6), sharex=True,
                           gridspec_kw=dict(height_ratios=[2.1, 1.0]))
    a0 = ax[0]
    a0.plot(df_geo["OD_mm"], df_geo["Q_cap"], "o-", color=C["blue"],
            ms=3, label=r"$Q_{cap}$")
    a0.plot(df_geo["OD_mm"], df_geo["Q_sonic"], "s--", color=C["red"],
            ms=3, label=r"$Q_{sonic}$")
    a0.plot(df_geo["OD_mm"], df_geo["Q_ent"], "^-.", color=C["green"],
            ms=3, label=r"$Q_{ent}$")
    a0.plot(df_geo["OD_mm"], df_geo["Q_env"], "-", color=C["gray"], lw=2.2,
            alpha=0.45, label="Amplop")
    a0.axhline(Q_PEAK / 1e3, color=C["purple"], ls=":", lw=1.2)
    a0.text(10.2, Q_PEAK / 1e3 * 1.12, "beban puncak pin", color=C["purple"],
            fontsize=7)
    a0.set_yscale("log"); a0.set_ylabel("Daya (kW)")
    a0.set_ylim(2, 200); a0.legend(ncol=2, fontsize=7, loc="lower right")
    a0.axvspan(17.0, 19.0, color=C["blue"], alpha=0.07)
    a1 = ax[1]
    a1.plot(df_geo["OD_mm"], df_geo["fuel_loss"], "d-", color=C["orange"], ms=3,
            label="Kehilangan volume UN (%)")
    a1.plot(df_geo["OD_mm"], df_geo["SF"], "s-", color=C["blue"], ms=3,
            label="SF amplop (-)")
    a1.axhline(1.3, color=C["red"], ls=":", lw=1)
    a1.axvspan(17.0, 19.0, color=C["blue"], alpha=0.07)
    a1.set_xlabel("Diameter luar pipa (mm)")
    a1.legend(fontsize=7, loc="upper left")
    fig.tight_layout(); fig.savefig(f"{FIG}/trade_geometri.pdf"); plt.close(fig)

    # --- Gambar 5.5: margin teras per cincin
    fig, ax = plt.subplots(figsize=(5.6, 3.1))
    x = np.arange(len(df_core))
    ax.bar(x - 0.27, df_core["SF_cap"], 0.26, color=C["blue"], label=r"SF kapiler")
    ax.bar(x, df_core["SF_sonic"], 0.26, color=C["red"], label=r"SF sonik")
    ax.bar(x + 0.27, df_core["SF_ent"], 0.26, color=C["green"],
           label=r"SF entrainment")
    ax.axhline(1.0, color=C["gray"], ls="--", lw=1)
    ax.set_xticks(x); ax.set_xticklabels([f"C{i+1}" for i in x])
    ax.set_xlabel("Cincin teras"); ax.set_ylabel("Faktor keselamatan (-)")
    ax.legend(ncol=3, fontsize=7); ax.set_ylim(0, df_core["SF_cap"].max() * 1.2)
    fig.tight_layout(); fig.savefig(f"{FIG}/margin_teras.pdf"); plt.close(fig)

    # --- Gambar 5.6: neraca massa
    fig, ax = plt.subplots(figsize=(5.4, 3.0))
    mb = df_mass[df_mass["item"] == "SUBTOTAL"]
    ax.bar(["Baseline\nNa-78K", "Retrofit\npipa kalor"],
           mb["massa"].values, color=[C["gray"], C["blue"]], width=0.5)
    for i, v in enumerate(mb["massa"].values):
        ax.text(i, v + 1, f"{v:.1f} kg", ha="center", fontsize=9)
    ax.set_ylabel("Massa subsistem pendingin (kg)")
    fig.tight_layout(); fig.savefig(f"{FIG}/neraca_massa.pdf"); plt.close(fig)

    # --- Gambar 5.7: sensitivitas
    fig, ax = plt.subplots(figsize=(5.6, 3.1))
    ax.plot(df_sens["T"], df_sens["SF"], "o-", color=C["blue"], label="SF amplop")
    ax.plot(df_sens["T"], df_sens["Ma"], "s--", color=C["red"], label="Mach maks.")
    ax.axhline(1.0, color=C["gray"], ls=":", lw=1)
    ax.set_xlabel("Temperatur uap evaporator (K)")
    ax.set_ylabel("SF (-) / Mach (-)"); ax.legend()
    ax2 = ax.twinx()
    ax2.plot(df_sens["T"], df_sens["T_fuel"], "^-", color=C["orange"],
             label=r"$T_{fuel,max}$")
    ax2.axhline(LIMIT_FUEL, color=C["orange"], ls=":", lw=1)
    ax2.set_ylabel("Temperatur puncak UN (K)")
    fig.tight_layout(); fig.savefig(f"{FIG}/sensitivitas.pdf"); plt.close(fig)




# ===========================================================================
# GAMBAR SKEMATIK DAN DIAGRAM ALIR
# ===========================================================================
def draw_schematic():
    from matplotlib.patches import Rectangle, FancyArrowPatch
    fig, ax = plt.subplots(figsize=(7.0, 2.1))
    y0, h = 0.0, 0.30
    L = 1.0
    le, la, lc = 0.25, 0.30, 0.45
    # selubung
    ax.add_patch(Rectangle((0, 0), L, h, facecolor="#d9e2f3",
                           edgecolor="black", lw=1.2, zorder=1))
    # sumbu (annulus)
    ax.add_patch(Rectangle((0, 0.055), L, 0.035, facecolor="#f4b183",
                           edgecolor="black", lw=0.8, zorder=2))
    ax.add_patch(Rectangle((0, h - 0.09), L, 0.035, facecolor="#f4b183",
                           edgecolor="black", lw=0.8, zorder=2))
    # inti uap
    ax.add_patch(Rectangle((0, 0.11), L, h - 0.22, facecolor="white",
                           edgecolor="black", lw=0.8, zorder=2))
    # batas zona
    for x in (le, le + la):
        ax.plot([x, x], [-0.02, h + 0.02], color=C["gray"], ls="--", lw=0.9)
    ax.text(le / 2, h + 0.055, "Evaporator", ha="center", fontsize=8)
    ax.text(le + la / 2, h + 0.055, "Adiabatik", ha="center", fontsize=8)
    ax.text(le + la + lc / 2, h + 0.055, "Kondensor", ha="center", fontsize=8)
    # aliran
    ax.annotate("", xy=(0.55, h - 0.02), xytext=(0.44, h - 0.02),
                arrowprops=dict(arrowstyle="-|>", color=C["red"], lw=1.4))
    ax.text(0.49, h + 0.005, r"uap Na ($\dot m_v$)", color=C["red"],
            fontsize=7, ha="center")
    ax.annotate("", xy=(0.40, 0.075), xytext=(0.30, 0.075),
                arrowprops=dict(arrowstyle="-|>", color=C["orange"], lw=1.4))
    ax.text(0.35, 0.10, r"cair ($\dot m_l$)", color=C["orange"], fontsize=7,
            ha="center")
    ax.annotate("", xy=(0.90, 0.075), xytext=(0.80, 0.075),
                arrowprops=dict(arrowstyle="-|>", color=C["orange"], lw=1.4))
    # kalor masuk / keluar
    ax.annotate("", xy=(le / 2, 0.005), xytext=(le / 2, -0.09),
                arrowprops=dict(arrowstyle="-|>", color=C["blue"], lw=1.6))
    ax.text(le / 2, -0.12, r"$Q_{in}$", color=C["blue"], ha="center", fontsize=8)
    ax.annotate("", xy=(le + la + lc / 2, -0.09),
                xytext=(le + la + lc / 2, 0.005),
                arrowprops=dict(arrowstyle="-|>", color=C["blue"], lw=1.6))
    ax.text(le + la + lc / 2, -0.12, r"$Q_{out}$", color=C["blue"], ha="center",
            fontsize=8)
    # keterangan radial
    ax.annotate("selubung (Mo-30Re)", xy=(0.02, 0.0), xytext=(0.02, -0.06),
                fontsize=6.5, color=C["gray"],
                arrowprops=dict(arrowstyle="-", color=C["gray"], lw=0.6))
    ax.annotate("sumbu komposit\n(arteri + kasa)", xy=(0.70, 0.055),
                xytext=(0.70, -0.075), fontsize=6.5, color="#b05a00",
                arrowprops=dict(arrowstyle="-", color="#b05a00", lw=0.6))
    ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.16, h + 0.10)
    ax.axis("off"); ax.grid(False)
    fig.savefig(f"{FIG}/skema_pipa_kalor.pdf"); plt.close(fig)


def draw_flowchart():
    from matplotlib.patches import FancyBboxPatch
    fig, ax = plt.subplots(figsize=(4.6, 6.4))
    steps = [
        ("Studi literatur & konsolidasi\ndata sekunder", "#d9e2f3"),
        ("Formulasi ROM 1D dua fasa\n(uap + sumbu + jaringan termal)", "#fbe5d6"),
        ("Verifikasi numerik\n(grid, neraca energi, sifat)", "#e2efda"),
        ("Kalibrasi 2 konstanta antarmuka\n& validasi leave-one-out (KRUSTY)", "#fff2cc"),
        ("Pemetaan amplop batas operasi\n(0g vs 1g, 900-1200 K)", "#e2efda"),
        ("Optimasi geometri pipa\n(OD vs margin vs volume UN)", "#fbe5d6"),
        ("Integrasi teras ISR 61 pin\n(daya, suhu, margin)", "#e2efda"),
        ("Analisis massa, daya parasitik,\n& toleransi kegagalan", "#fbe5d6"),
        ("Kesimpulan & rekomendasi desain", "#d9e2f3"),
    ]
    y = 1.0
    dy = 0.105
    for i, (txt, col) in enumerate(steps):
        yy = y - i * dy
        ax.add_patch(FancyBboxPatch((0.08, yy - 0.032), 0.84, 0.062,
                                    boxstyle="round,pad=0.006,rounding_size=0.01",
                                    facecolor=col, edgecolor="black", lw=0.7))
        ax.text(0.5, yy, txt, ha="center", va="center", fontsize=6.8)
        if i < len(steps) - 1:
            ax.annotate("", xy=(0.5, yy - dy + 0.036), xytext=(0.5, yy - 0.036),
                        arrowprops=dict(arrowstyle="-|>", color="black", lw=0.8))
    ax.set_xlim(0, 1); ax.set_ylim(0.02, 1.03)
    ax.axis("off"); ax.grid(False)
    fig.savefig(f"{FIG}/diagram_alir.pdf"); plt.close(fig)


# ===========================================================================
def main():
    print("== sifat ==");      df_prop = study_properties()
    print("== verifikasi =="); df_ver = study_verification()
    print("== KRUSTY ==");     df_kr, met = study_krusty()
    print("== amplop ==");     df_env = study_envelope()
    print("== aksial ==");     s_ax = study_axial()
    print("== geometri ==");   df_geo = study_geometry()
    print("== teras ==");      df_core = study_core()
    print("== kegagalan ==");  study_failure()
    print("== massa ==");      df_mass = study_mass()
    print("== sensitivitas =="); df_sens = study_sensitivity()
    print("== gambar ==")
    figures(df_prop, df_kr, df_env, s_ax, df_geo, df_core, df_mass, df_sens)
    draw_schematic()
    draw_flowchart()
    with open(f"{RES}/ringkasan.json", "w") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
