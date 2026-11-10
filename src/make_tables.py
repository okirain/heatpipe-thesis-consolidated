"""make_tables.py -- Menghasilkan tabel markdown dari results/*.csv"""
import json, os
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
TAB = os.path.join(ROOT, "tables")
os.makedirs(TAB, exist_ok=True)


def n(x, nd=2, thousands=True):
    if x is None:
        return "-"
    s = f"{x:,.{nd}f}" if thousands else f"{x:.{nd}f}"
    return s.replace(",", "#").replace(".", ",").replace("#", ".")


def sci(x, nd=2):
    m, e = f"{x:.{nd}e}".split("e")
    return f"{m.replace('.', ',')}$\\times10^{{{int(e)}}}$"


def write(name, lines):
    with open(os.path.join(TAB, name), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


# --------------------------------------------------------------- sifat
df = pd.read_csv(f"{RES}/tabel_sifat.csv")
sel = df[df["T"].isin([800, 900, 1000, 1100, 1200])]
L = ["| $T$ (K) | $P_{sat}$ (kPa) | $\\rho_l$ (kg/m³) | $\\rho_v$ (kg/m³) | $\\mu_l$ (µPa·s) | $h_{fg}$ (MJ/kg) | $\\sigma$ (mN/m) | $k_l$ (W/m·K) |",
     "|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|"]
for _, r in sel.iterrows():
    L.append(f"| {n(r['T'],0)} | {n(r['P_sat_kPa'],2)} | {n(r['rho_l'],0)} | {n(r['rho_v'],4)} | "
             f"{n(r['mu_l'],1)} | {n(r['h_fg']/1e3,2)} | {n(r['sigma'],1)} | {n(r['k_l'],1)} |")
L.append(": Sifat termofisika natrium jenuh pada rentang operasi reaktor antariksa. {#tbl-sifat}")
write("tab_sifat.md", L)
write("tab_sifat_plain.md", [x for x in L if not x.startswith(":")])

# --------------------------------------------------------------- verifikasi
df = pd.read_csv(f"{RES}/verifikasi.csv")
L = ["| Uji | Parameter | Model | Acuan | Dev. (%) |", "|:---|:---|:---:|:---:|:---:|"]
nama = {"konvergensi_grid": "Grid", "neraca_energi": "Energi",
        "Mach_desain": "Mach", "tekanan_cairan": "Cairan",
        "tekanan_uap": "Uap", "jangkar_sifat": "Sifat"}
par_map = {
    "mdot*hfg,Q": r"$\dot m h_{fg}$ vs $Q$",
    "Ma_maks pada Q_puncak": r"$Ma_{maks}$",
    "solver vs Darcy analitik": "vs Darcy",
    "solver vs Hagen-Poiseuille": "vs Hagen-Poiseuille",
    "T_didih pada 101,3 kPa [K]": r"$T_{didih}$ (K)",
    "rho_l(1000 K) [kg/m3]": r"$\rho_l$ 1000 K (kg/m³)",
    "sigma(1000 K) [N/m]": r"$\sigma$ 1000 K (N/m)",
    "h_fg(1000 K) [MJ/kg]": r"$h_{fg}$ 1000 K (MJ/kg)",
    "P_sat(1000 K) [kPa]": r"$P_{sat}$ 1000 K (kPa)",
}
for _, r in df.iterrows():
    if r["uji"] == "Mach_desain":
        continue
    dev = "--" if pd.isna(r["galat_rel"]) else n(r["galat_rel"], 2)
    par = str(r["parameter"]).replace('"', "")
    par = par_map.get(par, par)
    nd = 2 if abs(r["referensi"]) > 100 else 4
    L.append(f"| {nama.get(r['uji'], r['uji'])} | {par} | {n(r['nilai'],nd)} | {n(r['referensi'],nd)} | {dev} |")
L.append(": Hasil verifikasi numerik dan konsistensi model. {#tbl-verifikasi}")
write("tab_verifikasi.md", L)
write("tab_verifikasi_plain.md", [x for x in L if not x.startswith(":")])

# --------------------------------------------------------------- KRUSTY
df = pd.read_csv(f"{RES}/krusty_validasi.csv")
L = ["| Titik uji | $Q$ (W) | $T_{e,exp}$ (K) | $T_{c,exp}$ (K) | $\\Delta T_{exp}$ (K) | $\\Delta T_{hp}$ (K) | $\\Delta T_{pred}$ (K) | Galat LOO (%) |",
     "|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|"]
for _, r in df.iterrows():
    L.append(f"| {r['run']} | {n(r['Q'],1)} | {n(r['T_evap_exp'],2)} | {n(r['T_cond_exp'],2)} | "
             f"{n(r['dT_exp'],1)} | {n(r['dT_hp'],2)} | {n(r['dT_loo'],2)} | {n(r['galat_rel'],2)} |")
L.append(": Perbandingan prediksi model dengan data uji KRUSTY; $\\Delta T_{hp}$ adalah beda suhu intrinsik pipa dan $\\Delta T_{pred}$ adalah prediksi *leave-one-out*. {#tbl-krusty}")
write("tab_krusty.md", L)

# --------------------------------------------------------------- amplop
df = pd.read_csv(f"{RES}/amplop.csv")
g0 = df[(df["kondisi"] == "0g") & (df["T"].isin([900, 950, 1000, 1050, 1075, 1100, 1150]))]
L = ["| $T$ (K) | $Q_{cap}$ (kW) | $Q_{sonic}$ (kW) | $Q_{ent}$ (kW) | $Q_{boil}$ (kW) | $Q_{visc}$ (kW) | Pembatas |",
     "|:---:|:---:|:---:|:---:|:---:|:---:|:---:|"]
lab = {"cap": "Kapiler", "sonic": "Sonik", "entrainment": "Entrainment",
       "boiling": "Pendidihan", "viscous": "Viskos"}
for _, r in g0.iterrows():
    L.append(f"| {n(r['T'],0)} | {n(r['Q_cap']/1e3,2)} | {n(r['Q_sonic']/1e3,2)} | {n(r['Q_ent']/1e3,2)} | "
             f"{n(r['Q_boil']/1e3,1)} | {n(r['Q_visc']/1e3,1)} | {lab[r['pembatas']]} |")
L.append(": Batas-batas perpindahan kalor pipa acuan (OD 18 mm) pada kondisi mikrogravitasi. {#tbl-amplop}")
write("tab_amplop.md", L)

# --------------------------------------------------------------- trade
df = pd.read_csv(f"{RES}/trade_geometri.csv")
L = ["| OD (mm) | $D_v$ (mm) | $Q_{cap}$ (kW) | $Q_{sonic}$ (kW) | $Q_{ent}$ (kW) | SF (-) | Kehilangan UN (%) |",
     "|:---:|:---:|:---:|:---:|:---:|:---:|:---:|"]
for _, r in df.iterrows():
    mark = "**" if abs(r["OD_mm"] - 18) < 0.01 else ""
    L.append(f"| {mark}{n(r['OD_mm'],0)}{mark} | {n(r['D_v_mm'],1)} | {n(r['Q_cap'],2)} | {n(r['Q_sonic'],2)} | "
             f"{n(r['Q_ent'],2)} | {n(r['SF'],2)} | {n(r['fuel_loss'],1)} |")
L.append(": Pengaruh diameter luar pipa terhadap batas operasi pada 1075 K dan konsekuensinya pada volume bahan bakar. {#tbl-trade}")
write("tab_trade.md", L)

# --------------------------------------------------------------- teras
df = pd.read_csv(f"{RES}/teras.csv")
L = ["| Cincin | $N$ | $F_r$ | $Q_{pin}$ (kW) | $Q_{cin}$ (kW) | $T_{clad}$ (K) | $T_{UN}$ (K) | $\\Delta P_{req}$ (Pa) | $Ma$ | Margin (%) | SF |",
     "|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|"]
for _, r in df.iterrows():
    nm = str(r['cincin']).replace(" (pusat)", "").replace(" (tepi)", "").replace("Cincin ", "C")
    L.append(f"| {nm} | {int(r['n'])} | {n(r['fr'],2)} | {n(r['Q']/1e3,2)} | {n(r['Q_tot']/1e3,2)} | "
             f"{n(r['T_clad'],1)} | {n(r['T_fuel'],1)} | {n(r['dP_req'],0)} | {n(r['Ma'],3)} | "
             f"{n(r['margin_cap'],1)} | {n(r['SF'],2)} |")
tot = df["Q_tot"].sum()
L.append(f"| **Total** | **61** | **1,00** | **{n(tot/61/1e3,2)}** | **{n(tot/1e3,2)}** | "
         f"**{n(df['T_clad'].max(),1)}** | **{n(df['T_fuel'].max(),1)}** | **{n(df['dP_req'].max(),0)}** | "
         f"**{n(df['Ma'].max(),3)}** | **{n(df['margin_cap'].min(),1)}** | **{n(df['SF'].min(),2)}** |")
L.append(": Distribusi daya dan keadaan termohidrolika teras ISR 61 pin pada 500 kWth, $T_v$ = 1074 K, $g$ = 0; C1 adalah pin pusat dan C5 adalah cincin terluar. {#tbl-teras}")
write("tab_teras.md", L)

# --------------------------------------------------------------- massa
df = pd.read_csv(f"{RES}/neraca_massa.csv")
base = df[(df.kategori.str.contains("Baseline")) & (df.item != "SUBTOTAL")]
retr = df[(df.kategori.str.contains("Retrofit")) & (df.item != "SUBTOTAL")]
sb = base["massa"].sum(); sr = retr["massa"].sum()
L = ["| Komponen loop pompa Na-78K | Massa (kg) |", "|:---|:---:|"]
for _, r in base.iterrows():
    L.append(f"| {r['item']} | {n(r['massa'],2)} |")
L.append(f"| **Subtotal baseline** | **{n(sb,2)}** |")
L.append(": Estimasi massa subsistem pendingin teras pada desain baseline berpompa Na-78K. {#tbl-massa-base}")
write("tab_massa_base.md", L)

L = ["| Komponen sistem pipa kalor | Massa (kg) |", "|:---|:---:|"]
for _, r in retr.iterrows():
    L.append(f"| {r['item']} | {n(r['massa'],2)} |")
L.append(f"| **Subtotal retrofit** | **{n(sr,2)}** |")
L.append(": Estimasi massa subsistem pendingin teras pada desain retrofit 61 pipa kalor natrium. {#tbl-massa-retr}")
write("tab_massa_retr.md", L)

L = ["| Besaran | Nilai |", "|:---|:---:|",
     f"| Massa baseline | {n(sb,2)} kg |",
     f"| Massa retrofit | {n(sr,2)} kg |",
     f"| **Selisih bersih** | **$-$**{n(sb-sr,2)} kg ({n((sb-sr)/sb*100,1)}%) |"]
L.append(": Ringkasan selisih massa subsistem pendingin teras. {#tbl-massa}")
write("tab_massa.md", L)

# --------------------------------------------------------------- pompa
df = pd.read_csv(f"{RES}/daya_pompa.csv")
L = ["| Besaran | Nilai | Satuan |", "|:---|:---:|:---:|"]
for _, r in df.iterrows():
    L.append(f"| {r['besaran']} | {n(r['nilai'],3)} | {r['satuan']} |")
L.append(": Estimasi daya pompa elektromagnetik loop Na-78K pada daya teras 500 kWth. {#tbl-pompa}")
write("tab_pompa.md", L)

# --------------------------------------------------------------- bahan bakar
df = pd.read_csv(f"{RES}/bahan_bakar.csv")
L = ["| Besaran | Nilai | Satuan |", "|:---|:---:|:---:|"]
for _, r in df.iterrows():
    sat = str(r['satuan']).replace("cm3", "cm³")
    L.append(f"| {r['besaran']} | {n(r['nilai'],2)} | {sat} |")
L.append(": Konsekuensi geometri pipa 18 mm terhadap inventaris bahan bakar dan umur teras. {#tbl-bahan-bakar}")
write("tab_bahan_bakar.md", L)

# --------------------------------------------------------------- sensitivitas
df = pd.read_csv(f"{RES}/sensitivitas.csv")
L = ["| $T_v$ (K) | $Q_{env}$ (kW) | SF (-) | $Ma_{maks}$ (-) | $\\Delta P_{req}$ (Pa) | $T_{UN,maks}$ (K) |",
     "|:---:|:---:|:---:|:---:|:---:|:---:|"]
for _, r in df.iterrows():
    L.append(f"| {n(r['T'],0)} | {n(r['Q_env']/1e3,2)} | {n(r['SF'],2)} | {n(r['Ma'],3)} | "
             f"{n(r['dP_req'],0)} | {n(r['T_fuel'],1)} |")
L.append(": Sensitivitas kinerja pipa dan suhu bahan bakar terhadap temperatur operasi. {#tbl-sens}")
write("tab_sens.md", L)
write("tab_sens_plain.md", [x for x in L if not x.startswith(":")])

print("tabel selesai")