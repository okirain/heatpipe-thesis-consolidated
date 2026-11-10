
src = open("hp_model.py", encoding="utf-8").read().splitlines()

def block(start_pat, end_pat, label):
    i = next(k for k, l in enumerate(src) if start_pat in l)
    j = next(k for k in range(i + 1, len(src)) if end_pat in src[k])
    return label, src[i:j]

parts = [
    block("def h_fg(", "def rho_l", "# --- Sifat natrium: panas laten, tekanan uap jenuh, temperatur saturasi ---"),
    block("@dataclass", "def mdot_profile", "# --- Geometri pipa kalor dan properti turunannya ---"),
    block("def solve_axial(", "def dP_cap_max", "# --- Penyelesaian aksial uap dan cairan (integrasi mundur Runge-Kutta) ---"),
    block("def q_capillary(", "def envelope", "# --- Batas kapiler dengan metode bagi dua ---"),
    block("def annular_fuel_dT(", "def krusty_pipe", "# --- Konduksi radial pelet bahan bakar anular ---"),
]
out = ['"""Kutipan inti hp_model.py; berkas lengkap tersedia pada folder src."""', ""]
for label, lines in parts:
    out.append(label)
    body = list(lines)
    while body and not body[-1].strip():
        body.pop()
    out.extend(body)
    out += ["", "# ...", ""]
open("lampiran_kode.py", "w", encoding="utf-8").write("\n".join(out))
print(len(out), "baris")
