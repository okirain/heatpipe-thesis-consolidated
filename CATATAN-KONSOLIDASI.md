# Catatan Konsolidasi dan Audit Angka

Berkas ini bukan bagian dari skripsi. Isinya menjelaskan apa yang diambil dari
tujuh arsip sesi, apa yang diubah, dan mengapa angka pada skripsi konsolidasi
berbeda dari angka pada draf-draf sebelumnya.

## 1. Sumber yang Digabungkan

| Arsip sesi | Kata (draf) | Yang diambil | Yang tidak dipakai |
|:---|:---:|:---|:---|
| heatpipe-thesis-yuisa | 11.827 | Struktur Bab III paling lengkap (lima batas, korelasi sifat, jaringan termal), kerangka Bab V | Angka validasi dan margin (tidak dapat direproduksi) |
| heatpipe-thesis-jeremiah | 10.151 | Alur argumentasi Bab II, matriks perbandingan pompa vs pipa kalor | Tabel KRUSTY 1.200--4.950 W (itu daya total, bukan per pipa) |
| heatpipe-thesis-gamma | 11.104 | Pembahasan rezim batas per rentang temperatur, daftar lambang | Klaim MAPE 1,41% dan RMSE 0,042 K |
| heatpipe-thesis-beta | 9.623 | Gagasan dua tahap (lubang 7,65 mm gagal, lalu optimasi diameter) | Diameter 12,7 mm (gagal margin pada 11 kW), kalibrasi "MoS" |
| heatpipe-thesis-yukita | 10.190 | Pembahasan kecelakaan operasional dan skenario pelepasan pin | Profil 8 termokopel (tidak konsisten dengan R_th 0,005 K/W) |
| reduced-order-heatpipe | 3.214 | Daftar komponen neraca massa, skenario single-pipe failure | Klaim "choking pada 1,02 kW di 1000 K" (akibat korelasi P_sat salah) |
| final-project-s1 | -- | Templat UGM (header.tex, meta.tex, CSL, halaman depan) | Isi |

Total naskah konsolidasi: sekitar 11.300 kata prosa, 73 halaman A4, 13 tabel
bernomor, 9 gambar, dan 4 lampiran.

## 2. Pertentangan Angka Antar-Sesi dan Cara Menyelesaikannya

| Hal | Rentang klaim antar-sesi | Keputusan konsolidasi | Dasar |
|:---|:---|:---|:---|
| Diameter pipa | 12,7 / 15 / 15,88 / 18 mm | **18 mm** | Hasil sapuan: 16 mm memberi SF 1,02 (gagal), 18 mm memberi 1,37 |
| Temperatur kerja | 950 / 980 / 1050 / 1075 K | **1075 K (evaporator)** | SF ≥ 1,3 tercapai; T_UN 1249 K masih jauh di bawah batas |
| Batas pembatas | kapiler (semua sesi) | **entrainment** | Q_ent 15,22 kW < Q_sonic 29,41 kW < Q_cap 51,97 kW |
| MAPE validasi | 0,26--6,03% | **4,68% (leave-one-out)** | Protokol LOO; dua konstanta dilatih pada 4 titik, diuji pada titik ke-5 |
| R² validasi | 0,892--0,994 | **0,961** | Idem |
| Penghematan massa | 14,3--39,5 kg | **40,94 kg (63,0% subsistem)** | Neraca item per item; kisaran ketidakpastian 30--45 kg |
| Daya parasitik pompa | 2,5--18,5 kWe | **0,53 kWe** | Dihitung dari laju alir 2,08 kg/s, ΔP 30 kPa, efisiensi 15% |
| Bahan bakar | "tetap 155 kg" (yuisa, beta) | **-27,2% (112,8 kg)** | Aritmetika annulus: (32,65²-18,2²)/(32,65²-7,65²) |
| Umur teras | tidak dibahas | **17 -> 12,4 tahun** | Kesebandingan inventaris fisil pada daya tetap |

## 3. Kesalahan yang Diperbaiki

### 3.1 Korelasi tekanan uap jenuh

Draf-draf memakai

    ln P[MPa] = 11,9463 - 12633,73/T - 0,5 ln T

yang memberi 1,6 kPa pada 1000 K, sekitar 15 kali lebih rendah daripada nilai
termodinamika yang benar (24 kPa). Akibatnya densitas uap terlalu kecil, batas
sonik terlalu rendah, dan muncul kesimpulan keliru seperti "aliran tercekik pada
1 kW". Konsolidasi mengintegrasikan Clausius-Clapeyron dari titik didih normal
memakai h_fg(T) yang sama dengan model, sehingga P_sat dan h_fg konsisten.

### 3.2 Batas yang sebenarnya membatasi

Karena pipa harus masuk ke rongga pin, penampang uapnya kecil (14,6 mm) dan
kecepatan uapnya tinggi (133 m/s pada 11 kW). Pada geometri seperti itu, kriteria
Weber untuk entrainment tercapai lebih dulu daripada batas kapiler. Ini temuan
fisik, bukan sekadar angka: desain in-pin mengubah urutan batas yang lazim.

### 3.3 Daya pompa

Klaim 12,5--18,5 kWe terlalu besar satu orde. Daya hidraulik loop 500 kWth hanya
80 W; dengan efisiensi pompa EM 15% menjadi 534 W. Konsekuensinya, argumen
"penghematan listrik" melemah, dan argumen keselamatan menjadi alasan utama.

### 3.4 Massa dan bahan bakar

Beberapa draf mengklaim volume bahan bakar "tetap" sambil memperbesar rongga
tengah pin. Aritmetika annulus menunjukkan volume turun 27% pada tinggi aktif
tetap. Ini konsekuensi yang harus dinyatakan, bukan disembunyikan; skripsi
konsolidasi menyajikannya sebagai trade-off utama.

### 3.5 Ketidaksesuaian kode dan naskah

Pada arsip reduced-order-heatpipe, menjalankan `heatpipe_rom.py --mode
validate_krusty` menghasilkan ΔT_sim 6,69 / 11,29 / 13,60 K, sedangkan naskah
bab5 mencantumkan 11,29 / 14,79 / 16,60 K. Artinya sebagian angka pada naskah
tidak dihasilkan oleh kode yang disertakan. Karena itu seluruh angka skripsi
konsolidasi dihitung ulang dari satu model dan dapat direproduksi.

## 4. Yang Belum Diverifikasi (Penting)

1. **Data KRUSTY.** Lima kondisi tunak (312,5--1000 W) diambil dari arsip sesi,
   bukan diekstrak ulang dari NASA/TM-2018-219702. Sebelum skripsi diajukan,
   angka pada Tabel 5.2 dan Lampiran C.1 harus dicocokkan dengan tabel di
   laporan NASA tersebut. Jika berbeda, cukup ubah `KRUSTY_DATA` pada
   `src/run_all.py` lalu jalankan ulang; seluruh tabel dan gambar ikut berubah.
2. **Analisis neutronik.** Tidak dilakukan. Angka umur 12,4 tahun adalah
   perkiraan skala, bukan hasil perhitungan teras.
3. **Transien.** Model tunak; penyalaan, pembekuan, dan riwayat daya tidak
   tercakup.
4. **Kriteria entrainment.** Memakai kriteria Weber klasik; dasar empirisnya
   berasal dari percobaan air pada temperatur rendah, sehingga margin 1,37
   sebaiknya dibaca sebagai indikasi orde besaran.

## 5. Struktur Berkas

    heatpipe-thesis-consolidated/
      skripsi.qmd                  dokumen induk (include semua bab)
      _quarto.yml                  konfigurasi render PDF
      meta.tex                     judul, nama, NIM, pembimbing
      header.tex                   tipografi template DTNTF UGM
      01..06-*.qmd                 Bab I--VI
      front/                       halaman judul s.d. abstract
      lampiran-a..d.qmd            lampiran
      tables/                      tabel markdown hasil generate
      figures/                     gambar PDF hasil generate
      results/                     CSV dan ringkasan.json
      src/hp_model.py              model fisika
      src/run_all.py               seluruh studi + gambar
      src/make_tables.py           CSV -> tabel markdown
      _output/skripsi.pdf          hasil render

## 6. Cara Membangun Ulang

    # 1. hitung ulang seluruh angka (butuh numpy, pandas, matplotlib)
    cd src && python run_all.py && python make_tables.py

    # 2. render PDF (butuh Quarto + XeLaTeX)
    cd .. && quarto render skripsi.qmd

Model berjalan sekitar 60 detik untuk seluruh studi. Untuk mengubah geometri
atau temperatur desain, sunting konstanta di bagian atas `src/run_all.py`
(`T_DESIGN`, `T_COND`, `HP`) lalu ulangi langkah di atas.

## 7. Saran Reviewer

Jika naskah ini dibawa ke sidang, tiga hal berikut paling mungkin dipertanyakan
dan sudah dijawab di dalam naskah:

1. Mengapa MAPE 4,68% dan bukan di bawah 1% seperti draf lain? Jawaban ada di
   V.2: hampir seluruh ΔT terukur berasal dari antarmuka pengujian, dan angka
   rendah pada draf lain diperoleh dengan menyesuaikan model ke data yang sama
   yang dipakai untuk menguji.
2. Mengapa diameter 18 mm, bukan 12,7 mm yang kompatibel KRUSTY? Jawaban ada di
   V.5: 12,7 mm tidak mampu membawa 11 kW pada 1075 K (SF 0,48).
3. Mengapa entrainment, bukan kapiler? Jawaban ada di V.3: penampang uap kecil
   menaikkan kecepatan uap, sehingga kriteria Weber tercapai lebih dulu.
