import os
import json
import numpy as np
import pandas as pd
from scipy.signal import lfilter
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import filedialog

# -----------------------
# Dataset selection
# -----------------------
root = tk.Tk()
root.withdraw()
dataset_root = filedialog.askdirectory(title="Select dataset root folder")
if not dataset_root:
    raise RuntimeError("No dataset folder selected!")
print(f"Using dataset root: {dataset_root}")

# -----------------------
# Utility functions
# -----------------------
def load_sample(path):
    rx = np.load(os.path.join(path, "rx.npy"))
    with open(os.path.join(path, "meta.json"), "r") as f:
        meta = json.load(f)
    return rx, meta

def matched_filter_and_timing(rx, sps=8):
    mf = np.ones(sps)
    mf_out = lfilter(mf, 1, rx)
    mu = 0.0
    out, i = [], 0
    while i < len(mf_out) - sps:
        out.append(mf_out[i + int(mu)])
        mu += 1.0
        if mu >= sps:
            mu -= sps
        i += sps
    return np.array(out)

def bpsk_llr(symbols, noise_var=1.0):
    return 2 * np.real(symbols) / (noise_var + 1e-12)

def bpsk_demod(symbols):
    return (np.real(symbols) < 0).astype(int)

# -----------------------
# Phase 3: Convolutional (soft Viterbi)
# -----------------------
def viterbi_decode_soft(llrs):
    G = [0o133, 0o171]
    K = 7
    n_states = 2 ** (K - 1)
    n_bits = len(llrs) // 2

    next_state = np.zeros((n_states, 2), dtype=int)
    outputs = np.zeros((n_states, 2, 2), dtype=int)
    for state in range(n_states):
        for inp in [0, 1]:
            shift_reg = (state << 1) | inp
            out = [bin(shift_reg & g).count("1") % 2 for g in G]
            outputs[state, inp] = out
            next_state[state, inp] = shift_reg & (n_states - 1)

    path_metrics = np.full(n_states, np.inf)
    path_metrics[0] = 0
    paths = np.full((n_states, n_bits), -1, dtype=int)

    for t in range(n_bits):
        new_metrics = np.full(n_states, np.inf)
        new_paths = np.full((n_states, n_bits), -1, dtype=int)
        rx_llr = llrs[2 * t:2 * t + 2]
        for state in range(n_states):
            if path_metrics[state] < np.inf:
                for inp in [0, 1]:
                    ns = next_state[state, inp]
                    exp_out = outputs[state, inp]
                    branch_metric = -np.sum(rx_llr * (1 - 2 * np.array(exp_out)))
                    metric = path_metrics[state] + branch_metric
                    if metric < new_metrics[ns]:
                        new_metrics[ns] = metric
                        new_paths[ns] = paths[state]
                        new_paths[ns, t] = inp
        path_metrics = new_metrics
        paths = new_paths

    best_state = np.argmin(path_metrics)
    decoded = paths[best_state]
    return decoded[decoded >= 0]

# -----------------------
# Phase 3: Reed–Solomon (15,11)
# -----------------------
PRIM_POLY = 0x13
GF_M = 4
GF_Q = 1 << GF_M
GF_N = 15
GF_K = 11
RS_T = 2

GF_EXP = np.zeros(2*GF_N, dtype=np.int16)
GF_LOG = np.full(GF_Q, -1, dtype=np.int16)

def _gf_init():
    x = 1
    for i in range(GF_N):
        GF_EXP[i] = x
        if GF_LOG[x] == -1:
            GF_LOG[x] = i
        x <<= 1
        if x & GF_Q:
            x ^= PRIM_POLY
    for i in range(GF_N, 2*GF_N):
        GF_EXP[i] = GF_EXP[i - GF_N]
_gf_init()

def gf_mul(a, b):
    if a == 0 or b == 0:
        return 0
    return int(GF_EXP[(GF_LOG[a] + GF_LOG[b]) % GF_N])

def gf_div(a, b):
    if b == 0:
        raise ZeroDivisionError
    if a == 0:
        return 0
    return int(GF_EXP[(GF_LOG[a] - GF_LOG[b]) % GF_N])

def gf_pow(a, n):
    if a == 0:
        return 0
    return int(GF_EXP[(GF_LOG[a] * (n % GF_N)) % GF_N])

def poly_eval(poly, x):
    y = 0
    for c in poly:
        y = gf_mul(y, x) ^ c
    return y

def syndromes(cw, nsyn):
    return [poly_eval(cw, GF_EXP[i % GF_N]) for i in range(1, nsyn+1)]

def berlekamp_massey(S):
    C = [1] + [0]*len(S)
    B = [1] + [0]*len(S)
    L = 0; m = 1; b = 1
    for n in range(len(S)):
        d = S[n]
        for i in range(1, L+1):
            d ^= gf_mul(C[i], S[n - i])
        if d != 0:
            T = C.copy()
            coef = gf_div(d, b)
            for i in range(m, len(C)):
                C[i] ^= gf_mul(coef, B[i - m])
            if 2*L <= n:
                L = n + 1 - L
                B = T
                b = d
                m = 1
            else:
                m += 1
        else:
            m += 1
    while len(C) > 0 and C[-1] == 0:
        C.pop()
    return C, L

def chien_search(locator):
    roots = []
    for i in range(GF_N):
        xinvi = GF_EXP[(GF_N - i) % GF_N]
        if poly_eval(locator, xinvi) == 0:
            roots.append(i)
    return roots

def forney_magnitudes(S, locator, err_pos_pows):
    Omega = [0]*(len(locator) - 1 + len(S))
    for i, si in enumerate(S):
        for j, lj in enumerate(locator):
            Omega[i + j] ^= gf_mul(si, lj)
    Omega = Omega[:len(S)]
    deriv = []
    for i in range(1, len(locator), 2):
        deriv.append(locator[i])
    mags = []
    for pos in err_pos_pows:
        xinvi = GF_EXP[(GF_N - pos) % GF_N]
        num = 0
        for k, c in enumerate(Omega):
            num ^= gf_mul(c, gf_pow(xinvi, k))
        den = 0
        for k, c in enumerate(deriv):
            den ^= gf_mul(c, gf_pow(xinvi, k))
        mags.append(0 if den == 0 else gf_div(num, den))
    return mags

def bits_to_nibbles_msb(bits):
    bits = np.asarray(bits, dtype=np.uint8)
    n = len(bits) // 4
    out = []
    for i in range(n):
        b3,b2,b1,b0 = bits[4*i:4*i+4]
        out.append(int((b3<<3)|(b2<<2)|(b1<<1)|b0))
    return out

def nibbles_to_bits_msb(nibs):
    out = []
    for v in nibs:
        out += [(v>>3)&1, (v>>2)&1, (v>>1)&1, v&1]
    return np.array(out, dtype=np.int8)

def _rs1511_decode_blocks(sym_list):
    data_out = []
    valid_blocks = 0
    total_blocks = 0
    for i in range(0, len(sym_list) - (len(sym_list) % GF_N), GF_N):
        total_blocks += 1
        cw = sym_list[i:i+GF_N]
        S = syndromes(cw, 2*RS_T)
        if max(S) == 0:
            data_out.extend(cw[:GF_K])
            valid_blocks += 1
            continue
        locator, L = berlekamp_massey(S)
        if L == 0 or L > RS_T:
            data_out.extend(cw[:GF_K]); continue
        err_pos = chien_search(locator)
        if len(err_pos) < L:
            data_out.extend(cw[:GF_K]); continue
        mags = forney_magnitudes(S, locator, err_pos)
        cw_corr = cw[:]
        for idx, mag in zip([GF_N-1-p for p in err_pos], mags):
            if 0 <= idx < GF_N:
                cw_corr[idx] ^= mag
        if max(syndromes(cw_corr, 2*RS_T)) == 0:
            data_out.extend(cw_corr[:GF_K])
            valid_blocks += 1
        else:
            data_out.extend(cw[:GF_K])
    return data_out, valid_blocks, total_blocks

def rs1511_decode(bits):
    bits = np.asarray(bits, dtype=np.uint8)
    n_nibbles = len(bits) // 4
    if n_nibbles < GF_N:
        return np.zeros(0, dtype=np.int8)
    syms = bits_to_nibbles_msb(bits[:n_nibbles*4])
    data, valid_blocks, total_blocks = _rs1511_decode_blocks(syms)
    decoded_bits = nibbles_to_bits_msb(data)
    fer = 1 - (valid_blocks / total_blocks if total_blocks > 0 else 0)
    return decoded_bits, fer

# -----------------------
# Phase 4: Doppler CFO + Costas loop
# -----------------------
def estimate_cfo_kay(x, fs, K=8):
    N = len(x)
    K = min(K, N - 2)
    C = np.array([np.vdot(x[k:], x[:-k]) for k in range(1, K + 1)])
    phi = np.unwrap(np.angle(C))
    slope = np.polyfit(np.arange(1, K + 1), phi, 1)[0]
    return (fs / (2 * np.pi)) * slope

def apply_cfo_correction(x, fs, cfo_hz):
    n = np.arange(len(x))
    rot = np.exp(-1j * 2 * np.pi * cfo_hz * n / fs)
    return x * rot

def costas_loop(x, beta=0.01, alpha=0.001):
    out = np.zeros_like(x, dtype=complex)
    phase = 0.0; freq = 0.0
    for n, s in enumerate(x):
        rot = np.exp(-1j * phase)
        y = s * rot
        out[n] = y
        err = np.sign(np.real(y)) * np.imag(y)
        freq += beta * err
        phase += freq + alpha * err
    return out

# -----------------------
# Main processing
# -----------------------
def process_folder(folder, do_timing, do_cfo, coding=None):
    results = []
    for root, _, files in os.walk(folder):
        if "rx.npy" not in files or "meta.json" not in files:
            continue
        rx, meta = load_sample(root)
        clean_bits = np.array(meta.get("clean_bits", []))
        sps = int(meta.get("sps", 8))
        fs = int(meta.get("fs", 48000))
        if do_timing:
            rx = matched_filter_and_timing(rx, sps=sps)
        else:
            rx = rx[::sps]
        if do_cfo:
            cfo_hat = estimate_cfo_kay(rx, fs, K=8)
            rx = apply_cfo_correction(rx, fs, cfo_hat)
            rx = costas_loop(rx)

        fer = None
        if coding == "conv":
            llrs = bpsk_llr(rx)
            decoded = viterbi_decode_soft(llrs)
        elif coding == "rs":
            hard_bits = bpsk_demod(rx)
            decoded, fer = rs1511_decode(hard_bits)
        else:
            hard_bits = bpsk_demod(rx)
            decoded = hard_bits

        out_path = os.path.join(root, "decoded_bits.npy")
        np.save(out_path, decoded)

        min_len = min(len(clean_bits), len(decoded))
        ber = np.mean(clean_bits[:min_len] != decoded[:min_len]) if min_len > 0 else 1.0

        results.append({
            "dataset": os.path.basename(os.path.dirname(folder)),
            "folder": os.path.basename(folder),
            "sample": os.path.basename(root),
            "coding": coding or "none",
            "ber": ber,
            "fer": fer,  # <-- NEW COLUMN
            "decoded_len": len(decoded),
        })
    return results

# -----------------------
# Dataset roots
# -----------------------
roots = [
    (os.path.join(dataset_root, "phase2_snr"), False, False, None),
    (os.path.join(dataset_root, "phase1_timing"), True, False, None),
    (os.path.join(dataset_root, "phase4_doppler"), True, True, None),
    (os.path.join(dataset_root, "phase3_coding", "convolutional"), True, False, "conv"),
    (os.path.join(dataset_root, "phase3_coding", "reed_solomon"), True, False, "rs"),
]

# -----------------------
# Run all phases
# -----------------------
all_results = []
for folder, timing, cfo, coding in roots:
    if os.path.exists(folder):
        print(f"Processing {folder} ...")
        all_results.extend(process_folder(folder, timing, cfo, coding))

csv_path = os.path.join(dataset_root, "bpsk_full_pipeline_results.csv")
pd.DataFrame(all_results).to_csv(csv_path, index=False)
print(f"Results saved to {csv_path}")

# -----------------------
# Save paths (Desktop/Plots)
# -----------------------
desktop = os.path.join(os.path.expanduser("~"), "Desktop")
plot_dir = os.path.join(desktop, "Plots")
const_dir = os.path.join(plot_dir, "Constellations")
doppler_dir = os.path.join(plot_dir, "Doppler")
os.makedirs(const_dir, exist_ok=True)
os.makedirs(doppler_dir, exist_ok=True)

# -----------------------
# Constellation Plots (Phase 2)
# -----------------------
base_dir = os.path.join(dataset_root, "phase2_snr")
snr_folders = ["snr_0db", "snr_5db", "snr_10db", "snr_15db"]

for snr_folder in snr_folders:
    folder_path = os.path.join(base_dir, snr_folder, "sample_000")
    rx_file = os.path.join(folder_path, "rx.npy")
    if os.path.exists(rx_file):
        rx = np.load(rx_file)
        plt.figure(figsize=(6, 6))
        plt.scatter(rx.real, rx.imag, s=1, alpha=0.5, label=f"{snr_folder}")
        plt.scatter([-1, 1], [0, 0], c="red", marker="x", s=200, label="Ideal BPSK")
        plt.axhline(0, color="gray", linestyle="--", linewidth=0.5)
        plt.axvline(0, color="gray", linestyle="--", linewidth=0.5)
        plt.xlabel("In-phase (I)")
        plt.ylabel("Quadrature (Q)")
        plt.title(f"Constellation Diagram at {snr_folder}")
        plt.legend()
        plt.grid(True, linestyle="--", alpha=0.5)
        save_path = os.path.join(const_dir, f"constellation_{snr_folder}.png")
        plt.savefig(save_path, dpi=300)
        plt.close()

# -----------------------
# Doppler Compensation Plots (Phase 4)
# -----------------------
base_path = os.path.join(dataset_root, "phase4_doppler")

def doppler_correct(rx, freq_offset, fs):
    n = np.arange(len(rx))
    return rx * np.exp(-1j * 2 * np.pi * freq_offset * n / fs)

if os.path.exists(base_path):
    snr_folders = [f for f in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, f))]
    results = {}
    for snr_folder in snr_folders:
        snr_path = os.path.join(base_path, snr_folder)
        sample_folders = [f for f in os.listdir(snr_path) if os.path.isdir(os.path.join(snr_path, f))]
        for sample_folder in sample_folders:
            sample_path = os.path.join(snr_path, sample_folder)
            rx_path = os.path.join(sample_path, "rx.npy")
            if not os.path.exists(rx_path):
                continue
            rx = np.load(rx_path)
            fs = 1.0
            freq_offset = 0.01
            rx_uncorrected = rx
            rx_corrected = doppler_correct(rx, freq_offset, fs)
            if snr_folder not in results:
                results[snr_folder] = []
            results[snr_folder].append((rx_uncorrected, rx_corrected))
    for snr_folder, samples in results.items():
        plt.figure(figsize=(10, 5))
        rx_uncorrected, rx_corrected = samples[0]
        plt.subplot(1, 2, 1)
        plt.scatter(rx_uncorrected.real[:200], rx_uncorrected.imag[:200], s=5, alpha=0.5)
        plt.title(f"{snr_folder} - Before Doppler Correction")
        plt.xlabel("In-Phase"); plt.ylabel("Quadrature")
        plt.subplot(1, 2, 2)
        plt.scatter(rx_corrected.real[:200], rx_corrected.imag[:200], s=5, alpha=0.5)
        plt.title(f"{snr_folder} - After Doppler Correction")
        plt.xlabel("In-Phase"); plt.ylabel("Quadrature")
        plt.tight_layout()
        save_path = os.path.join(doppler_dir, f"doppler_{snr_folder}.png")
        plt.savefig(save_path, dpi=300)
        plt.close()
