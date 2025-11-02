import os
import matplotlib.pyplot as plt

def plot_binding_map(sequence, primers, results_dir):
    seq_len = len(sequence)
    plt.figure(figsize=(10,2))
    y = 1
    pairs = primers['PRIMER_PAIR_NUM_RETURNED']

    for i in range(pairs):
        start_l = primers[f'PRIMER_LEFT_{i}'][0]
        end_r   = primers[f'PRIMER_RIGHT_{i}'][0]
        plt.hlines(y, start_l, end_r, colors='royalblue', lw=4)
        plt.text(start_l, y+0.05, f"P{i+1}", size=8)
        y += 0.25

    plt.xlabel("Position (bp)")
    plt.title("Primer Binding Map")
    plt.xlim(0, seq_len)
    plt.ylim(0.8, 1.8)
    plt.tight_layout()
    out_png = os.path.join(results_dir, "primer_binding_map.png")
    plt.savefig(out_png, dpi=220)
    plt.close()
    print(f"🧬 Binding map saved to {out_png}")
