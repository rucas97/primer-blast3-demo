import os
import subprocess
from Bio import SeqIO, SearchIO
import primer3
import matplotlib.pyplot as plt
import seaborn as sns

# --- Configuration ---
INPUT_FASTA = 'example_input/target_sequence.fasta'
REFERENCE_FASTA = 'example_input/genome_reference.fasta'
RESULTS_DIR = 'results'
os.makedirs(RESULTS_DIR, exist_ok=True)


def read_fasta(path):
    """Read the first record from a FASTA file and return sequence string."""
    record = next(SeqIO.parse(path, "fasta"))
    return str(record.seq).upper().replace("\n", "").replace(" ", "")


def design_primers(sequence):
    """Design PCR primers using Primer3."""
    primers = primer3.bindings.design_primers(
        {
            "SEQUENCE_ID": "demo_target",
            "SEQUENCE_TEMPLATE": sequence,
            "SEQUENCE_INCLUDED_REGION": [0, len(sequence)],
        },
        {
            "PRIMER_OPT_SIZE": 20,
            "PRIMER_MIN_SIZE": 18,
            "PRIMER_MAX_SIZE": 27,
            "PRIMER_OPT_TM": 60.0,
            "PRIMER_MIN_TM": 57.0,
            "PRIMER_MAX_TM": 63.0,
            "PRIMER_MIN_GC": 35.0,
            "PRIMER_MAX_GC": 65.0,
            "PRIMER_PRODUCT_SIZE_RANGE": [[30, 300]],
            "PRIMER_NUM_RETURN": 5,
        },
    )
    return primers


def run_blast(primer_seq):
    """
    Run local BLAST+ to check primer specificity.
    Requires genome_reference.fasta indexed with makeblastdb.
    """
    primer_fasta = os.path.join(RESULTS_DIR, "temp_primer.fasta")
    out_xml = os.path.join(RESULTS_DIR, "temp_blast.xml")

    # Write primer FASTA
    with open(primer_fasta, "w") as f:
        f.write(">primer\n")
        f.write(primer_seq + "\n")

    # Build command (fully subprocess-style, no Bio.Application)
    blast_cmd = [
        "blastn",
        "-query", primer_fasta,
        "-db", REFERENCE_FASTA,
        "-task", "blastn-short",
        "-word_size", "7",
        "-evalue", "1000",
        "-reward", "2",
        "-penalty", "-3",
        "-outfmt", "5",
        "-max_target_seqs", "5",
        "-out", out_xml,
    ]

    try:
        subprocess.run(blast_cmd, check=True, capture_output=True, text=True)
        hits = []

        for qresult in SearchIO.parse(out_xml, "blast-xml"):
            for hit in qresult.hits[:3]:
                hsp = hit.hsps[0]
                alen = getattr(hsp, "alignment_length", None) or getattr(hsp, "align_length", None) or 0
                identity = (hsp.ident_num / alen * 100) if alen else 0
                start = hsp.hit_start + 1
                end = hsp.hit_end
                hits.append(
                    f"{hit.id} | identity={identity:.1f}% | range={start}-{end} | length={alen}"
                )

        if not hits:
            hits.append("No local hit detected")

        return hits
    except Exception as e:
        return [f"Local BLAST error: {e}"]


def write_results(res):
    """Write primer design details and BLAST results into text files."""
    blast_report = os.path.join(RESULTS_DIR, "blast_report.txt")
    out_path = os.path.join(RESULTS_DIR, "output_primers.txt")

    pairs = res.get("PRIMER_PAIR_NUM_RETURNED", 0)
    with open(out_path, "w") as f_primers, open(blast_report, "w") as f_blast:
        for i in range(pairs):
            left = res[f"PRIMER_LEFT_{i}_SEQUENCE"]
            right = res[f"PRIMER_RIGHT_{i}_SEQUENCE"]
            tm_left = res[f"PRIMER_LEFT_{i}_TM"]
            tm_right = res[f"PRIMER_RIGHT_{i}_TM"]
            gc_left = res[f"PRIMER_LEFT_{i}_GC_PERCENT"]
            gc_right = res[f"PRIMER_RIGHT_{i}_GC_PERCENT"]
            size = res[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"]

            # Primer info file
            f_primers.write(f"Primer pair {i + 1}\n")
            f_primers.write(f" Forward: {left} | Tm={tm_left:.2f} | GC%={gc_left:.2f}\n")
            f_primers.write(f" Reverse: {right} | Tm={tm_right:.2f} | GC%={gc_right:.2f}\n")
            f_primers.write(f" Product size: {size} bp\n\n")

            # Local BLAST analysis
            hits_f = run_blast(left)
            hits_r = run_blast(right)
            f_blast.write(f"Primer pair {i + 1}\n")
            f_blast.write(" Forward BLAST top hits:\n")
            for h in hits_f:
                f_blast.write(f"   {h}\n")
            f_blast.write(" Reverse BLAST top hits:\n")
            for h in hits_r:
                f_blast.write(f"   {h}\n")
            f_blast.write("-" * 50 + "\n\n")

    print(f"✅ Primer design + BLAST results saved:\n{out_path}\n{blast_report}")


def plot_binding_map(sequence, res):
    """Visualize primer binding sites on the target sequence."""
    plt.figure(figsize=(10, 2))
    seq_len = len(sequence)
    y = 1.0
    pairs = res.get("PRIMER_PAIR_NUM_RETURNED", 0)

    for i in range(pairs):
        start_left = res[f"PRIMER_LEFT_{i}"][0]
        end_right = res[f"PRIMER_RIGHT_{i}"][0]
        plt.hlines(y, start_left, end_right, color="royalblue", lw=4)
        plt.text(start_left, y + 0.05, f"P{i+1}", fontsize=8)
        y += 0.25

    plt.title("Primer Binding Map")
    plt.xlabel("Nucleotide position (bp)")
    plt.xlim(0, seq_len)
    plt.ylim(0.8, 1.8)
    plt.tight_layout()
    out_png = os.path.join(RESULTS_DIR, "primer_binding_map.png")
    plt.savefig(out_png, dpi=220)
    plt.close()
    print(f"🧬 Binding map saved to: {out_png}")


if __name__ == "__main__":
    seq = read_fasta(INPUT_FASTA)
    primers = design_primers(seq)
    write_results(primers)
    plot_binding_map(seq, primers)
