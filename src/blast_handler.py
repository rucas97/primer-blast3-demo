import subprocess
from Bio import SearchIO

def run_blast(primer_seq, ref_fasta, results_dir):
    """Run local BLAST+ with XML output and parse top hits."""
    primer_fasta = f"{results_dir}/temp_primer.fasta"
    out_xml      = f"{results_dir}/temp_blast.xml"

    with open(primer_fasta, "w") as f:
        f.write(">primer\n" + primer_seq + "\n")

    cmd = [
        "blastn", "-query", primer_fasta, "-db", ref_fasta,
        "-task", "blastn-short", "-word_size", "7",
        "-evalue", "1000", "-outfmt", "5",
        "-max_target_seqs", "5", "-out", out_xml
    ]
    subprocess.run(cmd, check=True)

    hits = []
    for q in SearchIO.parse(out_xml, "blast-xml"):
        for h in q.hits[:3]:
            hs = h.hsps[0]
            alen = getattr(hs, "alignment_length", getattr(hs, "align_length", 0))
            ident = (hs.ident_num / alen * 100) if alen else 0
            hits.append(f"{h.id} | identity={ident:.1f}% | len={alen}")
    return hits or ["No hit found"]
