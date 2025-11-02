def score_primer(f_tm, r_tm, prod_size, f_gc, r_gc):
    """Composite score 0‑1000 for ranking primer pairs."""
    tm_delta = abs(f_tm - r_tm)
    gc_mid   = abs((f_gc + r_gc)/2 - 50)
    size_pen = abs(prod_size - 150) / 3
    score    = 1000 - (tm_delta*10 + gc_mid*5 + size_pen)
    return round(max(score, 0), 1)
