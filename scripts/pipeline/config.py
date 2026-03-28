"""Default parameters for the soloLTRs pipeline."""

# --- Alignment ---
DEFAULT_THREADS = 8
MIN_MAPQ_ALIGN = 1       # keep multi-mapped within a family, discard MAPQ=0
MIN_MAPQ_DEPTH = 20      # filter for depth calculation; keeps internal region clean
                         # (MAPQ 12-19 reads are mostly false-positive in internal regions)

# BWA-MEM2 scoring for sensitive mode (optimised for diverged TE copies)
# Default bwa-mem2: -B 4 -T 30 -k 19  → misses copies >15% divergence
# Sensitive mode:   -B 2 -T 15 -k 13  → recovers copies up to ~25% divergence
#   -B 2: lower mismatch penalty (150bp read at 25% div: 105 match - 45×2 = 15 → passes -T 15)
#   -T 15: lower minimum alignment score threshold
#   -k 13: shorter minimum seed length (finds seeds in more diverged reads)
BWA_MISMATCH_PENALTY = 2   # -B
BWA_MIN_SCORE = 15          # -T
BWA_MIN_SEED_LEN = 13       # -k

# --- Depth calculation ---
EDGE_TRIM_BP = 50        # exclude this many bp at each region boundary
MIN_LTR_DEPTH = 10.0     # minimum mean LTR depth to call (below = low_coverage)
MIN_MAPPED_READS = 50    # minimum total mapped reads per family for confidence

# --- Ratio thresholds ---
R_MAX = 0.50             # R > 0.50 is theoretically impossible; clamp + report all_intact
R_MIN = 0.01             # R < 0.01 → near_all_solo (S/C > 48)
R_SOLO_THRESHOLD = 0.45  # binary PAV: solo_detected if R < this value
                         # R=0.45 ↔ S/C ≈ 0.22 (any meaningful solo LTR presence)

# --- Bootstrap ---
N_BOOTSTRAP = 1000
CI_ALPHA = 0.05          # 95% confidence intervals
