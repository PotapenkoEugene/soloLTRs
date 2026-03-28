PYTHON := python3
LIBDIR := data/lib
T0DIR  := data/sim/t0
T1DIR  := data/sim/t1
T2DIR  := data/sim/t2

# ──────────────────────────────────────────────────
# Library (must run before any tier)
# ──────────────────────────────────────────────────
lib: $(LIBDIR)/ltr_families.fa

$(LIBDIR)/ltr_families.fa:
	$(PYTHON) scripts/sim/make_lib.py --outdir $(LIBDIR)

# ──────────────────────────────────────────────────
# Tier 0: unit tests (~10 KB, seconds)
# ──────────────────────────────────────────────────
t0: lib
	$(PYTHON) scripts/sim/make_t0.py --libdir $(LIBDIR) --outdir $(T0DIR)

# ──────────────────────────────────────────────────
# Tier 1: 1 Mb genome (~50 MB output)
# ──────────────────────────────────────────────────
t1: lib
	$(PYTHON) scripts/sim/make_t1.py --libdir $(LIBDIR) --outdir $(T1DIR)

t1-reads: t1
	@command -v art_illumina >/dev/null 2>&1 || \
	  (echo "art_illumina not found — run: nix shell nixpkgs#art" && exit 1)
	art_illumina -ss HS25 \
	  -i $(T1DIR)/genome/t1_genome.fa \
	  -l 150 -f 30 -m 500 -s 50 -p \
	  -o $(T1DIR)/reads/t1_reads
	@echo "Reads written to $(T1DIR)/reads/"

# ──────────────────────────────────────────────────
# Tier 2: 50 Mb genome, 5 families × 4 divergence bins
# ──────────────────────────────────────────────────
t2: lib
	$(PYTHON) scripts/sim/make_t2.py --libdir $(LIBDIR) --outdir $(T2DIR)

t2-reads-30x: t2
	@command -v art_illumina >/dev/null 2>&1 || \
	  (echo "art_illumina not found — run: nix shell nixpkgs#art" && exit 1)
	art_illumina -ss HS25 \
	  -i $(T2DIR)/genome/t2_genome.fa \
	  -l 150 -f 30 -m 500 -s 50 -p \
	  -o $(T2DIR)/reads_30x/t2_reads_30x

t2-reads-10x: t2
	@command -v art_illumina >/dev/null 2>&1 || \
	  (echo "art_illumina not found — run: nix shell nixpkgs#art" && exit 1)
	art_illumina -ss HS25 \
	  -i $(T2DIR)/genome/t2_genome.fa \
	  -l 150 -f 10 -m 500 -s 50 -p \
	  -o $(T2DIR)/reads_10x/t2_reads_10x

# ──────────────────────────────────────────────────
# Convenience targets
# ──────────────────────────────────────────────────
all-genomes: t0 t1 t2

all-reads: t1-reads t2-reads-30x t2-reads-10x

# ──────────────────────────────────────────────────
# Evaluate (run after pipeline produces predictions)
# ──────────────────────────────────────────────────
evaluate-t0-%:
	$(PYTHON) scripts/sim/evaluate.py \
	  --truth $(T0DIR)/$*/t0_$*.gff3 \
	  --pred $(PRED_GFF3) \
	  --tolerance $(or $(TOL),50)

evaluate-t1:
	$(PYTHON) scripts/sim/evaluate.py \
	  --truth $(T1DIR)/ground_truth/t1_ground_truth.gff3 \
	  --pred $(PRED_GFF3) \
	  --tolerance $(or $(TOL),50) \
	  --report results/eval_t1.tsv

evaluate-t2:
	$(PYTHON) scripts/sim/evaluate.py \
	  --truth $(T2DIR)/ground_truth/t2_ground_truth.gff3 \
	  --pred $(PRED_GFF3) \
	  --tolerance $(or $(TOL),50) \
	  --report results/eval_t2.tsv

# ──────────────────────────────────────────────────
# Sanity check: verify ground truth counts
# ──────────────────────────────────────────────────
check-t0:
	@echo "=== T0 clean ===" && grep -c "solo_LTR" $(T0DIR)/clean/t0_clean.gff3 || true
	@echo "=== T0 diverged ===" && grep -c "solo_LTR" $(T0DIR)/diverged/t0_diverged.gff3 || true
	@echo "=== T0 intact_nearby ===" && grep -c "solo_LTR\|LTR_retrotransposon" $(T0DIR)/intact_nearby/t0_intact_nearby.gff3 || true
	@echo "=== T0 edge_cases ===" && grep -c "solo_LTR" $(T0DIR)/edge_cases/t0_edge_cases.gff3 || true

check-t1:
	@echo "=== T1 solo LTRs ===" && grep -c "solo_LTR" $(T1DIR)/ground_truth/t1_ground_truth.gff3
	@echo "=== T1 intact ===" && grep -c "LTR_retrotransposon" $(T1DIR)/ground_truth/t1_ground_truth.gff3
	@echo "=== T1 genome size ===" && grep -v ">" $(T1DIR)/genome/t1_genome.fa | tr -d '\n' | wc -c

check-t2:
	@echo "=== T2 solo LTRs by bin ===" && grep "solo_LTR" $(T2DIR)/ground_truth/t2_ground_truth.gff3 | grep -oP "divergence_bin=[^;]+" | sort | uniq -c

.PHONY: lib t0 t1 t1-reads t2 t2-reads-30x t2-reads-10x all-genomes all-reads \
        check-t0 check-t1 check-t2 evaluate-t1 evaluate-t2
