PYTHON   := python3
LIBDIR   := data/lib
T0DIR    := data/sim/t0
T1DIR    := data/sim/t1
T2DIR    := data/sim/t2
PIPLIB   := data/pipeline/lib
PIPRES   := data/pipeline/results

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
	@GENOME_SIZE=$$(grep -v ">" $(T1DIR)/genome/t1_genome.fa | tr -d '\n' | wc -c); \
	N_PAIRS=$$(( GENOME_SIZE * 30 / 300 )); \
	echo "Simulating $$N_PAIRS read pairs at 30x from T1 genome..."; \
	mkdir -p $(T1DIR)/reads; \
	wgsim -N $$N_PAIRS -1 150 -2 150 -d 500 -s 50 -e 0.01 -r 0.0 -R 0.0 -S 42 \
	  $(T1DIR)/genome/t1_genome.fa \
	  $(T1DIR)/reads/t1_reads1.fq \
	  $(T1DIR)/reads/t1_reads2.fq 2>&1 | tail -2
	@echo "Reads written to $(T1DIR)/reads/"

# ──────────────────────────────────────────────────
# Tier 2: 50 Mb genome, 5 families × 4 divergence bins
# ──────────────────────────────────────────────────
t2: lib
	$(PYTHON) scripts/sim/make_t2.py --libdir $(LIBDIR) --outdir $(T2DIR)

t2-reads-30x: t2
	@GENOME_SIZE=$$(grep -v ">" $(T2DIR)/genome/t2_genome.fa | tr -d '\n' | wc -c); \
	N_PAIRS=$$(( GENOME_SIZE * 30 / 300 )); \
	echo "Simulating $$N_PAIRS read pairs at 30x from T2 genome..."; \
	mkdir -p $(T2DIR)/reads_30x; \
	wgsim -N $$N_PAIRS -1 150 -2 150 -d 500 -s 50 -e 0.01 -r 0.0 -R 0.0 -S 42 \
	  $(T2DIR)/genome/t2_genome.fa \
	  $(T2DIR)/reads_30x/t2_reads_30x1.fq \
	  $(T2DIR)/reads_30x/t2_reads_30x2.fq 2>&1 | tail -2

t2-reads-10x: t2
	@GENOME_SIZE=$$(grep -v ">" $(T2DIR)/genome/t2_genome.fa | tr -d '\n' | wc -c); \
	N_PAIRS=$$(( GENOME_SIZE * 10 / 300 )); \
	echo "Simulating $$N_PAIRS read pairs at 10x from T2 genome..."; \
	mkdir -p $(T2DIR)/reads_10x; \
	wgsim -N $$N_PAIRS -1 150 -2 150 -d 500 -s 50 -e 0.01 -r 0.0 -R 0.0 -S 42 \
	  $(T2DIR)/genome/t2_genome.fa \
	  $(T2DIR)/reads_10x/t2_reads_10x1.fq \
	  $(T2DIR)/reads_10x/t2_reads_10x2.fq 2>&1 | tail -2

# ──────────────────────────────────────────────────
# Pipeline: prepare library
# ──────────────────────────────────────────────────
prepare: lib
	@command -v bwa-mem2 >/dev/null 2>&1 || \
	  (echo "bwa-mem2 not found — run: nix shell nixpkgs#bwa-mem2" && exit 1)
	@command -v samtools >/dev/null 2>&1 || \
	  (echo "samtools not found — run: nix shell nixpkgs#samtools" && exit 1)
	python3 -m scripts.pipeline.cli prepare \
	  --ltr $(LIBDIR)/ltr_only.fa \
	  --internal $(LIBDIR)/internal_only.fa \
	  --out-dir $(PIPLIB)
	@echo "Library prepared: $(PIPLIB)/"

# ──────────────────────────────────────────────────
# Pipeline: run on T1 reads
# ──────────────────────────────────────────────────
run-t1: prepare t1-reads
	@mkdir -p $(PIPRES)/t1
	python3 -m scripts.pipeline.cli run \
	  --lib-dir $(PIPLIB) \
	  --r1 $(T1DIR)/reads/t1_reads1.fq \
	  --r2 $(T1DIR)/reads/t1_reads2.fq \
	  --out $(PIPRES)/t1 \
	  --sample-id t1_test \
	  --bootstrap

# ──────────────────────────────────────────────────
# Validate T1: check S/C ≈ 10.0 (50 solo / 5 intact)
# ──────────────────────────────────────────────────
validate-t1: run-t1
	@echo "=== T1 Validation (expected S/C ≈ 10.0 for RIRE2) ==="
	@awk -F'\t' 'NR>1 && $$1=="RIRE2" { \
	  sc=$$3; r=$$2; \
	  printf "  RIRE2: R=%.4f  S/C=%.4f  call=%s  PAV=%s  confidence=%s\n", r, sc, $$4, $$5, $$12; \
	  if (sc >= 9.0 && sc <= 11.0) print "  PASS: S/C within 10%% of expected 10.0"; \
	  else print "  FAIL: S/C out of range [9.0, 11.0]" }' \
	  $(PIPRES)/t1/ratios.tsv

# ──────────────────────────────────────────────────
# Pipeline: run on T2 reads (30x, 10x)
# ──────────────────────────────────────────────────
run-t2-30x: prepare t2-reads-30x
	@mkdir -p $(PIPRES)/t2_30x
	python3 -m scripts.pipeline.cli run \
	  --lib-dir $(PIPLIB) \
	  --r1 $(T2DIR)/reads_30x/t2_reads_30x1.fq \
	  --r2 $(T2DIR)/reads_30x/t2_reads_30x2.fq \
	  --out $(PIPRES)/t2_30x \
	  --sample-id t2_30x \
	  --bootstrap

run-t2-10x: prepare t2-reads-10x
	@mkdir -p $(PIPRES)/t2_10x
	python3 -m scripts.pipeline.cli run \
	  --lib-dir $(PIPLIB) \
	  --r1 $(T2DIR)/reads_10x/t2_reads_10x1.fq \
	  --r2 $(T2DIR)/reads_10x/t2_reads_10x2.fq \
	  --out $(PIPRES)/t2_10x \
	  --sample-id t2_10x \
	  --bootstrap

# ──────────────────────────────────────────────────
# Validate T2: check S/C ≈ 6.67 per family
# ──────────────────────────────────────────────────
validate-t2-30x: run-t2-30x
	@echo "=== T2 30x Validation ==="
	@echo "    True S/C = 6.67 (40 solo / 6 intact per family)"
	@echo "    Detectable S/C ≈ 3.3–4.5 (only 0–15%% divergence copies map at MAPQ≥1)"
	@echo "    Pass criteria: all families detected, S/C in [2.5, 6.0], D_INT ≈ 180x"
	@echo ""
	@awk -F'\t' 'NR>1 { \
	  sc=$$3; r=$$2; dltr=$$6; dint=$$7; \
	  if (sc != "NA") { \
	    det_ok = ($$5 == "1") ? "detected" : "NOT_detected"; \
	    sc_ok = (sc+0 >= 2.5 && sc+0 <= 6.0) ? "PASS" : "FAIL"; \
	    int_ok = (dint+0 >= 150 && dint+0 <= 210) ? "INT_ok" : "INT_warn"; \
	    printf "  %-10s R=%.4f  S/C=%.2f  %s  %s  D_LTR=%.0f  D_INT=%.0f  %s\n", \
	      $$1, r, sc, sc_ok, det_ok, dltr, dint, int_ok \
	  } else printf "  %-10s %s\n", $$1, $$4 }' \
	  $(PIPRES)/t2_30x/ratios.tsv

validate-t2-10x: run-t2-10x
	@echo "=== T2 10x Validation ==="
	@echo "    Pass criteria: all families detected, S/C in [2.5, 6.0], D_INT ≈ 60x"
	@awk -F'\t' 'NR>1 { \
	  sc=$$3; r=$$2; dltr=$$6; dint=$$7; \
	  if (sc != "NA") { \
	    det_ok = ($$5 == "1") ? "detected" : "NOT_detected"; \
	    sc_ok = (sc+0 >= 2.5 && sc+0 <= 6.0) ? "PASS" : "FAIL"; \
	    int_ok = (dint+0 >= 40 && dint+0 <= 80) ? "INT_ok" : "INT_warn"; \
	    printf "  %-10s R=%.4f  S/C=%.2f  %s  %s  D_LTR=%.0f  D_INT=%.0f  %s\n", \
	      $$1, r, sc, sc_ok, det_ok, dltr, dint, int_ok \
	  } else printf "  %-10s %s\n", $$1, $$4 }' \
	  $(PIPRES)/t2_10x/ratios.tsv

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
        check-t0 check-t1 check-t2 evaluate-t1 evaluate-t2 \
        prepare run-t1 validate-t1 run-t2-30x run-t2-10x validate-t2-30x validate-t2-10x
