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

# ──────────────────────────────────────────────────
# Semi-synthetic benchmark
# Phase A: annotation (run once per real species)
# Phase B: catalog -> eval (automated, repeatable)
# ──────────────────────────────────────────────────
# Override SS_SPECIES to switch species, e.g.:
#   make semisyn-merge SS_SPECIES=arabidopsis SS_PREFIX=at
SS_SPECIES   := t2_mock
SS_PREFIX    := sp          # short species ID prefix for locus IDs
SS_CATALOG   := data/semisyn/catalogs/$(SS_SPECIES)_catalog.tsv
SS_GENOME    := data/sim/t2/genome/t2_genome.fa
SS_EXTRACTED := data/semisyn/extracted/$(SS_SPECIES)
SS_SAMPLES   := data/semisyn/samples/$(SS_SPECIES)
SS_RESULTS   := data/semisyn/results/$(SS_SPECIES)
SS_EVAL      := data/semisyn/eval/$(SS_SPECIES)
SS_LIB       := data/semisyn/lib/$(SS_SPECIES)
SS_COVERAGES := 5 10 30
SS_NSAMPLES  := 8
SS_BG_GC     := 0.43
SS_THREADS   := 16

# Phase A annotation paths (set SS_ANNOTATIONS to your tool output directory)
SS_ANNOTATIONS  := data/semisyn/annotations/$(SS_SPECIES)

# ── Per-species overrides ──────────────────────────────────────────────────
# When SS_SPECIES=arabidopsis, EDTA names outputs after the input genome file
ifeq ($(SS_SPECIES),arabidopsis)
  SS_PREFIX      := at
  SS_GENOME      := data/semisyn/genomes/arabidopsis/TAIR10.fa
  SS_BG_GC       := 0.36
  SS_EDTA_INTACT := $(SS_ANNOTATIONS)/TAIR10.fa.mod.EDTA.intact.gff3
  SS_EDTA_ANNO   := $(SS_ANNOTATIONS)/TAIR10.fa.mod.EDTA.TEanno.gff3
  SS_EDTA_LIB    := $(SS_ANNOTATIONS)/TAIR10.fa.mod.EDTA.TElib.fa
  # soloLTRseeker needs simple chromosome names; use EDTA's .mod genome (already in annotations dir)
  SS_SLTR_GENOME := $(SS_ANNOTATIONS)/TAIR10.fa.mod
else
  SS_EDTA_INTACT := $(SS_ANNOTATIONS)/$(SS_SPECIES).EDTA.intact.gff3
  SS_EDTA_ANNO   := $(SS_ANNOTATIONS)/$(SS_SPECIES).EDTA.TEanno.gff3
  SS_EDTA_LIB    := $(SS_ANNOTATIONS)/$(SS_SPECIES).EDTA.TElib.fa
  SS_SLTR_GENOME := $(SS_GENOME)
endif
SS_L4LTR_DIR    := $(SS_ANNOTATIONS)/look4ltrs
SS_SLTR_SOLOS   := $(SS_ANNOTATIONS)/soloLTRseeker_solos.gff3

# ── Phase A targets (manual, once per real species) ──────────────────────

# Parse EDTA intact.gff3 + TEanno.gff3 -> per-tool TSV
semisyn-parse-edta:
	$(PYTHON) scripts/semisyn/parse_edta.py \
	  --intact-gff3 $(SS_EDTA_INTACT) \
	  --anno-gff3 $(SS_EDTA_ANNO) \
	  --species-prefix $(SS_PREFIX) \
	  --out $(SS_ANNOTATIONS)/edta_loci.tsv

# Parse Look4LTRs results directory -> per-tool TSV
semisyn-parse-look4ltrs:
	$(PYTHON) scripts/semisyn/parse_look4ltrs.py \
	  --results-dir $(SS_L4LTR_DIR) \
	  --species-prefix $(SS_PREFIX) \
	  --out $(SS_ANNOTATIONS)/look4ltrs_loci.tsv

# Parse soloLTRseeker GFF3 (+ EDTA intact for intact calls) -> per-tool TSV
semisyn-parse-soloLTRseeker:
	$(PYTHON) scripts/semisyn/parse_soloLTRseeker.py \
	  --solo-gff3 $(SS_SLTR_SOLOS) \
	  --intact-gff3 $(SS_EDTA_INTACT) \
	  --intact-format edta \
	  --species-prefix $(SS_PREFIX) \
	  --out $(SS_ANNOTATIONS)/soloLTRseeker_loci.tsv

# 3-tool consensus merge -> final catalog TSV
# Look4LTRs is optional — pass --look4ltrs only if the TSV exists (e.g. not available for arabidopsis)
SS_LOOK4LTRS_TSV := $(SS_ANNOTATIONS)/look4ltrs_loci.tsv
semisyn-merge: semisyn-parse-edta semisyn-parse-soloLTRseeker
	$(PYTHON) scripts/semisyn/merge_annotations.py \
	  --edta            $(SS_ANNOTATIONS)/edta_loci.tsv \
	  $(if $(wildcard $(SS_LOOK4LTRS_TSV)),--look4ltrs $(SS_LOOK4LTRS_TSV),) \
	  --soloLTRseeker   $(SS_ANNOTATIONS)/soloLTRseeker_loci.tsv \
	  --species-prefix  $(SS_PREFIX) \
	  --species         $(SS_SPECIES) \
	  --out $(SS_CATALOG)

# Build per-species consensus library from EDTA TElib.fa
semisyn-lib: semisyn-merge
	$(PYTHON) scripts/semisyn/build_consensus_lib.py \
	  --edta-lib $(SS_EDTA_LIB) \
	  --catalog  $(SS_CATALOG) \
	  --out-dir  $(SS_LIB)

# Generate mock catalog from T2 ground truth (dev/test only)
semisyn-mock-catalog:
	$(PYTHON) scripts/semisyn/mock_catalog_from_t2.py \
	  --truth data/sim/t2/ground_truth/t2_ground_truth.json \
	  --regions data/pipeline/lib/regions.tsv \
	  --out $(SS_CATALOG)

# Extract per-locus sequences from genome
semisyn-extract: $(SS_CATALOG)
	$(PYTHON) scripts/semisyn/extract_sequences.py \
	  --catalog $(SS_CATALOG) \
	  --genome $(SS_GENOME) \
	  --flank 700 \
	  --out-dir $(SS_EXTRACTED)

# Partition loci into N samples with target S/C ratios
semisyn-partition: semisyn-extract
	$(PYTHON) scripts/semisyn/partition_samples.py \
	  --catalog $(SS_CATALOG) \
	  --n-samples $(SS_NSAMPLES) \
	  --seed 42 \
	  --species $(SS_SPECIES) \
	  --out $(SS_SAMPLES)/manifest.json

# Compose virtual genomes per sample
semisyn-genomes: semisyn-partition
	$(PYTHON) scripts/semisyn/compose_genome.py \
	  --manifest $(SS_SAMPLES)/manifest.json \
	  --loci-fa $(SS_EXTRACTED)/loci_with_flanking.fa \
	  --loci-meta $(SS_EXTRACTED)/loci_metadata.json \
	  --bg-gc $(SS_BG_GC) \
	  --out-dir $(SS_SAMPLES)

# Simulate reads at multiple coverages per sample
semisyn-reads: semisyn-genomes
	@for sample_dir in $(SS_SAMPLES)/sample_*/; do \
	  sid=$$(basename $$sample_dir); \
	  genome=$$sample_dir/genome.fa; \
	  gsize=$$(grep -v ">" $$genome | tr -d '\n' | wc -c); \
	  if [ "$$gsize" -lt 1000 ]; then \
	    echo "  $$sid: genome too small ($$gsize bp), skipping"; continue; fi; \
	  for cov in $(SS_COVERAGES); do \
	    npairs=$$(( gsize * cov / 300 )); \
	    mkdir -p $$sample_dir/reads_$${cov}x; \
	    wgsim -N $$npairs -1 150 -2 150 -d 500 -s 50 -e 0.01 -r 0.0 -R 0.0 -S 42 \
	      $$genome \
	      $$sample_dir/reads_$${cov}x/r1.fq \
	      $$sample_dir/reads_$${cov}x/r2.fq 2>/dev/null; \
	    echo "  $$sid $${cov}x: $$npairs pairs"; \
	  done; \
	done

# Prepare consensus library for pipeline (uses synthetic lib for mock)
semisyn-prepare: prepare
	@mkdir -p $(SS_LIB)
	@cp data/pipeline/lib/collapsed_consensus.fa $(SS_LIB)/collapsed_consensus.fa
	@cp data/pipeline/lib/regions.tsv $(SS_LIB)/regions.tsv
	@cp -r data/pipeline/lib/collapsed_consensus.fa.* $(SS_LIB)/ 2>/dev/null || true
	@echo "Using synthetic consensus library for mock benchmark"
	$(PYTHON) scripts/semisyn/mock_catalog_from_t2.py \
	  --truth data/sim/t2/ground_truth/t2_ground_truth.json \
	  --regions data/pipeline/lib/regions.tsv \
	  --out $(SS_CATALOG) 2>/dev/null || true

# Run pipeline on all samples x coverages
# NOTE: depends only on reads; caller must ensure library is prepared first
# (semisyn-all calls semisyn-prepare first; semisyn-all-real calls semisyn-prepare-real)
semisyn-run: semisyn-reads
	@for sample_dir in $(SS_SAMPLES)/sample_*/; do \
	  sid=$$(basename $$sample_dir); \
	  for cov in $(SS_COVERAGES); do \
	    reads_dir=$$sample_dir/reads_$${cov}x; \
	    [ -f "$$reads_dir/r1.fq" ] || continue; \
	    outdir=$(SS_RESULTS)/$${sid}_$${cov}x; \
	    mkdir -p $$outdir; \
	    python3 -m scripts.pipeline.cli run \
	      --lib-dir $(SS_LIB) \
	      --r1 $$reads_dir/r1.fq \
	      --r2 $$reads_dir/r2.fq \
	      --out $$outdir \
	      --sample-id $${sid}_$${cov}x \
	      --bootstrap 2>/dev/null; \
	    echo "  $$sid $${cov}x: done"; \
	  done; \
	done

# Evaluate: compare pipeline output vs ground truth
semisyn-eval:
	$(PYTHON) scripts/semisyn/evaluate_semisyn.py \
	  --manifest $(SS_SAMPLES)/manifest.json \
	  --results-dir $(SS_RESULTS) \
	  --coverages $(shell echo $(SS_COVERAGES) | tr ' ' ',') \
	  --out-dir $(SS_EVAL)

# Full Phase B pipeline (catalog -> eval) — mock/t2_mock only
semisyn-all: semisyn-mock-catalog semisyn-extract semisyn-partition \
             semisyn-genomes semisyn-reads semisyn-prepare semisyn-run semisyn-eval

# Full Phase B pipeline for real species (catalog already built by Phase A)
semisyn-all-real: semisyn-extract semisyn-partition semisyn-genomes \
                  semisyn-reads semisyn-prepare-real semisyn-run semisyn-eval

# ── Tool installation (one-time setup) ───────────────────────────────────

LOOK4LTRS_BIN := tools/Look4LTRs/build/look4ltrs
SLTR_REPO_URL  := https://github.com/estpr/soloLTRseeker.git
L4LTR_REPO_URL := https://github.com/BioinformaticsToolsmith/Look4LTRs.git

# Pull EDTA Docker image (bundles LTR_retriever, RepeatMasker, LTRharvest)
semisyn-install-edta:
	docker pull oushujun/edta:2.0.0
	@docker run --rm oushujun/edta:2.0.0 EDTA.pl --help 2>&1 | head -3

# Clone and build Look4LTRs from source (C++17, no external deps)
$(LOOK4LTRS_BIN):
	mkdir -p tools
	git clone $(L4LTR_REPO_URL) tools/Look4LTRs
	cd tools/Look4LTRs && mkdir -p build && cd build && cmake .. -DCMAKE_C_COMPILER=/usr/bin/gcc-13 -DCMAKE_CXX_COMPILER=/usr/bin/g++-13 && make -j$$(nproc)
	@echo "Look4LTRs built: $(LOOK4LTRS_BIN)"

semisyn-install-look4ltrs: $(LOOK4LTRS_BIN)

# Clone soloLTRseeker and build Docker image (bundles BLAST+, cd-hit, EMBOSS water, BEDTools)
semisyn-install-soloLTRseeker:
	@[ -d tools/soloLTRseeker ] || git clone $(SLTR_REPO_URL) tools/soloLTRseeker
	cp scripts/semisyn/soloLTRseeker.Dockerfile tools/soloLTRseeker/Dockerfile
	docker build -t sololtrs:latest tools/soloLTRseeker
	@echo "soloLTRseeker Docker image ready"

# Install all 3 annotation tools
semisyn-install-tools: semisyn-install-edta semisyn-install-look4ltrs semisyn-install-soloLTRseeker

# ── Genome download ───────────────────────────────────────────────────────

SS_GENOME_URL_AT := https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz

semisyn-genome-arabidopsis: data/semisyn/genomes/arabidopsis/TAIR10.fa.fai

data/semisyn/genomes/arabidopsis/TAIR10.fa.fai: data/semisyn/genomes/arabidopsis/TAIR10.fa
	samtools faidx $<

data/semisyn/genomes/arabidopsis/TAIR10.fa:
	mkdir -p data/semisyn/genomes/arabidopsis
	curl -k -L $(SS_GENOME_URL_AT) | gunzip > $@
	@echo "Downloaded TAIR10: $$(wc -c < $@) bytes"
	@echo "Chromosomes: $$(grep -c '>' $@)"

# ── Phase A: run annotation tools ────────────────────────────────────────

# EDTA: copy genome into annotations dir (EDTA writes output next to input)
semisyn-run-edta: $(SS_GENOME)
	mkdir -p $(SS_ANNOTATIONS)
	@[ -f $(SS_ANNOTATIONS)/$$(basename $(SS_GENOME)) ] || \
	  cp $(SS_GENOME) $(SS_ANNOTATIONS)/
	docker run --rm \
	  -v $$(realpath $(SS_ANNOTATIONS)):/work -w /work \
	  oushujun/edta:2.0.0 \
	  EDTA.pl \
	    --genome /work/$$(basename $(SS_GENOME)) \
	    --species others \
	    --anno 1 \
	    --threads $(SS_THREADS)
	@echo "EDTA done. Key outputs:"
	@ls -lh $(SS_ANNOTATIONS)/*.EDTA.intact.gff3 $(SS_ANNOTATIONS)/*.EDTA.TEanno.gff3 \
	         $(SS_ANNOTATIONS)/*.EDTA.TElib.fa 2>/dev/null || true

# Look4LTRs: can run in parallel with EDTA
semisyn-run-look4ltrs: $(LOOK4LTRS_BIN) $(SS_GENOME)
	mkdir -p $(SS_L4LTR_DIR)
	$(LOOK4LTRS_BIN) \
	  -f $$(dirname $$(realpath $(SS_GENOME))) \
	  -o $$(realpath $(SS_L4LTR_DIR))/ \
	  -pa $(SS_THREADS)
	@echo "Look4LTRs done. Outputs:"
	@ls -lh $(SS_L4LTR_DIR)/*.gff3 2>/dev/null || ls -lh $(SS_L4LTR_DIR)/ 2>/dev/null

# soloLTRseeker: must run after EDTA (needs EDTA intact GFF3 as annotation input)
# CLI: soloLTRseeker ann_file fasta_file
# Output: written to sample_soloLTR/output_DATE_N/ subdirectory, named sample_wga_soloLTR.gff3
# Note: SS_SLTR_GENOME must have simple chromosome names (single-token FASTA headers)
semisyn-run-soloLTRseeker: $(SS_GENOME)
	@[ -f $(SS_EDTA_INTACT) ] || \
	  (echo "ERROR: EDTA intact GFF3 not found at $(SS_EDTA_INTACT) — run semisyn-run-edta first" && exit 1)
	@[ -f $(SS_SLTR_GENOME) ] || \
	  (echo "ERROR: soloLTRseeker genome not found at $(SS_SLTR_GENOME)" && exit 1)
	mkdir -p $(SS_ANNOTATIONS)
	@# soloLTRseeker writes output to sample_soloLTR/output_DATE_N/ subdirectory
	docker run --rm \
	  -v $$(realpath $(SS_ANNOTATIONS)):/work -w /work \
	  sololtrs:latest \
	    $$(basename $(SS_EDTA_INTACT)) $$(basename $(SS_SLTR_GENOME))
	@# Find output GFF3 (named sample_wga_soloLTR.gff3) and copy to canonical path
	@outfile=$$(find $(SS_ANNOTATIONS) -name "*_wga_soloLTR.gff3" -newer $(SS_EDTA_INTACT) 2>/dev/null | head -1); \
	  if [ -n "$$outfile" ]; then \
	    cp "$$outfile" "$(SS_SLTR_SOLOS)"; \
	    echo "soloLTRseeker done: $(SS_SLTR_SOLOS)"; \
	    echo "  Source: $$outfile"; \
	  else \
	    echo "ERROR: soloLTR GFF3 not found. Check logs in $(SS_ANNOTATIONS)/*/output_*/out.log"; \
	    exit 1; \
	  fi

# ── Phase B: prepare real-species library (replaces semisyn-prepare for real genomes) ──

semisyn-prepare-real: $(SS_LIB)/ltr_only.fa
	@command -v bwa-mem2 >/dev/null 2>&1 || \
	  (echo "bwa-mem2 not found — run inside: nix shell nixpkgs#bwa-mem2" && exit 1)
	$(PYTHON) -m scripts.pipeline.cli prepare \
	  --ltr $(SS_LIB)/ltr_only.fa \
	  --internal $(SS_LIB)/internal_only.fa \
	  --out-dir $(SS_LIB)
	@echo "Library prepared: $(SS_LIB)/"

$(SS_LIB)/ltr_only.fa: semisyn-lib
	@test -f $@ || (echo "ERROR: semisyn-lib did not produce $@" && exit 1)

# ── Catalog sanity check ─────────────────────────────────────────────────

semisyn-validate-catalog:
	@echo "=== Catalog: $(SS_CATALOG) ==="
	@[ -f $(SS_CATALOG) ] || (echo "ERROR: catalog not found" && exit 1)
	@echo "Total loci (excluding header):"
	@tail -n +2 $(SS_CATALOG) | wc -l
	@echo "Solo LTRs:"
	@awk -F'\t' '$$2=="solo_LTR"' $(SS_CATALOG) | wc -l
	@echo "Intact LTR-RTs:"
	@awk -F'\t' '$$2=="intact_LTR_RT"' $(SS_CATALOG) | wc -l
	@echo "Superfamilies:"
	@awk -F'\t' 'NR>1{print $$4}' $(SS_CATALOG) | sort | uniq -c | sort -rn
	@echo "n_tools distribution:"
	@awk -F'\t' 'NR>1{print $$15}' $(SS_CATALOG) | sort | uniq -c
	@echo "Per-family S/C:"
	@awk -F'\t' 'NR>1{fam=$$3; type=$$2; \
	    s[fam]+=(type=="solo_LTR"); i[fam]+=(type=="intact_LTR_RT")} \
	  END{for(f in s) if(i[f]>0) \
	    printf "  %-30s  S=%d  I=%d  S/C=%.2f\n", f, s[f], i[f], s[f]/i[f]}' \
	  $(SS_CATALOG) | sort

.PHONY: lib t0 t1 t1-reads t2 t2-reads-30x t2-reads-10x all-genomes all-reads \
        check-t0 check-t1 check-t2 evaluate-t1 evaluate-t2 \
        prepare run-t1 validate-t1 run-t2-30x run-t2-10x validate-t2-30x validate-t2-10x \
        semisyn-parse-edta semisyn-parse-look4ltrs semisyn-parse-soloLTRseeker \
        semisyn-merge semisyn-lib \
        semisyn-mock-catalog semisyn-extract semisyn-partition semisyn-genomes \
        semisyn-reads semisyn-prepare semisyn-run semisyn-eval semisyn-all semisyn-all-real \
        semisyn-install-edta semisyn-install-look4ltrs semisyn-install-soloLTRseeker \
        semisyn-install-tools \
        semisyn-genome-arabidopsis \
        semisyn-run-edta semisyn-run-look4ltrs semisyn-run-soloLTRseeker \
        semisyn-prepare-real semisyn-validate-catalog
