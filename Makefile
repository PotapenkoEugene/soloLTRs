PYTHON   := python3
THREADS  := 8

# ──────────────────────────────────────────────────────────────────────────────
# soloLTRs: Cossu M/U pipeline
#
# Goal: reproduce and then accelerate Cossu et al. 2017 (GBE, DOI: 10.1093/gbe/evx260)
#   Phase 1 — Reproduce: BLAST-based tag search + BWA ALN tract mapping
#   Phase 2 — Accelerate: replace BLAST with even faster alternative if needed
# ──────────────────────────────────────────────────────────────────────────────

# Paths
COSSU_REF      := references/cossu
AT_PARALOGS    := $(COSSU_REF)/Supporting_Data_S1/ANGIOSPERMS/complete_Arabidopsis.fa
AT_READS_R1    := data/reads/arabidopsis/ERR171441_1.fastq.gz
AT_READS_R2    := data/reads/arabidopsis/ERR171441_2.fastq.gz
AT_OUT         := data/cossu/arabidopsis

# ──────────────────────────────────────────────────
# Cossu M/U pipeline — Arabidopsis (Phase 1)
# ──────────────────────────────────────────────────

## mu-arabidopsis: run full M/U pipeline on Arabidopsis ERR171441 R1+R2 reads
mu-arabidopsis: $(AT_READS_R1) $(AT_READS_R2)
	mkdir -p $(AT_OUT)
	$(PYTHON) -m scripts.pipeline.cli all \
	  --paralogs $(AT_PARALOGS) \
	  --reads    $(AT_READS_R1) $(AT_READS_R2) \
	  --out-dir  $(AT_OUT) \
	  --threads  $(THREADS)

## mu-prepare-arabidopsis: extract tags only (no reads needed)
mu-prepare-arabidopsis:
	mkdir -p $(AT_OUT)
	$(PYTHON) -m scripts.pipeline.cli prepare \
	  --paralogs $(AT_PARALOGS) \
	  --out-dir  $(AT_OUT)

## mu-run-arabidopsis: run pipeline (assumes prepare already done)
mu-run-arabidopsis: $(AT_READS_R1) $(AT_READS_R2)
	$(PYTHON) -m scripts.pipeline.cli run \
	  --paralogs $(AT_PARALOGS) \
	  --reads    $(AT_READS_R1) $(AT_READS_R2) \
	  --out-dir  $(AT_OUT) \
	  --threads  $(THREADS)

## mu-clean-arabidopsis: remove all intermediate files (keep only results.tsv)
mu-clean-arabidopsis:
	rm -f $(AT_OUT)/reads.fa $(AT_OUT)/blast_hits.tsv
	rm -rf $(AT_OUT)/reads_db
	rm -f $(AT_OUT)/tracts.fa $(AT_OUT)/tracts.bam $(AT_OUT)/tracts.bam.bai
	rm -f $(AT_OUT)/hits.tsv

.PHONY: mu-arabidopsis mu-prepare-arabidopsis mu-run-arabidopsis mu-clean-arabidopsis

# ──────────────────────────────────────────────────
# Arabidopsis genome (TAIR10) — still needed for Phase 1
# ──────────────────────────────────────────────────

SS_GENOME_URL_AT := https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
AT_GENOME        := data/genomes/arabidopsis/TAIR10.fa
AT_GENOME_FAI    := $(AT_GENOME).fai

genome-arabidopsis: $(AT_GENOME_FAI)

$(AT_GENOME_FAI): $(AT_GENOME)
	samtools faidx $<

$(AT_GENOME):
	mkdir -p data/genomes/arabidopsis
	curl -k -L $(SS_GENOME_URL_AT) | gunzip > $@
	@echo "Downloaded TAIR10: $$(wc -c < $@) bytes"
	@echo "Chromosomes: $$(grep -c '>' $@)"

# ──────────────────────────────────────────────────
# Annotation tools (keep: EDTA output for arabidopsis still in
# data/semisyn/annotations/arabidopsis/ and may be reused)
# ──────────────────────────────────────────────────

ANNOTATIONS     := data/semisyn/annotations
SS_THREADS      := 16

semisyn-run-edta:
	@[ -n "$(SS_SPECIES)" ] || (echo "Usage: make semisyn-run-edta SS_SPECIES=arabidopsis" && exit 1)
	mkdir -p $(ANNOTATIONS)/$(SS_SPECIES)
	@[ -f $(ANNOTATIONS)/$(SS_SPECIES)/$$(basename $(SS_GENOME)) ] || \
	  cp $(SS_GENOME) $(ANNOTATIONS)/$(SS_SPECIES)/
	docker run --rm \
	  -v $$(realpath $(ANNOTATIONS)/$(SS_SPECIES)):/work -w /work \
	  oushujun/edta:2.0.0 \
	  EDTA.pl \
	    --genome /work/$$(basename $(SS_GENOME)) \
	    --species others \
	    --anno 1 \
	    --threads $(SS_THREADS)

.PHONY: genome-arabidopsis semisyn-run-edta mu-arabidopsis mu-prepare-arabidopsis mu-run-arabidopsis mu-clean-arabidopsis
