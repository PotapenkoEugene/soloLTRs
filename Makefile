PYTHON   := python3

# ──────────────────────────────────────────────────────────────────────────────
# soloLTRs: Cossu M/U pipeline
#
# Goal: reproduce and then accelerate Cossu et al. 2017 (GBE, DOI: 10.1093/gbe/evx260)
#   Phase 1 — Reproduce: RepeatMasker-based M/U pipeline
#   Phase 2 — Accelerate: replace RepeatMasker with fast aligner
# ──────────────────────────────────────────────────────────────────────────────

# Paths
COSSU_REF   := references/cossu
DATA_COSSU  := data/cossu

# ──────────────────────────────────────────────────
# Cossu M/U pipeline (TODO: implement in Phase 1)
# ──────────────────────────────────────────────────

# mu-prepare: extract START/END tags + 20-nt tracts from paralogs
# mu-run:     RepeatMasker read search → tract extraction → BWA ALN → M/U count
# mu-fast:    same pipeline with fast aligner replacing RepeatMasker

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

.PHONY: genome-arabidopsis semisyn-run-edta
