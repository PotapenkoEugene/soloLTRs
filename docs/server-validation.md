# Server Plan: Full Multi-Species Validation (Phase 3)

Run the soloLTRs M/U pipeline (minimap2 engine) on all 11 Cossu et al. 2017 species to reproduce Table S2/S3 S/C ratios. Arabidopsis is included here for completeness — its reads and validated result are in the repo if you want to re-run, but the result (S/C=0.862 vs paper 0.865) is already confirmed.

## Species overview

| # | Species | Key | SRA | Reads (Gbp) | Coverage | Expected S/C |
|---|---------|-----|-----|------------|----------|-------------|
| 1 | Arabidopsis thaliana | arabidopsis | ERR171441 | 1.18 | 7.56x | 0.86 |
| 2 | Amborella trichopoda | amborella | SRR954506 | 1.4 | 1.61x | 1.38 |
| 3 | Brachypodium distachyon | brachypodium | SRR494350 | 1.6 | 4.51x | 1.31 |
| 4 | Oryza sativa | oryza | SRR1450197 | 1.4 | 2.86x | 1.78 |
| 5 | Populus trichocarpa | populus | SRR1302874 | 1.46 | 3.02x | 0.95 |
| 6 | Vitis vinifera | vitis | SRR931846 | 1.4 | 3.37x | 1.28 |
| 7 | Zea mays | zea | ERR361368 | 1.4 | 0.53x | 0.15 |
| 8 | Picea abies | picea_abies | ERR268359 | 3.94 | 0.20x | 0.18 |
| 9 | Picea glauca | picea_glauca | SRR1259601 | 4.30 | 0.22x | 0.24 |
| 10 | Pinus taeda | pinus_taeda | SRR1049544 | 5.04 | 0.23x | 0.13 |
| 11 | Pinus lambertiana | pinus_lambertiana | SRR2027098 | 5.04 | 0.16x | 0.27 |

**Total download:** ~28 GB compressed. **Disk needed:** ~80 GB (downloads + uncompressed tmp + results).

All reads originally from ENA. Both SRR and ERR accessions work with `fasterq-dump`.

---

## Prerequisites

### Software

```bash
# sra-tools (fasterq-dump)
nix shell nixpkgs#sratoolkit          # NixOS / nix
# OR: conda install -c bioconda sra-tools
# OR: sudo apt install sra-toolkit   # Ubuntu

# Pipeline tools (minimap2, bwa, samtools, python3)
# Already in the nix environment — or install individually:
nix shell nixpkgs#minimap2 nixpkgs#bwa nixpkgs#samtools nixpkgs#python3
```

### Repo

```bash
# Clone fresh:
git clone <repo-url> soloLTRs && cd soloLTRs

# OR sync from your desktop (skip large read files):
rsync -av --exclude='data/reads/' desktop:~/Documents/Projects/soloLTRs/ soloLTRs/
cd soloLTRs
```

### Conifer paralog FASTAs

Already committed at `references/cossu/Supporting_Data_S1/CONIFERS/`. If for some reason missing:

```bash
make prepare-conifers
```

---

## Step 1: Download reads

Each species is independent — parallelise freely. The `make download-{species}` targets are idempotent (skip if `{SRA}_1.fastq.gz` already exists).

### All at once (if bandwidth/disk allows)

```bash
make -j4 download-all
```

### Species by species (recommended for sequential / limited disk)

```bash
# Angiosperms (~1.4–1.6 Gbp each, ~3–5 GB compressed)
make download-arabidopsis
make download-amborella
make download-brachypodium
make download-oryza
make download-populus
make download-vitis
make download-zea

# Conifers (~3.9–5.0 Gbp each, ~8–12 GB compressed)
make download-picea_abies
make download-picea_glauca
make download-pinus_taeda
make download-pinus_lambertiana
```

Each target runs:
```
fasterq-dump {SRA} --outdir data/reads/{species}/ --split-files --threads 8 --temp data/reads/{species}/
gzip data/reads/{species}/{SRA}_1.fastq  data/reads/{species}/{SRA}_2.fastq
```

**Temp space note:** `fasterq-dump` writes uncompressed FASTQ to `--temp` before gzipping. Peak disk per species = ~3× the final compressed size. To use a different temp dir (e.g. `/scratch`):

```bash
# Override for a single species manually:
fasterq-dump ERR268359 --outdir data/reads/picea_abies/ --split-files --threads 8 --temp /scratch/tmp/
gzip data/reads/picea_abies/ERR268359_1.fastq data/reads/picea_abies/ERR268359_2.fastq
```

---

## Step 2: Run pipeline

Each species is independent. Start pipeline runs as soon as each species' reads are ready — no need to wait for all downloads.

### All at once (after all downloads)

```bash
make mu-all
```

### Species by species

```bash
make mu-arabidopsis-mm2
make mu-amborella-mm2
make mu-brachypodium-mm2
make mu-oryza-mm2
make mu-populus-mm2
make mu-vitis-mm2
make mu-zea-mm2
make mu-picea_abies-mm2
make mu-picea_glauca-mm2
make mu-pinus_taeda-mm2
make mu-pinus_lambertiana-mm2
```

Each run uses minimap2 engine, 8 threads. Expected runtimes:
- Angiosperms (1.4–1.6 Gbp reads): ~1–2 min each
- Zea mays (1.4 Gbp but large M/U counts): ~2–3 min
- Conifers (3.9–5.0 Gbp reads): ~5–15 min each

Results land at `data/cossu/{species}-mm2/results.tsv`.

### Interleave downloads and runs (recommended)

```bash
# In separate terminals / tmux panes:
# Pane 1: sequential downloads
for sp in arabidopsis amborella brachypodium oryza populus vitis zea picea_abies picea_glauca pinus_taeda pinus_lambertiana; do
    make download-$sp
    make mu-${sp}-mm2   # run pipeline immediately after each download
done

# Pane 2: monitor validation as results arrive
watch -n 30 make mu-validate
```

---

## Step 3: Validate

```bash
make mu-validate
```

This runs `scripts/pipeline/validate.py` against `references/cossu/expected_ratios.tsv` (tolerance 15%) and prints:

```
Cossu et al. 2017 S/C validation (tolerance=15%)
-----------------------------------------------------------------------------------------------
species            sra           expected_SC  observed_SC  abs_error    rel_error%  result
-----------------------------------------------------------------------------------------------
arabidopsis        ERR171441           0.860        0.862      0.002          0.2%  PASS
amborella          SRR954506           1.380        ?            ?              ?     ?
...
-----------------------------------------------------------------------------------------------
Summary written to data/cossu/validation_summary.tsv
```

Also writes `data/cossu/validation_summary.tsv` for further analysis.

---

## Step 4: Copy results back to desktop

Only copy the small results files — not the reads.

```bash
# From desktop:
for sp in arabidopsis amborella brachypodium oryza populus vitis zea picea_abies picea_glauca pinus_taeda pinus_lambertiana; do
    mkdir -p data/cossu/${sp}-mm2
    rsync -av server:soloLTRs/data/cossu/${sp}-mm2/results.tsv  data/cossu/${sp}-mm2/
    rsync -av server:soloLTRs/data/cossu/${sp}-mm2/timings.tsv  data/cossu/${sp}-mm2/
done

# Then run final validation locally:
make mu-validate
```

---

## Troubleshooting

### fasterq-dump: temp space error

```bash
# Use a larger partition for temp:
fasterq-dump {SRA} --outdir data/reads/{species}/ --split-files --threads 8 --temp /scratch/
gzip data/reads/{species}/{SRA}_*.fastq
```

### Single-end reads (only `_1.fastq.gz` produced)

Some accessions may be single-end. If `_2.fastq.gz` is missing, run the pipeline with one reads file:

```bash
python3 -m scripts.pipeline.cli all \
  --paralogs {PARALOGS} \
  --reads data/reads/{species}/{SRA}_1.fastq.gz \
  --out-dir data/cossu/{species}-mm2 --threads 8 --engine minimap2
```

### S/C deviates >15% from expected

Check the TOTAL row in `results.tsv`. Families with M+U < 50 are noisy at low coverage (common for conifers at 0.16–0.23x). The TOTAL is what counts. If TOTAL also deviates, check:

1. Are reads actually paired-end? (`zcat {SRA}_1.fastq.gz | head -8`)
2. Do tag extraction stats look reasonable? (`cat data/cossu/{species}-mm2/timings.tsv`)
3. Is M+U in the right ballpark vs paper? (e.g. Picea abies should have ~146K total tracts)

### Pipeline caching

Each step checks for existing output and skips if present. To re-run from scratch:

```bash
rm -rf data/cossu/{species}-mm2
make mu-{species}-mm2
```
