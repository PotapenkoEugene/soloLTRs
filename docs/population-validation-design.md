# Research Design: Validating M/U Pipeline as Population-Scale Recombination Proxy

## Context

The minimap2-accelerated Cossu M/U pipeline estimates solo:intact LTR ratios (S/C) from WGS reads (~42s/sample for Arabidopsis). The goal: use per-sample S/C as a **recombination proxy across population samples at low coverage** -- populations/ecotypes with higher historical recombination should show higher S/C.

**Gap**: No study has operationalized solo:intact ratios from WGS reads as a population-scale recombination metric. The biology is established (Vitte & Panaud 2003, Du et al. 2012, Cossu 2017) but the inference framework doesn't exist.

**Focus**: Per-sample S/C variation across populations (each individual -> one S/C value -> compare distributions).

---

## Tier 1: Coverage Sensitivity -- "How Low Can You Go?"

**Rationale**: Before any population work, we need to know the minimum coverage for reliable per-sample S/C. The M/U ratio is theoretically coverage-independent, but variance increases as tract counts shrink.

### Experiment Design

**Controlled downsampling of Arabidopsis ERR171441** (existing data, 30.9M reads, ~20x):
- Downsample to: 15x, 10x, 5x, 2x, 1x, 0.5x using `seqtk sample` with fixed seed
- 10 replicates per coverage level (different random seeds) for variance estimation
- Run M/U pipeline at each level -> S/C per family + TOTAL

### Expected Outputs
- **S/C vs coverage curve** with 95% CI (bootstrap from replicate variance)
- **Per-family M+U tract counts** vs coverage -> minimum tract threshold for reliable estimates
- **Recommendation**: minimum coverage for TOTAL S/C within +/-10% of true value, and per-family threshold

### Key Questions Answered
- At 5x, is TOTAL S/C stable? (Cossu used 0.16-7.56x, suggesting yes)
- Do per-family estimates break before TOTAL? At what M+U threshold?
- Does BWA ALN l=15 calibration hold across coverages?

---

## Tier 2: Assembly-Level Ground Truth -- "Is the Pipeline Accurate Per-Individual?"

**Rationale**: Before correlating S/C with recombination, we need to confirm the pipeline accurately estimates solo:intact ratios for individual genomes (not just the reference).

### Dataset A: Arabidopsis 72 Pan-Genome Accessions

| Field | Value |
|-------|-------|
| BioProject | PRJEB62038 |
| Paper | Jiao & Schneeberger 2024 |
| N accessions | 72 |
| Data type | Complete assemblies (PacBio HiFi + ONT) + Illumina short reads |
| Genome size | 135 Mb |

**Why ideal**: Both assembly-level ground truth AND short reads exist for the **same individuals**. The pipeline is already validated on Arabidopsis (S/C=0.862 for ERR171441).

**Validation**:
1. Run EDTA on each assembly -> count solo vs intact LTRs per family -> true S/C per accession
2. Run M/U pipeline on Illumina reads for same accessions -> predicted S/C
3. Pearson/Spearman correlation of true vs predicted S/C across 72 accessions
4. Per-family breakdown: which families are well-estimated, which are noisy?

**Success criterion**: r > 0.8 between assembly-derived and pipeline-derived S/C

### Dataset B: Maize NAM 26 Assemblies (Hard Test)

| Field | Value |
|-------|-------|
| BioProject | PRJNA697943 |
| Paper | Hufford et al. 2021 (Science) |
| N accessions | 26 |
| Data type | Complete assemblies with TE annotations already done |
| Genome size | 2.3 Gb |

**Why useful**: Large genome, very low solo:intact ratio (S/C=0.15). Tests the method at its theoretical limit. TE annotations already exist per assembly -- no need to run EDTA.

**Challenge**: Need to identify matching Illumina accessions for each NAM line (likely in HapMap3 PRJNA389800 or WiDiv PRJNA661271).

---

## Tier 3: Population-Scale Recombination Correlation -- "Does S/C Track Recombination?"

**Focus**: Per-sample S/C variation. Each individual gets one TOTAL S/C value. Compare S/C distributions across populations/ecotypes/species that differ in known recombination characteristics.

### Dataset C: Arabidopsis 1001 Genomes -- Best Starting Point

| Field | Value |
|-------|-------|
| BioProject | PRJNA273563 |
| Paper | 1001 Genomes Consortium |
| N accessions | 1,135 |
| Coverage | Variable, typically 20-40x |
| Genome size | 135 Mb |

**Independent recombination data**:
- Rowan et al. 2019 (eLife): high-resolution crossover maps from >2,000 meioses
- Choi et al. 2013 (Nature Genetics): genome-wide recombination maps
- Salome et al. 2012: fine-scale recombination rate variation
- Quadrana et al. 2016 (Nature Genetics): TE insertion polymorphisms across 211 accessions

**Validation approach**:
- Run M/U on all 1,135 accessions -> per-sample S/C distribution
- Group accessions by ecotype, geography, or known recombination modifiers (e.g., FANCM knockouts increase crossovers)
- Test: do accessions from populations/regions with known higher recombination show higher S/C?
- Cost: ~1,135 x 42s = ~13 hours (trivially parallelizable)

**Potential confounders**: TE family composition varies between accessions (Quadrana 2016). S/C could reflect TE insertion history rather than recombination rate per se. Control: compute S/C per-family and check if the pattern holds within families.

### Dataset D: Brachypodium 320 Natural Accessions

| Field | Value |
|-------|-------|
| BioProject | PRJEB73379 |
| Paper | Stritt lab (2024), builds on Stritt et al. 2020 |
| N accessions | 320 |
| Coverage | ~25-30x |
| Genome size | 272 Mb |

**Why valuable**:
- Stritt et al. 2020 (New Phytologist, doi:10.1111/nph.16308) already computed per-family solo:intact ratios using assembly-based methods for 32 families
- Coverage-based TE quantification already done by authors -- direct comparison target
- Known population structure (Mediterranean, Turkish, etc.)
- Active TE families: Angela copia half-life ~66 kyr (fastest documented TE removal in plants)

**Validation**:
- Run M/U on 320 accessions -> per-sample S/C
- Compare per-family S/C against Stritt 2020 published values
- Test: does S/C vary across populations with different demographic histories?

### Dataset E: Rice 3K Genomes -- Direct Recombination Evidence

| Field | Value |
|-------|-------|
| BioProject | PRJEB6180 |
| Paper | 3K Rice Genomes Consortium |
| N accessions | 3,024 (use subset: 100-500) |
| Coverage | ~14x average |
| Genome size | 389 Mb |

**Key advantage**: Du et al. 2012 (Plant Cell) **directly demonstrated** solo LTR frequency increases with local recombination rate in rice. This is the most direct published evidence for the biological premise.

**Validation**:
- Run M/U on 100-500 accessions stratified by subspecies (indica vs japonica -- known to differ in recombination patterns)
- Test: indica vs japonica S/C distribution differences
- Recombination maps: Si et al. 2015, Huang et al. 2009

### Dataset F: Maize HapMap3 -- The Low-Coverage Hard Test

| Field | Value |
|-------|-------|
| BioProject | PRJNA389800 |
| Paper | Bukowski et al. 2018 |
| N accessions | 277 |
| Coverage | ~7x average |
| Genome size | 2.3 Gb |

**Why include**: 7x coverage is realistic for population-scale studies. Combined with NAM genetic maps (McMullen et al. 2009; Rodgers-Melnick et al. 2015), this tests S/C as recombination proxy at realistic low coverage in a challenging genome.

**Challenge**: Maize S/C = 0.15 (very low). Signal-to-noise may be poor. Recent TE amplification dominates. This is explicitly a stress test.

---

## Per-Region Analysis (Future -- Not For Implementation Now)

**Biology note for later consideration**:

The M/U method as currently implemented gives one S/C value per LTR-RT **family** per sample. Since different families insert preferentially in different genomic compartments (e.g., Gypsy in pericentromeric heterochromatin, Copia in gene-rich euchromatin), per-family S/C already carries implicit regional information.

A more explicit per-region approach could:
- Bin genome into windows (e.g., 1 Mb)
- Estimate S/C per window by assigning families to their predominant genomic region
- Correlate with recombination rate per window
- Expected pattern (from literature): centromeric regions show S/C << 1 (suppressed UR), euchromatic arms show S/C >> 1

Key references for regional variation:
- Rice CEN8: solo:intact = 0.9:1 in centromeres vs 2.2:1 in euchromatic arms (Ma & Bennetzen 2006)
- Vitte & Panaud 2003: solo:intact varies from ~0.2 pericentromeric to >2.0 euchromatin in rice
- Cossu 2017 two-regime model: small genomes have high UR/high S/C; large genomes have suppressed UR/high gene conversion

This would require mapping LTR-RT families to genomic coordinates (from reference annotations), which the current read-level pipeline doesn't do. Could be a Phase 2 extension.

---

## LD-Based Orthogonal Validation (Future Consideration)

Tools like **LDhelmet** or **pyrho** estimate population-scaled recombination rate (rho) from LD patterns in SNP data. For any population dataset:
1. Call SNPs (bcftools) -> estimate rho per genomic window
2. Run M/U pipeline -> S/C per sample
3. Correlate population-level rho with population-level S/C

No existing study has combined TE-based and LD-based recombination estimates. This would be the strongest possible validation but is a separate analysis project.

---

## Recommended Execution Priority

| Priority | Tier | Dataset | Rationale |
|----------|------|---------|-----------|
| 1 | Tier 1 | ERR171441 downsampling | Cheapest, uses existing data, answers fundamental coverage question |
| 2 | Tier 2A | Arabidopsis 72 pan-genome | Assembly-level truth, same species as Tier 1 |
| 3 | Tier 3C | Arabidopsis 1001 Genomes | Population-scale test, same species, 1,135 samples |
| 4 | Tier 3D | Brachypodium 320 | Independent species, Stritt ground truth exists |
| 5 | Tier 3E | Rice 3K subset | Du et al. recombination correlation |
| 6 | Tier 2B + 3F | Maize NAM + HapMap3 | Hardest test, largest genome, lowest S/C |

---

## Key Accessions Summary

| Dataset | BioProject | Species | N | Coverage | Ground Truth Type |
|---------|-----------|---------|---|----------|-------------------|
| Pan-genome 72 | PRJEB62038 | A. thaliana | 72 | High | Assembly-level solo:intact |
| 1001 Genomes | PRJNA273563 | A. thaliana | 1,135 | 20-40x | Recombination maps + TE PAV |
| Stritt 320 | PRJEB73379 | B. distachyon | 320 | ~25-30x | Published per-family S/C |
| 3K Rice | PRJEB6180 | O. sativa | 3,024 | ~14x | Du et al. recomb correlation |
| HapMap3 | PRJNA389800 | Z. mays | 277 | ~7x | NAM recombination maps |
| NAM 26 | PRJNA697943 | Z. mays | 26 | Assembly | TE annotations per assembly |

---

## Gaps in Current Knowledge Base

- No notes on LD-based recombination tools (LDhat, LDhelmet, pyrho)
- No formal coverage sensitivity analysis done (T2 5x/10x still TODO)
- No plan connecting S/C to population-level recombination inference
- Arabidopsis 72 pan-genome (PRJEB62038) and Brachypodium 320 (PRJEB73379) not yet in project notes
- Quadrana et al. 2016 TE polymorphism data not captured
