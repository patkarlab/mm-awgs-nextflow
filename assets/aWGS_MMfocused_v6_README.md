# Adaptive WGS panel v6 — myeloma SV+CNV+SNV focused

## What to load on the P2i

Two files go into MinKNOW's adaptive sampling configuration:

1. **Reference FASTA**: `GCF_009914755.1_T2T-CHM13v2.0_genomic.fa`
   - T2T-CHM13v2.0 (RefSeq/NCBI version)
   - Located on the server at: `/goast/hemat_data/references/T2T/`
   - Uses NC_060925.1, NC_060926.1, ... naming convention

2. **Target BED**: `aWGS_MMfocused_v6_t2t_NC.bed`
   - 38 regions, 24.16 Mb, 0.776% of genome
   - MD5: `ce1819710755e243d7048d84c2ab068a`
   - Uses matching NC_ contig naming

These two files must be loaded together. Do NOT mix the NC_-named BED with a chr-style FASTA or vice versa — MinKNOW won't be able to match the contig names.

## What changed from v5 → v6

The v5 panel (23.94 Mb, 28 regions) was lab-driven additions of MM IGH translocation partners. v6 adds 10 MM driver genes identified by cross-referencing the v5 gene list against four contemporary MM driver-gene catalogs (Walker MGP 2018, Lohr 2014, Maura 2019, Bolli 2025) and trims the TP53 window from ±1 Mb to ±500 kb.

### Two design decisions

**1. TP53 trim: ±1 Mb (v5) → ±500 kb (v6)**

CN-LOH analysis at 17p was the original justification for the v4 ±1 Mb TP53 flank. In practice at adaptive sampling depths (15–20× on-target), the SNP density per kb is insufficient to reliably distinguish CN-LOH from copy-neutral diploidy across the 2.5 Mb v5 window. The diagnostic value lives in the central ~1 Mb (TP53 gene body plus deletion breakpoint resolution for del(17p)); the outer flanks added panel size without analytic payoff. ±500 kb retains del(17p) breakpoint resolution and reclaims ~1.5 Mb of panel budget for higher-yield additions.

**2. Top 10 MM driver genes added** (focal ±50 kb, ~1.7 Mb total)

Cross-referenced against the v5 panel using Walker et al. MGP (15 significantly mutated genes: IRF4, KRAS, NRAS, MAX, HIST1H1E, RB1, EGR1, TP53, TRAF3, FAM46C, DIS3, BRAF, LTB, CYLD, FGFR3), Lohr 2014 (11 significantly mutated genes), Maura 2019 (55–61 driver genes), and Bolli 2025 (62 genes with significantly recurrent point mutations). Top 10 by frequency that were not already in v5:

| Gene | Freq | Pathway | Rationale |
|---|---|---|---|
| DIS3 | ~11% NDMM | RNA exosome | TSG; poor outcome marker (Bolli 2025) |
| TRAF3 | ~5% | NF-κB | Common MM TSG |
| PRDM1 | ~5% | Plasma cell biology | TSG; BLIMP-1 master regulator |
| ATM | ~4–5% | DNA damage response | UK Myeloma XI: combined ATM/ATR signal associated with poor PFS/OS; PARPi-relevant |
| CYLD | ~4% | NF-κB | TSG (Walker MGP top 15) |
| H1-4 (HIST1H1E) | ~4% | Chromatin | Walker MGP top 15 |
| MAX | ~3–4% | MYC pathway | Walker MGP top 15 |
| EGR1 | ~3–4% | Plasma cell biology | Walker MGP top 15; favorable prognosis marker |
| LTB | ~3% | NF-κB | Walker MGP top 15 |
| ATR | ~3% | DNA damage response | Complements ATM; UK Myeloma XI |

All 10 added as focal (±50 kb) regions. Adaptive sampling for these genes targets SNV detection and focal CNV; none has documented MM-specific BCRs that would justify a wider window.

**Genes considered but deferred to v7**: KMT2C, ZFHX4, SP140 (next tier, ~2.5% frequency); KMT2B, IKZF3 (~2%). Adding all 5 would add ~0.6 Mb. Defer until the n is large enough to make the bookkeeping worthwhile.

## v6 regions added or modified

| Region | NC_ coords | chr coords | Size | Provenance | Notes |
|---|---|---|---|---|---|
| ATR | NC_060927.1:145,146,518-145,376,002 | chr3:145,146,518-145,376,002 | 229 kb | v6_new | DDR; ±50 kb flank |
| EGR1 | NC_060929.1:138,941,706-139,045,531 | chr5:138,941,706-139,045,531 | 104 kb | v6_new | Plasma cell TF; ±50 kb flank |
| H1-4 | NC_060930.1:25,974,464-26,075,251 | chr6:25,974,464-26,075,251 | 101 kb | v6_new | Histone; ±50 kb flank |
| LTB | NC_060930.1:31,383,609-31,485,477 | chr6:31,383,609-31,485,477 | 102 kb | v6_new | NF-κB; ±50 kb flank |
| PRDM1 | NC_060930.1:107,118,602-107,335,801 | chr6:107,118,602-107,335,801 | 217 kb | v6_new | BLIMP-1 TSG; ±50 kb flank (RefSeq gene span is 117 kb, includes alternative TSSs) |
| ATM | NC_060935.1:108,180,361-108,426,596 | chr11:108,180,361-108,426,596 | 246 kb | v6_new | DDR; ±50 kb flank |
| DIS3 | NC_060937.1:71,923,642-72,053,383 | chr13:71,923,642-72,053,383 | 130 kb | v6_new | RNA exosome TSG; ±50 kb flank |
| MAX | NC_060938.1:59,160,964-59,357,534 | chr14:59,160,964-59,357,534 | 197 kb | v6_new | MYC pathway; ±50 kb flank |
| TRAF3 | NC_060938.1:96,963,796-97,197,747 | chr14:96,963,796-97,197,747 | 234 kb | v6_new | NF-κB TSG; ±50 kb flank |
| CYLD | NC_060940.1:56,489,726-56,649,616 | chr16:56,489,726-56,649,616 | 160 kb | v6_new | NF-κB TSG; ±50 kb flank |
| TP53 (re-flanked) | NC_060941.1:7,072,543-8,091,594 | chr17:7,072,543-8,091,594 | 1.02 Mb | v5_modified | Was ±1 Mb in v5 (2.52 Mb); now ±500 kb (1.02 Mb). Still covers TNFSF12 (chr17:7.45 Mb) incidentally |

## Panel composition by purpose

| Flank | Genes |
|---|---|
| Focal ±50 kb (SNV, focal CNV, or <1% rare-partner gene-body coverage) | CDKN2C, NRAS, KRAS, BRAF, IRF4, XBP1, RB1, TENT5C, FCRL4/FCRL5, BCL2, CXCR4, LRRK2, MAP3K14, PAX5, TNFAIP8, TXNDC5, **ATR**, **EGR1**, **H1-4**, **LTB**, **PRDM1**, **ATM**, **DIS3**, **MAX**, **TRAF3**, **CYLD** |
| ±250 kb (medium) | FGFR3, NSD2 |
| ±500 kb (translocation BCR) | CCND1, MAFB, CCND3, CCND2, MAFA |
| ±750 kb (wider BCR) | MAF + WWOX (merged) |
| **±500 kb (LOH + breakpoint resolution; also covers TNFSF12 for free)** | **TP53** |
| ±2.5 Mb (MYC) | MYC |
| Locus + 200 kb (clipped to contig end where needed; IGL covers IGLL5 for free) | IGH, IGK, IGL |

## Analytical division of labor

| Variant class | Detection method |
|---|---|
| IGH/IGK/IGL translocations to canonical partners | On-target adaptive sampling + Sniffles2 + CuteSV (+ Severus, nanomonsv) |
| IGH translocations to rare partners (TXNDC5, FCRL4/5, BCL2, etc.) | Same callers + MM-specific annotation dictionary |
| MYC translocations to IGH/IGK/IGL/TXNDC5 | On-target adaptive sampling spanning both sides of the BND |
| SNVs at MM driver genes (now including DIS3, TRAF3, PRDM1, ATM, CYLD, H1-4, MAX, EGR1, LTB, ATR) | On-target adaptive sampling + Clair3 or DeepVariant |
| Focal CNVs and del(17p) breakpoints | On-target adaptive sampling + per-region coverage analysis |
| Large-scale CNVs (1q21 amp, monosomy 13, hyperdiploidy, broad 17p del) | Off-target reads + ichorCNA |

**Note (v6 change)**: CN-LOH at 17p (TP53) has been retired as a v6 analytical target — SNP density at adaptive sampling depths is insufficient. del(17p) breakpoint resolution is preserved in the ±500 kb window.

## Wet-lab parameters (no change from v5 run)

- **DNA shearing**: 21 passes of 26G needle — gave on-target read N50 of 8–11 kb on prior runs, which sits in the published sweet spot for cancer SV adaptive sampling. Keep as-is.
- **Pre-library QC**: Femto Pulse or TapeStation — confirm visible peak at 8–12 kb fragment size.
- **Library prep**: standard SQK-LSK114 protocol (or current equivalent). For multiplexed runs (3-plex recommended), use Native Barcoding Kit 24 V14 (SQK-NBD114.24) plus EXP-NBA114.
- **Adaptive sampling mode**: ENRICH (accept on-target reads, reject off-target).

## Target sequencing depth

Aim for ≥15–20× on-target coverage per sample. With ~24 Mb target and expected enrichment factor of 10–15×, this means roughly 30–40 Gb total sequencing output per sample at 1-plex. At 3-plex multiplexing on a P2 flowcell, expect ~10–12× per sample.

If a sample has lower input DNA mass, prioritize longer sequencing time over higher pore-loading concentration — adaptive sampling rejects more reads from concentrated libraries.

## Per-sample expected output

After sequencing, each sample's BAM file should show on the analysis server:
- Mean on-target coverage: ≥15× target, 5× minimum to call high-confidence SVs
- On-target read N50: 8–12 kb
- Adaptive sampling enrichment factor: 10–15× (on-target reads per Gb sequenced ÷ off-target reads per Gb)

A diagnostic script (`qc_v6_cohort.sh`) is provided to run on each new sample's BAM to confirm these QC metrics against the v6 panel.

## File checksums

```
aWGS_MMfocused_v6_t2t_NC.bed   MD5: ce1819710755e243d7048d84c2ab068a
aWGS_MMfocused_v6_t2t_chr.bed  MD5: b9ce72ca3b1ba1294b8d49bad0b7dab2
```

The `_NC` version is for MinKNOW (matches the NCBI T2T FASTA on the P2i).
The `_chr` version is for downstream analysis on the server (matches the chr1-renamed T2T FASTA).

## Complete genome regions covered (v6)

| Region | NC_ coordinates | chr coordinates | Size | MM relevance |
|---|---|---|---|---|
| CDKN2C | NC_060925.1:50,790,624-50,904,515 | chr1:50,790,624-50,904,515 | 114 kb | del(1p32) — high-risk MM |
| NRAS | NC_060925.1:114,665,929-114,778,216 | chr1:114,665,929-114,778,216 | 112 kb | SNV, 5–10% of MM |
| TENT5C | NC_060925.1:117,565,263-117,687,599 | chr1:117,565,263-117,687,599 | 122 kb | 17% SNV (Bolli 2018); rare IGH partner |
| FCRL4+FCRL5 (IRTA1+IRTA2) | NC_060925.1:156,600,413-156,785,121 | chr1:156,600,413-156,785,121 | 185 kb | 1–2% IGH partner (Mayo) |
| IGK locus | NC_060926.1:88,666,370-90,990,947 | chr2:88,666,370-90,990,947 | 2.32 Mb | Light chain translocation anchor |
| CXCR4 | NC_060926.1:136,508,831-136,612,630 | chr2:136,508,831-136,612,630 | 104 kb | <1% IGH partner (focal) |
| **ATR** | **NC_060927.1:145,146,518-145,376,002** | **chr3:145,146,518-145,376,002** | **229 kb** | **~3% MM SNV; DDR; UK Myeloma XI prognostic** |
| FGFR3 / NSD2 | NC_060928.1:1,541,772-2,230,701 | chr4:1,541,772-2,230,701 | 689 kb | t(4;14) partner |
| TNFAIP8 | NC_060929.1:119,739,607-119,970,525 | chr5:119,739,607-119,970,525 | 231 kb | <1% IGH partner |
| **EGR1** | **NC_060929.1:138,941,706-139,045,531** | **chr5:138,941,706-139,045,531** | **104 kb** | **~3–4% MM SNV; Walker MGP top 15; favorable prognosis** |
| IRF4 | NC_060930.1:200,136-319,771 | chr6:200,136-319,771 | 120 kb | Rare t(IRF4;14); plasma cell TF; lenalidomide target |
| TXNDC5 | NC_060930.1:7,700,266-7,829,537 | chr6:7,700,266-7,829,537 | 129 kb | 4.9% MM partner (CoMMpass); 3rd most common MYC partner |
| **H1-4** | **NC_060930.1:25,974,464-26,075,251** | **chr6:25,974,464-26,075,251** | **101 kb** | **~4% MM SNV; chromatin; Walker MGP top 15** |
| **LTB** | **NC_060930.1:31,383,609-31,485,477** | **chr6:31,383,609-31,485,477** | **102 kb** | **~3% MM SNV; NF-κB; Walker MGP top 15** |
| CCND3 | NC_060930.1:41,263,493-42,378,294 | chr6:41,263,493-42,378,294 | 1.11 Mb | t(6;14) partner |
| **PRDM1** | **NC_060930.1:107,118,602-107,335,801** | **chr6:107,118,602-107,335,801** | **217 kb** | **~5% MM SNV; BLIMP-1 master regulator** |
| BRAF | NC_060931.1:141,977,505-142,289,131 | chr7:141,977,505-142,289,131 | 312 kb | SNV (V600E) — MEK inhibitor target |
| MYC | NC_060932.1:126,362,888-131,370,405 | chr8:126,362,888-131,370,405 | 5.01 Mb | MYC translocation BCR |
| MAFA | NC_060932.1:144,077,736-145,080,401 | chr8:144,077,736-145,080,401 | 1.00 Mb | <1% Ig translocation partner |
| PAX5 | NC_060933.1:36,806,835-37,107,984 | chr9:36,806,835-37,107,984 | 301 kb | <1% IGH partner |
| CCND1 | NC_060935.1:69,158,031-70,171,351 | chr11:69,158,031-70,171,351 | 1.01 Mb | t(11;14) partner |
| **ATM** | **NC_060935.1:108,180,361-108,426,596** | **chr11:108,180,361-108,426,596** | **246 kb** | **~4–5% MM SNV; DDR; PARPi-relevant** |
| CCND2 | NC_060936.1:3,780,521-4,812,135 | chr12:3,780,521-4,812,135 | 1.03 Mb | <1% Ig translocation partner |
| KRAS | NC_060936.1:25,026,496-25,172,152 | chr12:25,026,496-25,172,152 | 146 kb | SNV, 15–25% of MM |
| LRRK2 | NC_060936.1:40,127,355-40,371,422 | chr12:40,127,355-40,371,422 | 244 kb | <1% IGH partner |
| RB1 | NC_060937.1:47,474,085-47,752,182 | chr13:47,474,085-47,752,182 | 278 kb | del(13q)/monosomy 13 marker |
| **DIS3** | **NC_060937.1:71,923,642-72,053,383** | **chr13:71,923,642-72,053,383** | **130 kb** | **~11% MM SNV; RNA exosome TSG; poor outcome (Bolli 2025)** |
| **MAX** | **NC_060938.1:59,160,964-59,357,534** | **chr14:59,160,964-59,357,534** | **197 kb** | **~3–4% MM SNV; MYC pathway** |
| **TRAF3** | **NC_060938.1:96,963,796-97,197,747** | **chr14:96,963,796-97,197,747** | **234 kb** | **~5% MM SNV; NF-κB TSG** |
| IGH locus | NC_060938.1:99,639,469-101,161,492 | chr14:99,639,469-101,161,492 | 1.52 Mb | Primary Ig translocation anchor (3' clipped to chr14 telomere) |
| **CYLD** | **NC_060940.1:56,489,726-56,649,616** | **chr16:56,489,726-56,649,616** | **160 kb** | **~4% MM SNV; NF-κB TSG** |
| WWOX + MAF (merged) | NC_060940.1:83,955,399-86,407,071 | chr16:83,955,399-86,407,071 | 2.45 Mb | t(14;16) breakpoint cluster (WWOX intron 8) + MAF gene |
| **TP53 (re-flanked)** | **NC_060941.1:7,072,543-8,091,594** | **chr17:7,072,543-8,091,594** | **1.02 Mb** | **del(17p) breakpoint resolution; ±500 kb flank (was ±1 Mb in v5)** |
| MAP3K14 (NIK) | NC_060941.1:46,066,902-46,220,784 | chr17:46,066,902-46,220,784 | 154 kb | <1% IGH partner; NF-κB pathway |
| BCL2 (focal) | NC_060942.1:63,276,497-63,575,122 | chr18:63,276,497-63,575,122 | 299 kb | <1% IGH partner; venetoclax-eligibility marker |
| MAFB | NC_060944.1:41,917,334-42,920,722 | chr20:41,917,334-42,920,722 | 1.00 Mb | t(14;20) partner |
| IGL locus (also covers IGLL5) | NC_060946.1:22,239,629-23,545,823 | chr22:22,239,629-23,545,823 | 1.31 Mb | Light chain translocation anchor; incidentally covers IGLL5 (chr22:23.31 Mb, 18% SNV frequency) |
| XBP1 | NC_060946.1:29,205,848-29,311,869 | chr22:29,205,848-29,311,869 | 106 kb | Plasma cell master TF; non-Ig MYC hijacking partner |

Total: 24.16 Mb across 38 regions (0.776% of T2T-CHM13v2.0 genome).

## v3 → v4 → v5 → v6 evolution summary

| Version | Size | Regions | % genome | Key change |
|---|---|---|---|---|
| v3 | 24.24 Mb | 23 | 0.782% | T2T migration + BCR expansion; recovered all 3 FISH+ translocations from prior run |
| v4 | 22.17 Mb | 19 | 0.711% | Removed lymphoma/Waldenström regions; added WWOX/IRF4/XBP1; expanded TP53 ±1 Mb |
| v5 | 23.94 Mb | 28 | 0.768% | Lab-driven additions of 9 documented IGH partners (Tier 1+2+3); Tier 4 novel partners excluded |
| **v6** | **24.16 Mb** | **38** | **0.776%** | **Added 10 top-mutated MM driver genes (DIS3, TRAF3, PRDM1, ATM, CYLD, H1-4, MAX, EGR1, LTB, ATR); trimmed TP53 from ±1 Mb to ±500 kb; retired CN-LOH at 17p as an analytical target** |

## Bold = added or modified in v6

## Build provenance

This panel was built deterministically from v5 + the T2T-CHM13v2.0 NCBI RefSeq GFF using `build_v6_panel.py`. The 27 v5-retained regions are byte-identical to the corresponding v5 BED rows; the 11 modified/new regions (TP53 + 10 new genes) come from RefSeq gene-span coordinates plus the per-gene flanks documented above.

Build artifacts on the server:
- Build script: `~/inbox/from_claude/build_v6_panel.py` (or wherever the user filed it)
- Build directory: `/goast/nikhil_awgs_testing/panel/v6_build/`
  - `aWGS_MMfocused_v6_BUILD_STAMP.txt` — inputs, MD5s, warnings
  - `aWGS_MMfocused_v6_build_report.tsv` — per-region provenance with raw gene-body coordinates
