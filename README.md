Variant-Aware Off-Target Screening for Antisense Oligonucleotides

A pipeline that screens candidate antisense oligonucleotides (ASOs) against the human genome and then asks a question most screens skip: does this off-target profile change in people who carry common genetic variants?

Target gene: TTR (transthyretin), the target of the approved ASO inotersen.

The problem

An ASO is a ~20 nt oligonucleotide designed to bind one transcript and trigger RNase H1-mediated cleavage. Twenty bases is not a lot of information: near-matches occur throughout a 3 Gb genome, and binding at an unintended site silences a gene the drug was never meant to touch.

Standard screens align candidates against the reference genome. But a reference genome is one sample. A site carrying 2 mismatches in GRCh38 may carry only 1 in someone homozygous for a common variant — making an ASO acceptably specific for most people and not for them. This pipeline quantifies that.

Pipeline
TTR transcript (GENCODE)
        │
        ▼  tile_candidates.py      20-mers, 1 nt step, design filters
   candidates.fa
        │
        ▼  bowtie2 -a              end-to-end, gapless, ≤3 mismatches
   hits.bam
        │
        ▼  parse_alignments.py     mismatch positions in ASO coordinates
   hits_table.tsv ──────────────▶ offtarget_sites.bed
        │                                  │
        │                                  ▼  bcftools -R
        │                          gnomAD v4.1 SNVs (PASS, AF>0.001)
        │                                  │
        ▼  variant_rescore.py ◀────────────┘
   variant_rescored.tsv            REF vs ALT duplex quality
        │
        ▼  score_and_rank.py       composite risk, gene context
   candidates_ranked.tsv
        │
        ▼  make_figures.py
   figures/
Key design decisions

Genome alignment, not transcriptome. Gapmer ASOs act on nuclear pre-mRNA, so intronic off-targets are real. Aligning to GRCh38 also makes variant intersection a coordinate join rather than a transcript→genome conversion with its attendant off-by-one errors.

Mismatch position, not just count. RNase H1 requires a well-formed duplex across its footprint, which sits over the oligo's central region. A mismatch at position 1 or 20 leaves cleavage largely intact; one at position 10 can abolish it. A 2-mismatch off-target with both mismatches terminal is more dangerous than a 1-mismatch off-target with the mismatch central. Scoring tracks both.

Strand convention made explicit. The ASO is the reverse complement of its target site. Minus-strand alignments store the reverse complement in SAM SEQ, so mismatch positions must be flipped back to ASO coordinates (aso_pos = len - 1 - stored_pos). Getting this wrong mirrors every position around the oligo centre — invisible in aggregate statistics, wrong in every detail. Unit-tested in scripts/test_mismatch_positions.py.

Bowtie2 forced into exhaustive short-read mode. -N 1 -L 10 -i C,1,0 for sensitive seeding, --mp 6,6 with --score-min L,-18,0 to admit exactly ≤3 mismatches, and --rdg/--rfg 100,100 to forbid gaps.

Validation

variant_rescore.py independently recomputes every alignment's mismatch count from the reference FASTA and compares it against the aligner's NM tag.

Validation: 2282 sites agree with NM tag, 0 disagree

Agreement across all 2,282 alignments confirms the genome build, coordinate handling, strand convention, and reverse-complement logic are correct.

Two further internal checks fall out of the data rather than being enforced:

Adjacent tiles (overlapping by 19/20 nt) hit the same locus with the mismatch position sliding by exactly one per tile.
Where a variant removes a mismatch, the variant's ASO position equals the removed mismatch position — computed independently.
Results (chr18 development run)

259 candidates from 597 tiled windows (43% pass rate; CpG content was the leading rejection reason, as expected for GC-rich TTR).

Metric	Value
Alignments	2,282
Perfect (0 mm) off-targets	0
1-mismatch off-targets	6, at 3 distinct loci
Variant–site pairs scored	377
Variant-gained sites	20
Variant-protected sites	319

Off-target gene context. The six 1-mismatch off-targets resolve to three loci:

chr18:22,931,244 — RBBP8 (CtIP), a tumour suppressor required for homologous recombination. Hit by candidates 575/576/577 with mismatches at ASO positions 16/17/18 — all outside the RNase H footprint, so cleavage-competent. Rejected.
chr18:72,771,912 — NETO1, brain-expressed. Mismatch at position 15. Flagged.
chr18:5,358,040 — intergenic. Lowest concern.

Variant-driven risk. The strongest case is ASO_TTR_00575 at chr18:11,259,602 (A>C, AF 0.572). The ALT is the major allele: GRCh38 carries the minor allele here, so this site is 2-mismatch for most of the population while a reference-only screen reports 3. Reference screening systematically underestimates it.

Variants destroy chance similarity ~16× more often than they create it (319 protected vs 20 gained) — most single-base changes make a near-match worse.

Ranking. 242 candidates scored, 17 excluded as junction-spanning. Lead candidate ASO_TTR_00450 has zero off-targets on chr18. The bottom three ranks are 575/576/577, recovered independently by the composite score. Ten candidates touch RBBP8 with scores spanning 3.1–255.9 — the same gene at radically different risk, separated by mismatch position.

Limitations

Stated plainly, because they bound what the results mean:

SNVs only. Indels shift local coordinates and require haplotype reconstruction. 17 of 252 filtered variants (7%) are indels; they are counted and reported, not scored.
Splice junctions invisible. 17 candidates span exon-exon boundaries and have no contiguous genomic match. End-to-end genomic alignment cannot evaluate them. This matters: nusinersen works by binding across a splice site. Transcriptome-level alignment would be needed.
Single variants, not haplotypes. Each variant is scored independently. Co-occurring variants in linkage could have combined effects.
No expression or chromatin filter. An off-target in a gene not expressed in the dosed tissue is less concerning than this scoring implies.
Three-mismatch hits dominate by count and are largely noise. 1,915 of 2,029 off-targets are 3-mismatch; most would not cleave efficiently. They are down-weighted, not excluded.
AF > 0.001 threshold excludes rare variants that matter to individuals even if not to populations.
Thermodynamics not modelled. Scoring uses mismatch count and position as a proxy for duplex stability rather than computing ΔG.
Reproducing
bash
# Environment (versions pinned: numpy 2.x breaks pandas built against 1.x)
conda create -p ./envs/aso -c conda-forge -c bioconda \
  python=3.11 "numpy<2" "pandas<2.2" matplotlib seaborn pysam biopython bowtie

# Cluster tools from modules
module load bowtie2 samtools bedtools bcftools

References: GRCh38 (no-alt), GENCODE v44, gnomAD v4.1 genomes.

Contig naming differs between sources — GENCODE and UCSC use chr18, ClinVar uses 18. Unify with bcftools annotate --rename-chrs before any intersection. A silent zero-overlap looks like "no variants found" rather than an error.

See workflow/ for SLURM submission scripts.

Repository layout
scripts/
  tile_candidates.py           candidate generation + design filters
  parse_alignments.py          BAM → tidy hit table, ASO-relative positions
  test_mismatch_positions.py   unit test for the strand flip
  variant_rescore.py           REF vs ALT duplex quality
  score_and_rank.py            composite risk score
  make_figures.py              figures 1–5
workflow/
  align_production.sh          genome-wide alignment (SLURM)
  get_gnomad.sh                reference download (SLURM)
config/
  config.yaml                  all paths and parameters
results/                       outputs (not tracked)
