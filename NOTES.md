
## Key finding (chr18 dev run)
6 single-mismatch off-targets at 3 distinct loci:
- chr18:22931244 RBBP8 (CtIP, tumour suppressor, HR repair) - hit by ASO_TTR_00575/576/577,
  mismatches at pos 16/17/18 = outside RNase H central window = cleavage competent. REJECT.
- chr18:72771912 NETO1 (brain-expressed) - ASO_TTR_00569, pos 15, terminal. FLAG.
- chr18:5358040 intergenic - ASO_TTR_00106/107, pos 7/8. Lowest concern.
Zero perfect (0mm) off-targets on chr18.

## Day 3 results (chr18 dev run)
- variant_rescore validation: 2282/2282 sites agree with aligner NM tag (0 disagree)
- 377 variant-site pairs scored: 319 protected, 38 neutral, 20 gained
- Top variant-gained: ASO_TTR_00575 chr18:11259602 A>C AF=0.572 (ALT is major allele;
  reference underestimates this site for most of the population), 3mm -> 2mm
- Ranking: 242 scored, 17 excluded as junction-spanning
- Lead candidate: ASO_TTR_00450, zero off-targets on chr18
- Bottom 3: ASO_TTR_00575/576/577, all RBBP8 hitters
- 10 candidates touch RBBP8; scores span 3.1-255.9 (mismatch position drives the spread)
