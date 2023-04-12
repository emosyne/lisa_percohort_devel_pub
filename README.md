# Genetic effects of tissue-specific enhancers in schizophrenia and hypertrophic cardiomyopathy
## Doctoral thesis by Emanuele Felice Osimo for Imperial College London
The link to the thesis will be included once publicly deposited.
Contact me at eosimo at ic.ac.uk

### Chapter Chapter 3 - Schizophrenia heritability from partitioned PRSs

This repository contains the code for work contained in Chapter 3. This is a Nextflow pipeline, in which I develop ‘partitioned’ polygenic risk scores, or PRSs where two
genomic compartments (e.g., the tissue-specific enhancers and the residual compartments) are
considered separately for polygenic risk scoring, and then disease heritability is calculated
separately for the original GWAS for the condition (h2
SNP ), as well as for the partitioned
PRSs (h2
pPRS). I also test if multiplying SNP-disease association measure odds ratio for
enhancer-based SNPs by either an effect size of the tissue-specific enhancer, or by its neural
expression, can increase disease h2
pPRS. To do so, this chapter consists of analyses which start from:

- A base GWAS: in this case, either the leave-one-out wave 3 schizophrenia GWAS by
Trubetskoy et al., 2022, or the latest HCM GWAS as by Tadros et al., 2023.
- A target population, either the sample left out in the schizophrenia leave-one-out analysis,
or the UK Biobank.
- A list of genomic coordinates, forming an enhancer-based genomic partition (these
were generated in Chapter 2).

Then, for each base (GWAS), target (population) and enhancer list combination,
the base and target datasets are QCed, and three genomic partitions are created: a 
GWAS partition (called original GWAS), and its two subsets, the tissue-specific enhancer partition
(+100bps at each side), and the residual partition. The three partitions are then clumped
as follows:

- The original GWAS partition is clumped individually, forming the original clumped
GWAS partition.
- The tissue-specific enhancer and the residual partitions are clumped together, after prioritising
enhancer-based SNPs over residual SNPs. The two list are then split back,
creating non-overlapping partitions in terms of LD.

The next step is then computing partitioned PRSs for the condition – schizophrenia
or HCM – for each of the three genomic partitions – using a fixed C+T threshold of 0.5.
For the tissue-specific enhancer partition, in addition to calculating the standard PRS (which
makes use of OR and p-value from the base GWAS), I also calculated modified versions
of the score, utilising OR multipliers derived from AR+C-derived enhancer tissue-specific
expression data. The resulting PRSs are:

- The original clumped GWAS PRS, including all SNPs (after clumping) and using the
GWAS-derived original ORs.
- A residual partition, equivalent to the original base, after subtraction of tissue-specific
enhancer SNPs, and clumping the two lists together – but prioritising enhancer SNPs.
- Three tissue-specific enhancer (TS_ENH) partitions, resulting from all SNPs within 100bps
of the specific enhancer list being considered. This partition is calculated as three distinct
‘versions’: using the original OR for each SNP; using a modified OR enhanced
using the AR+C effect size measure (ES); using a modified OR enhanced using tissuespecific
enhancer expression values.
- Partitioned PRSs (pPRS), as combined logistic models. I compute three pPRSs: one is
a simple additive model, dx ∼ residual + TS_ENH; another one using the additive
model plus interactions, dx ∼ residual × TS_ENH, to account for linear term interactions;
a final one using the additive model plus interactions, plus the quadratic terms:
dx ∼ residual × TS_ENH + residual2 + TS_ENH2.


The output of the analysis is presented as a number of plots, allowing the reader to
compare the adjusted amount of variance explained for the condition – called the coefficient
of determination – for each partition, and for models including more than one partition.
Please see Figure 3.2 for a graphical summary of the methods for this chapter.

Run this pipeline with:

```
nextflow run main.nf -profile servername
```
