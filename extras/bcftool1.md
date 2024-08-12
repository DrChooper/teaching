## bcftool command explanation
The command `bcftools stats -s - Av.vcf | grep "PSC"` is used to generate and filter variant statistics from a VCF (Variant Call Format) file (`Av.vcf`). Let's break down the components of this command:

1. **`bcftools stats -s - Av.vcf`**:
   - **`bcftools stats`**: This is a command from the `bcftools` suite used to generate statistics about variants in a VCF file. It calculates various statistics, including the number of reference homozygous (nRefHom), non-reference homozygous (nNonRefHom), heterozygous (nHets), and other metrics.
   - **`-s -`**: The `-s -` option tells `bcftools` to generate statistics for all samples in the VCF file. The `-` is used as a placeholder indicating that no specific samples are being selected, so all samples will be included.
   - **`Av.vcf`**: This is the input VCF file containing variant data for which the statistics are being generated.

2. **`| grep "PSC"`**:
   - **`|`**: This is a pipe operator that takes the output of the `bcftools stats` command and passes it as input to the next command, which is `grep`.
   - **`grep "PSC"`**: `grep` is a command-line utility used to search for specific patterns in text. In this case, it filters the output of `bcftools stats` to show only the lines that contain the string "PSC".

3. **`PSC` Line Explanation**:
   - The `PSC` line in the output stands for **Per-Sample Counts**. 
   
   
#### more to output sections
The output of `bcftools stats` is organized into multiple sections, each identified by a specific code (e.g., "SN," "TSTV," "PSC"). Here's a brief overview:

1. **SN - Summary Numbers**: Provides a total count of different variant types (SNPs, indels, etc.) across the entire dataset.

2. **TSTV - Transition/Transversion Statistics**: Reports the ratio of transitions (purine-purine or pyrimidine-pyrimidine changes) to transversions (purine-pyrimidine changes), a key quality metric.

3. **ST - Substitution Types**: Categorizes specific base substitutions, helping to understand mutation patterns.

4. **PSC - Per-Sample Counts**: Offers detailed statistics for each sample, including counts of homozygous and heterozygous sites, transitions, transversions, indels, and coverage depth.

5. **AF - Allele Frequency Spectrum**: Analyzes allele frequencies across the population, indicating the prevalence of different variants.

6. **IDD - Indel Distribution**: Focuses on the lengths and distribution of indels, providing insights into the indel spectrum.

7. **QUAL - Quality Score Distribution**: Displays the distribution of variant quality scores, reflecting the confidence in the calls.

8. **DP - Depth Distribution**: Reports the distribution of sequencing depth across variants, important for coverage assessment.

Each section provides targeted insights into different aspects of the variant data, making it easier to assess the quality and characteristics of the dataset.
   
---    
### PSC Output Explanation

This line includes counts related to the variant calls for each sample in the VCF file. Specifically, it reports:

- **`[2]id`**: An identifier for the sample (usually 0 or another number).
- **`[3]sample`**: The sample name or path to the BAM file associated with the sample.
- **`[4]nRefHom`**: The number of sites where the sample is homozygous for the reference allele (no variants).
- **`[5]nNonRefHom`**: The number of sites where the sample is homozygous for a non-reference allele (variants present).
- **`[6]nHets`**: The number of sites where the sample is heterozygous (one reference allele and one variant allele).
- **`[7]nTransitions`**: The number of transition mutations (e.g., A↔G or C↔T).
- **`[8]nTransversions`**: The number of transversion mutations (e.g., A↔C, A↔T, G↔C, G↔T).
- **`[9]nIndels`**: The number of insertions or deletions (indels).
- **`[10]average depth`**: The average depth of coverage for the sample.
- **`[11]nSingletons`**: The number of variants that are unique to this sample (singletons).
- **`[12]nHapRef`**: The number of sites where the sample has a haploid reference allele.
- **`[13]nHapAlt`**: The number of sites where the sample has a haploid alternate allele.
- **`[14]nMissing`**: The number of sites where the data is missing for this sample.

## Explanation of the Output and Contamination Concerns

When analyzing chloroplast genomes, which are haploid and present in a single copy within an organism, the expectation when aligning sequences from the same species is to observe little to no variation. Here’s a detailed explanation of the output and why it might raise concerns about potential contamination:

### Expected Findings in a Chloroplast Genome from the Same Species

- **Haploid Nature**: The chloroplast genome is haploid, so there should be only one allele at each position—matching the reference genome (nRefHom).
- **Absence of Variation**: Since all samples come from the same species, you would expect very few, if any, variants. Ideally, the counts for nNonRefHom (non-reference homozygous) and nHets (heterozygous sites) should be zero.
- **Consistency Across Samples**: All samples should display similar or identical variant statistics if they are from the same species and there is no contamination.

### Observations in the Output

| PSC | [2]id | [3]sample                                                        | [4]nRefHom | [5]nNonRefHom | [6]nHets | [7]nTransitions | [8]nTransversions | [9]nIndels | [10]average depth | [11]nSingletons | [12]nHapRef | [13]nHapAlt | [14]nMissing |
|-----|-------|------------------------------------------------------------------|------------|---------------|----------|-----------------|-------------------|------------|-------------------|-----------------|-------------|-------------|--------------|
| PSC | 0     | /mnt/s-ws/everyone/SCIEM401/Module_6_Variants/bams/1.sorted.bam  | 0          | 0             | 0        | 0               | 0                 | 0          | 0.0               | 0               | 355         | 79          | 0            |
| PSC | 0     | /mnt/s-ws/everyone/SCIEM401/Module_6_Variants/bams/10.sorted.bam | 0          | 0             | 0        | 0               | 0                 | 0          | 0.0               | 0               | 372         | 62          | 0            |
| PSC | 0     | /mnt/s-ws/everyone/SCIEM401/Module_6_Variants/bams/11.sorted.bam | 0          | 0             | 0        | 0               | 0                 | 0          | 0.0               | 0               | 364         | 70          | 0            |
| PSC | 0     | /mnt/s-ws/everyone/SCIEM401/Module_6_Variants/bams/12.sorted.bam | 0          | 0             | 0        | 0               | 0                 | 0          | 0.0               | 0               | 344         | 89          | 1            |
| PSC | 0     | /mnt/s-ws/everyone/SCIEM401/Module_6_Variants/bams/13.sorted.bam | 0          | 0             | 0        | 0               | 0                 | 0          | 0.0               | 0               | 338         | 96          | 0            |
| PSC | 0     | /mnt/s-ws/everyone/SCIEM401/Module_6_Variants/bams/14.sorted.bam | 0          | 0             | 0        | 0               | 0                 | 0          | 0.0               | 0               | 319         | 115         | 0            |
| PSC | 0     | /mnt/s-ws/everyone/SCIEM401/Module_6_Variants/bams/15.sorted.bam | 0          | 0             | 0        | 0               | 0                 | 0          | 0.0               | 0               | 300         | 134         | 0            |
| PSC | 0     | /mnt/s-ws/everyone/SCIEM401/Module_6_Variants/bams/16.sorted.bam | 0          | 0             | 0        | 0               | 0                 | 0          | 0.0               | 0               | 163         | 241         | 30           |
| PSC | 0     | /mnt/s-ws/everyone/SCIEM401/Module_6_Variants/bams/2.sorted.bam  | 0          | 0             | 0        | 0               | 0                 | 0          | 0.0               | 0               | 347         | 87          | 0            |
| PSC | 0     | /mnt/s-ws/everyone/SCIEM401/Module_6_Variants/bams/3.sorted.bam  | 0          | 0             | 0        | 0               | 0                 | 0          | 0.0               | 0               | 303         | 131         | 0            |
| PSC | 0     | /mnt/s-ws/everyone/SCIEM401/Module_6_Variants/bams/4.sorted.bam  | 0          | 0             | 0        | 0               | 0                 | 0          | 0.0               | 0               | 279         | 155         | 0            |
| PSC | 0     | /mnt/s-ws/everyone/SCIEM401/Module_6_Variants/bams/5.sorted.bam  | 0          | 0             | 0        | 0               | 0                 | 0          | 0.0               | 0               | 366         | 68          | 0            |
| PSC | 0     | /mnt/s-ws/everyone/SCIEM401/Module_6_Variants/bams/6.sorted.bam  | 0          | 0             | 0        | 0               | 0                 | 0          | 0.0               | 0               | 307         | 125         | 2            |
| PSC | 0     | /mnt/s-ws/everyone/SCIEM401/Module_6_Variants/bams/7.sorted.bam  | 0          | 0             | 0        | 0               | 0                 | 0          | 0.0               | 0               | 333         | 101         | 0            |
| PSC | 0     | /mnt/s-ws/everyone/SCIEM401/Module_6_Variants/bams/8.sorted.bam  | 0          | 0             | 0        | 0               | 0                 | 0          | 0.0               | 0               | 402         | 32          | 0            |
| PSC | 0     | /mnt/s-ws/everyone/SCIEM401/Module_6_Variants/bams/9.sorted.bam  | 0          | 0             | 0        | 0               | 0                 | 0          | 0.0               | 0               | 306         | 125         | 3            |

### Analysis of the Results
1. **Zero Counts in nRefHom and nNonRefHom**: All samples show zero counts for both nRefHom (reference homozygous) and nNonRefHom (non-reference homozygous) - indication that the sample matches the reference genome.

2. **Zero Counts in nHets**: As expected for a haploid genome, there are zero heterozygous sites across all samples.

3. **Non-Zero Counts in Singleton, HapRef, and HapAlt**:
   - **Singletons**: The presence of singleton variants (variants that are unique to one sample) is unexpected in samples from the same species and may indicate sequencing errors or contamination.
   - **nHapRef and nHapAlt**: These counts reflect the presence of reference and alternate alleles, which is normal. However, the discrepancies in singleton and missing data counts suggest potential issues with the data.

4. **Missing Data**: The non-zero counts in the nMissing column for some samples indicate that data is missing for certain sites. This could be due to low coverage or poor-quality sequencing and could also be a sign of contamination.

### Conclusion

Given that these samples are all from the same species, the output should have shown consistent and similar values across all samples with most sites categorized as nRefHom. The observed zero counts for nRefHom and nNonRefHom across all samples, along with the presence of singleton variants and missing data, suggest that there may be contamination or sequencing artifacts in the dataset. This could result from cross-sample contamination or technical issues during sequencing or data processing.

## Cleaned output
After removing contaminating sequences, the cleaned BAM files were analyzed to compare the variant calls. The table below shows the per-sample counts (`PSC`) from the cleaned data.

### Cleaned Output Table

| PSC | [2]id | [3]sample                                                        | [4]nRefHom | [5]nNonRefHom | [6]nHets | [7]nTransitions | [8]nTransversions | [9]nIndels | [10]average depth | [11]nSingletons | [12]nHapRef | [13]nHapAlt | [14]nMissing |
|-----|-------|------------------------------------------------------------------|------------|---------------|----------|-----------------|-------------------|------------|-------------------|-----------------|-------------|-------------|--------------|
| PSC | 0     | /mnt/s-ws/everyone/SCIEM401/Module_6_Variants/clean_bams/1.sorted.bam  | 0          | 0             | 0        | 0               | 0                 | 0          | 00.0              | 0               | 120         | 152         | 0            |
| PSC | 0     | /mnt/s-ws/everyone/SCIEM401/Module_6_Variants/clean_bams/10.sorted.bam | 0          | 0             | 0        | 0               | 0                 | 0          | 00.0              | 0               | 126         | 146         | 0            |
| PSC | 0     | /mnt/s-ws/everyone/SCIEM401/Module_6_Variants/clean_bams/11.sorted.bam | 0          | 0             | 0        | 0               | 0                 | 0          | 00.0              | 0               | 141         | 131         | 0            |
| PSC | 0     | /mnt/s-ws/everyone/SCIEM401/Module_6_Variants/clean_bams/12.sorted.bam | 0          | 0             | 0        | 0               | 0                 | 0          | 00.0              | 0               | 138         | 133         | 1            |
| PSC | 0     | /mnt/s-ws/everyone/SCIEM401/Module_6_Variants/clean_bams/13.sorted.bam | 0          | 0             | 0        | 0               | 0                 | 0          | 00.0              | 0               | 97          | 175         | 0            |
| PSC | 0     | /mnt/s-ws/everyone/SCIEM401/Module_6_Variants/clean_bams/14.sorted.bam | 0          | 0             | 0        | 0               | 0                 | 0          | 00.0              | 0               | 123         | 149         | 0            |
| PSC | 0     | /mnt/s-ws/everyone/SCIEM401/Module_6_Variants/clean_bams/15.sorted.bam | 0          | 0             | 0        | 0               | 0                 | 0          | 00.0              | 0               | 88          | 184         | 0            |
| PSC | 0     | /mnt/s-ws/everyone/SCIEM401/Module_6_Variants/clean_bams/16.sorted.bam | 0          | 0             | 0        | 0               | 0                 | 0          | 00.0              | 0               | 119         | 143         | 10           |
| PSC | 0     | /mnt/s-ws/everyone/SCIEM401/Module_6_Variants/clean_bams/2.sorted.bam  | 0          | 0             | 0        | 0               | 0                 | 0          | 00.0              | 0               | 121         | 151         | 0            |
| PSC | 0     | /mnt/s-ws/everyone/SCIEM401/Module_6_Variants/clean_bams/3.sorted.bam  | 0          | 0             | 0        | 0               | 0                 | 0          | 00.0              | 0               | 78          | 194         | 0            |
| PSC | 0     | /mnt/s-ws/everyone/SCIEM401/Module_6_Variants/clean_bams/4.sorted.bam  | 0          | 0             | 0        | 0               | 0                 | 0          | 00.0              | 0               | 77          | 195         | 0            |
| PSC | 0     | /mnt/s-ws/everyone/SCIEM401/Module_6_Variants/clean_bams/5.sorted.bam  | 0          | 0             | 0        | 0               | 0                 | 0          | 00.0              | 0               | 131         | 141         | 0            |
| PSC | 0     | /mnt/s-ws/everyone/SCIEM401/Module_6_Variants/clean_bams/6.sorted.bam  | 0          | 0             | 0        | 0               | 0                 | 0          | 00.0              | 0               | 82          | 190         | 0            |
| PSC | 0     | /mnt/s-ws/everyone/SCIEM401/Module_6_Variants/clean_bams/7.sorted.bam  | 0          | 0             | 0        | 0               | 0                 | 0          | 00.0              | 0               | 115         | 157         | 0            |
| PSC | 0     | /mnt/s-ws/everyone/SCIEM401/Module_6_Variants/clean_bams/8.sorted.bam  | 0          | 0             | 0        | 0               | 0                 | 0          | 00.0              | 0               | 135         | 137         | 0            |
| PSC | 0     | /mnt/s-ws/everyone/SCIEM401/Module_6_Variants/clean_bams/9.sorted.bam  | 0          | 0             | 0        | 0               | 0                 | 0          | 00.0              | 0               | 87          | 184         | 1            |

#### Analysis

The comparison of the cleaned data shows the following:

- **nRefHom** and **nNonRefHom**: Both columns remain at 0, as expected, given the nature of the data and the clean-up process.
- **nHapRef** and **nHapAlt**: These columns reflect the number of haploid reference and alternate alleles, with variation observed between samples.
- **nMissing**: Only a few samples have missing data, indicating that most data was successfully recovered after cleaning.

This cleaned data output suggests that the removal of contaminants has led to a more consistent and reliable set of variant calls across the samples.



