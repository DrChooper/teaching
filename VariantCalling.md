
# Variant Calling

## Content
# Table of Contents

- [Mapping Reads to Reference using bbmap](#mapping-reads-to-reference-using-bbmap)
- [Variant Calling with Original Data](#variant-calling-with-original-data)
- [Quality Control and Clean up](#quality-control-and-clean-up)
- [Variant Calling with Cleaned Data](#variant-calling-with-cleaned-data)
- [Variant Querying and Annotation](#variant-querying-and-annotation)

## Introduction
This workflow details the analysis of chloroplast DNA sequencing data, encompassing read mapping, variant calling, quality control, contamination removal, and annotation to identify genetic variations across the genome.

### Mapping Reads to reference using bbmap

**Objective:**
- To align sequence reads to a reference genome and generate a Sequence Alignment/Map (SAM) file.

**Background:**
- The reads have already been trimmed to remove adapter sequences.
- We are using `bbmap` as the alignment software. More information about `bbmap` can be found [here](https://sourceforge.net/projects/bbmap/).

**Step-by-Step Instructions:**

1. **Set the Source Directory:**
    - Set the source directory to the provided data (sequence reads)

    ```bash
    sourcedir=/mnt/s-ws/everyone/SCIEM401/Module_6_Variants
    ```
    - Ensure you are in your folder (check with `pwd`) which is something like `/mnt/s-ws/s-99`.
    - make sure your folder contains the `Av.cp.final.fasta` file (reference genome)
    
2. **Map Reads to the Reference Genome:**
    - Use `bbmap` to align the reads and generate BAM files.

    ```bash
    for sample in {1..2}
    do
    bbmap.sh in1=$sourcedir/$sample.cp.R1.trimmed.fq.gz in2=$sourcedir/$sample.cp.R2.trimmed.fq.gz out=$sample.bam ref=Av.cp.final.fasta mappedonly=t
    done
    ```

3. **Sort and Index BAM Files:**
    - Use `samtools` to sort and index the BAM files. This step improves the efficiency of downstream analysis.

    ```bash
    for sample in {1..2}
    do
    samtools sort -T tmp -o $sample.sorted.bam $sample.bam
    samtools index $sample.sorted.bam
    done
    ```

4. **Comparing Sizes of Sorted and Unsorted BAM Files**
    - Use the `ls` command to list the files in the directory and check their sizes.
    
    ```bash
    ls -lh
    ```

    - Compare the sizes of the sorted and unsorted BAM files.
    - Are the sorted files generally smaller or larger?
    *Sorted BAM files are typically smaller due to better compression of the sorted data.*

5. **Remove Unsorted BAM Files:**
    - Use the `rm` command to remove the unsorted BAM files.

    ```bash
    rm [0-9].bam
    ```

### Variant Calling with Original Data
- We use `bcftools` to call variants. This process will take several minutes.

- The output file is in VCF format, which is described in detail [here](http://samtools.github.io/hts-specs/). For large files, a binary (compressed) format called BCF can be used, but for this module, we'll stick to VCF format as the files are small.

- By default, `bcftools` assumes the samples are from a diploid organism and thus may be heterozygous for variants. However, since our samples are of chloroplast DNA, which like bacterial DNA, can be safely considered uniform, we treat these samples as effectively haploid using the `--ploidy 1` option.

```bash
bcftools mpileup -Ou -f Av.cp.final.fasta $sourcedir/bams/*.sorted.bam | bcftools call --ploidy 1 -mv -Ov -o Av.vcf
```

- To obtain a summary of the results, we use `bcftools stats`. By passing the output through `grep`, we can view the per-sample counts:

```bash
bcftools stats -s - Av.vcf | grep "PSC"
```

- The summary output will show columns such as `nHapRef` (same as reference) and `nHapAlt` (different from reference) since we specified the samples as haploid. If the samples were diploid, the variants would be counted in the `nRefHom` (homozygous reference), `nNonRefHom` (homozygous variant), and `nHets` (heterozygous) columns.

Here is an example of what the output might look like:

```bash
# PSC, Per-sample counts. Note that the ref/het/hom counts include only SNPs, for indels see PSI. The rest include both SNPs and indels.
# PSC   [2]id   [3]sample                                                       [4]nRefHom      [5]nNonRefHom   [6]nHets        [7]nTransitions  [8]nTransversions        [9]nIndels      [10]average depth       [11]nSingletons [12]nHapRef    [13]nHapAlt      [14]nMissing
PSC     0       /mnt/s-ws/everyone/SCIEM401/Module_6_Variants/bams/1.sorted.bam 0             0              0              0                0                0             0.0                 0              355          79            0
PSC     0       /mnt/s-ws/everyone/SCIEM401/Module_6_Variants/bams/10.sorted.bam 0             0              0              0                0                0             0.0                 0              372          62            0
PSC     0       /mnt/s-ws/everyone/SCIEM401/Module_6_Variants/bams/11.sorted.bam 0             0              0              0                0                0             0.0                 0              364          70            0
...

```

- Analyze the results to answer the following questions:
  - Which sample shows the least variation from the reference?
  - Which sample shows the most variation?
  - Note that one sample shows a surprisingly large number of variations. This might indicate contamination of the sample with sequencing reads from another plant.

### Quality Control and Clean up

**For example: Removal of Arabidopsis Contamination**

In this step, we will utilize `bbmap` to eliminate reads that map to the Arabidopsis thaliana chloroplast genome. Arabidopsis is a common contaminant in sequencing experiments due to its widespread use in research labs. By using the Arabidopsis thaliana chloroplast genome as the reference, we instruct `bbmap` to output only the unmapped reads (`outu=`), effectively filtering out the contaminating sequences.

```bash
for sample in 6
do
  bbmap.sh \
    in1=$sourcedir/$sample.cp.R1.trimmed.fq.gz \
    in2=$sourcedir/$sample.cp.R2.trimmed.fq.gz \
    ref=/mnt/s-ws/everyone/NC_000932.fa \
    outu1=$sample.cp.R1.trimmed.clean.fq.gz \
    outu2=$sample.cp.R2.trimmed.clean.fq.gz \
    minratio=0.9 maxindel=3 bwr=0.16 bw=12 fast minhits=2 \
    qtrim=r trimq=10 untrim idtag printunmappedcount kfilter=25 \
    maxsites=1 k=14
done
```
**Output Interpretation**
After running the `bbmap` command, you will see a summary of the results. Here is an example of what the output might look like:

```
java -ea -Xmx23776m -Xms23776m -cp /usr/share/java/bbmap.jar align2.BBMap build=1 overwrite=true fastareadlen=500 in1=/mnt/s-ws/everyone/SCIEM401/Module_6_Variants/6.cp.R1.trimmed.fq.gz in2=/mnt/s-ws/everyone/SCIEM401/Module_6_Variants/6.cp.R2.trimmed.fq.gz ref=/mnt/s-ws/everyone/NC_000932.fa outu1=6.cp.R1.trimmed.clean.fq.gz outu2=6.cp.R2.trimmed.clean.fq.gz minratio=0.9 maxindel=3 bwr=0.16 bw=12 fast minhits=2 qtrim=r trimq=10 untrim idtag printunmappedcount kfilter=25 maxsites=1 k=14
...
Mapping:                7.669 seconds.
Reads/sec:              51999.18
kBases/sec:             7650.01

   ------------------   Results   ------------------   

Genome:                 1
Key Length:             14
Max Indel:              3
Minimum Score Ratio:    0.9
Mapping Mode:           normal
Reads Used:             398798  (58670342 bases)

Pairing data:           pct pairs       num pairs       pct bases          num bases

mated pairs:              3.0868%            6155         3.1153%            1827756
bad pairs:                0.0015%               3         0.0015%                900
insert size avg:          618.40
unmapped:                88.6118%          353382        90.3187%           52990282
...
```

**Key Metrics:**

- **Reads Used**: This is the total number of reads that were processed.
- **Mapped Reads**: The percentage of reads that were successfully aligned to the reference genome. For example, `7.5251%` for Read 1 and `6.9514%` for Read 2.
- **Unmapped Reads**: The percentage of reads that did not align to the reference genome. This is what we are interested in, as these are the reads that were not contaminated with Arabidopsis DNA. For example, `88.6118%` of reads were unmapped.

**Questions to Consider:**

1. **What proportion of reads mapped to the Arabidopsis genome?**
   - This can be found in the "unmapped" percentage. Higher unmapped percentages indicate fewer reads were contaminated with Arabidopsis DNA.

2. **Did the percentage of unmapped reads meet your expectations?**
   - You should expect a high percentage of unmapped reads if the contamination was minimal.

3. **Why is it important to remove contaminant reads before further analysis?**
   - Contaminant reads can skew the results of downstream analyses, such as variant calling and gene expression studies.

By understanding these metrics and their implications, you can better interpret the quality and integrity of your sequencing data.


### Variant Calling with Cleaned Data

We have already performed the steps of mapping, sorting, indexing, and variant calling with the cleaned reads. The cleaned BAM files are available in `$sourcedir/clean_bams`.

**Steps and Output**

1. **Call Variants with the Clean BAM Files**

   The command used to call variants with the clean BAM files is as follows:

   ```bash
   bcftools mpileup -Ou -f Av.cp.final.fasta $sourcedir/clean_bams/*.sorted.bam | bcftools call --ploidy 1 -mv -Ov -o Av_clean.vcf
   ```

2. **Compare the Number of Variants Before and After Cleaning**

   To compare the number of variants before and after cleaning, the following commands were used:

   ```bash
   wc -l Av.vcf
   wc -l Av_clean.vcf
   ```

**Questions and Analysis**

1. **Has Cleaning the Reads Reduced the Number of Apparent Variants?**

   - **Yes, the number of apparent variants has been reduced.** The original VCF file (`Av.vcf`) had 463 lines, whereas the cleaned VCF file (`Av_clean.vcf`) has 301 lines. This reduction indicates that many of the variants detected in the original analysis were likely due to contamination.

2. **Why is it Important to Remove Contaminant Reads Before Variant Calling?**

   - **Accuracy**: Contaminant reads can introduce false-positive variants, leading to incorrect conclusions.
   - **Data Integrity**: Cleaned data ensures that the variants called are truly representative of the organism being studied, not contaminants.

3. **What Can We Infer from the Reduction in Variant Count?**

   - The significant reduction in the number of variants after cleaning suggests that the original data had a considerable amount of contamination. This emphasizes the importance of data preprocessing steps like read cleaning to achieve reliable results.

Examine the output to understand the impact of read cleaning on variant calling. By comparing the number of variants before and after cleaning, you can see the importance of this preprocessing step in achieving accurate results.

### Variant Querying and Annotation

In this section, we will focus on querying, filtering, and counting variants. We will filter variants to retain only SNPs or indels and only 'high-quality' variants. Follow the steps below to perform these tasks and analyze the results.

#### Querying and Filtering Variants

1. **Filter Variants**: Retain only high-quality SNPs (Single Nucleotide Polymorphisms) with a quality score greater than 30 and a depth greater than 20.

    ```bash
    bcftools filter -i 'TYPE="snp" && QUAL>30 && DP>20' Av_clean.vcf > Av.hiQ.vcf
    ```

2. **Count Variants**: Use `bcftools stats` to count and summarize the variants. The `grep "PSC"` command filters the output to show only per-sample counts.

    ```bash
    bcftools stats -s - Av.hiQ.vcf | grep "PSC"
    ```

3. **Analyze Outlier Sample**: Check if the cleaning process has resolved issues with any outlier samples by comparing the variant counts before and after filtering. [more to this here](extras/cleaning_snps.md)

4. **Query Specific Data**: Output a list of SNPs and how they differ from the reference sequence.

    ```bash
    bcftools query -i 'TYPE="snp"' -f '%CHROM %POS %REF %ALT\n' Av_clean.vcf
    ```

### Variant Querying and Annotation

In this section, we will focus on querying, filtering, and counting variants. We will filter variants to retain only SNPs or indels and only 'high-quality' variants. Follow the steps below to perform these tasks and analyze the results.

#### Querying and Filtering Variants

1. **Filter Variants**: Retain only high-quality SNPs (Single Nucleotide Polymorphisms) with a quality score greater than 30 and a depth greater than 20.

    ```bash
    bcftools filter -i 'TYPE="snp" && QUAL>30 && DP>20' Av_clean.vcf > Av.hiQ.vcf
    ```

2. **Count Variants**: Use `bcftools stats` to count and summarize the variants. The `grep "PSC"` command filters the output to show only per-sample counts.

    ```bash
    bcftools stats -s - Av.hiQ.vcf | grep "PSC"
    ```

3. **Analyze Outlier Sample**: Check if the cleaning process has resolved issues with any outlier samples by comparing the variant counts before and after filtering.

4. **Query Specific Data**: Output a list of SNPs and how they differ from the reference sequence.

    ```bash
    bcftools query -i 'TYPE="snp"' -f '%CHROM %POS %REF %ALT\n' Av_clean.vcf
    ```

#### Analysis and Interpretation

**Variant Distribution**:
- Compare the results before and after cleaning to determine if the cleaning process resolved issues with outlier samples.

**SNP Details**:
- List the specific SNPs and how they differ from the reference sequence.

**Gene Variation**:
- Identify which genes contain the most variants.
- Analyze the results to infer where most sequence variation occurs in the chloroplast genome.

#### Bedtools Intersect Analysis

Use the following command to count the number of variants in each gene:

```bash
bedtools intersect -c -a chloe/Av.cp.final.gff3 -b Av.hiQ.vcf | grep "CDS" | sort -nk 10
```

#### Sample Output Interpretation:

```plaintext
Av.cp.final     Chloe   CDS     103009  103176  ...   ID=rpl32.CDS.1;Parent=rpl32     0
Av.cp.final     Chloe   CDS     103669  104616  ...   ID=ccsA.CDS.1;Parent=ccsA       0
...
Av.cp.final     Chloe   CDS     108094  113754  ...   ID=ycf1.CDS.1;Parent=ycf1       5
Av.cp.final     Chloe   CDS     38139   40391   ...   ID=psaA.CDS.1;Parent=psaA       5
```

- **No Variants**: Most genes, such as `rpl32` and `ccsA`, contain no variants.
- **Few Variants**: Some genes, like `atpF`, have a small number of variants.
- **Most Variants**: Genes like `ycf1` and `psaA` contain more variants, with up to 5 variants each.

#### Questions for Students

1. **Variant Distribution**:
   - Compare the `nHapRef` and `nHapAlt` columns before and after cleaning. Has the number of variants decreased, indicating that the cleaning process was effective?

2. **SNP Details**:
   - Use the following command to list SNPs and their differences from the reference sequence:
     ```bash
     bcftools query -i 'TYPE="snp"' -f '%CHROM %POS %REF %ALT\n' Av_clean.vcf
     ```

3. **Gene Variation**:
   - Identify which genes, based on the Bedtools intersect output, contain the most variants.

4. **Variation Location**:
   - Most variants occur in non-coding regions, as coding sequences are generally under stronger evolutionary constraints.

🌟 **Woohoo, you've made it through my pracs! Congratulations!** 🎉