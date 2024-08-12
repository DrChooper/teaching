
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
    - Use the `rm` command to remove the unsorted BAM files.The sorted and unsorted files contain the same content just differently packed and marked so you don;t need both files. You want to keep tha smaller file and remove the larger files. You can do this using regular expession patters:

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
[the output explained here](extras/bcftool1.md)

- To obtain a summary of a particular section, we use `bcftools stats`. By passing the output through `grep`, we can view the per-sample counts (PSC). 

```bash
bcftools stats -s - Av.vcf | grep "PSC"
```
However there are other sections ([explanation and output here](extras/bcftool1.md))


- The summary output will show columns such as `nHapRef` (same as reference) and `nHapAlt` (different from reference) since we specified the samples as haploid. If the samples were diploid, the variants would be counted in the `nRefHom` (homozygous reference), `nNonRefHom` (homozygous variant), and `nHets` (heterozygous) columns.

Here is an example of what the output might look like and how to interpret what you see: [output here](extras/bcftool1.md)

The conclusion is that there is too much 
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

The `wc -l` command shows the number of lines in each file:

- **Av.vcf** has 463 lines.
- **Av_clean.vcf** has 301 lines.

This indicates that after cleaning out the Arabidopsis sequences, 162 lines were removed from the original VCF file.

Now you can look at the cleaned file:([explanation and output here](extras/bcftool1.md))

```bash
bcftools stats -s - Av_clean.vcf | grep "PSC"
```

*There was some confusion about less variants and more variants after filtering. Let's clarify and reconcile the points that seem contradictory:*

#### Summary of Points:
1. **Reduction in Line Count (Apparent Variants):**
   - After cleaning the reads, the VCF file has fewer lines (301) compared to the original (463). This reduction suggests that many of the variants detected in the original data were due to contamination.

2. **nHapRef and nHapAlt Changes:**
   - **nHapRef**: The number of haploid reference alleles has decreased.
   - **nHapAlt**: The number of haploid alternate alleles has increased.

#### Reconciling the Points:
The reduction in the overall number of variants (as evidenced by the line count decrease) and the changes in nHapRef and nHapAlt can be understood as follows:

- **Contamination Removal and Variant Reduction**: 
   - The initial data contained contaminating sequences that were falsely identified as variants or reference alleles. Cleaning these contaminants removed many false-positive variants, leading to a decrease in the total number of variants (`nHapRef` and overall line count).

- **Improved Variant Detection Post-Cleaning**:
   - After cleaning, the data is more representative of the true chloroplast genome, which may have allowed for a more accurate identification of real variants that were previously masked by contaminants. As a result, **while the overall number of variants decreased**, the **proportion of true variants** (nHapAlt) relative to the total identified variants has increased.

#### Detailed Explanation:
1. **Initial State with Contamination**:
   - Contaminants may have aligned to regions that were either misidentified as reference alleles (increasing nHapRef) or introduced spurious variants (increasing the overall variant count).

2. **Post-Cleaning**:
   - Cleaning removed these contaminants, thus reducing the overall variant count. However, without the noise introduced by contaminants, the genuine variants (nHapAlt) became clearer and more pronounced, even if fewer in total number.

#### Conclusion:
- **Why are there fewer lines in the VCF file after cleaning?**
   - Because many of the apparent variants in the original file were false positives introduced by contamination.

- **Why is nHapAlt higher after cleaning?**
   - The cleaning process removed the noise, leading to a clearer identification of true variants. Though the total number of variants decreased, the proportion of true alternate alleles (nHapAlt) relative to the cleaned data has increased.

This process highlights the critical importance of removing contaminants to obtain accurate and meaningful variant calling results. By ensuring that the data reflects only the organism of interest, we get a more reliable representation of its genetic makeup.



### Variant Querying and Annotation

**We will now move on with the cleaned data output.** In this section, we will focus on querying, filtering, and counting variants. We will filter variants to retain only SNPs or indels and only 'high-quality' variants. Follow the steps below to perform these tasks and analyze the results.

#### Querying and Filtering Variants

1. **Filter Variants**: Retain only high-quality SNPs (Single Nucleotide Polymorphisms) with a quality score greater than 30 and a depth greater than 20.

    ```bash
    bcftools filter -i 'TYPE="snp" && QUAL>30 && DP>20' Av_clean.vcf > Av.hiQ.vcf
    wc -l Av.hiQ.vcf
    ```

   Comparison of variant counts:
   - Original VCF (Av.vcf): 463 variants
   - Cleaned VCF (Av_clean.vcf): 301 variants
   - High-quality VCF (Av.hiQ.vcf): 250 variants

   The cleaning process has reduced the number of variants, suggesting that some contaminant or low-quality variants have been removed.

2. **Count Variants**: Use `bcftools stats` to count and summarize the variants. The `grep "PSC"` command filters the output to show only per-sample counts.

    ```bash
    bcftools stats -s - Av.hiQ.vcf | grep "PSC"
    ```

  - Do we have more true variants now?
  - Which specimen is the most close or distant from our sample *Samples Sorted by Difference from Reference (which is our sample)*
   - **Most different:** `3.sorted.bam` with `145` alternate haploid alleles.
   - **Least different:** `5.sorted.bam` with `109` alternate haploid alleles.


3. **Analyze Outlier Sample**: Check if the cleaning process has resolved issues with any outlier samples by comparing the variant counts before and after filtering.
Here the overview of the 3 stages of the dataset that we looked at:
[overview of outputs](extras/overview_vcf.md)


4. **Query Specific Data**: Output a list of SNPs and how they differ from the reference sequence.
How to: Use the `bcftools` to extract SNPs from cleaned or hQC file and then compare

    ```bash
    bcftools query -i 'TYPE="snp"' -f '%CHROM %POS %REF %ALT\n' Av_clean.vcf
    ```

[check here](extras/snp_compare.md)

5. **Gene Variation**:
- Identify which genes contain the most variants.
- Analyze the results to infer where most sequence variation occurs in the chloroplast genome.

#### Interpreting the `bedtools intersect` Output

```bash
bedtools intersect -c -a chloe/Av.cp.final.gff3 -b Av.hiQ.vcf | grep "CDS" | sort -nk 10
```

This command is a powerful way to analyze how variants intersect with coding sequences (CDS) in a genome. Here's a breakdown of what the command does and how to interpret the output:

#### Command Breakdown:
1. **`bedtools intersect -c`**: 
   - `bedtools intersect` is a tool that finds overlaps between two sets of genomic features.
   - The `-c` option counts the number of overlaps between features in the first file (`-a`) and those in the second file (`-b`).

2. **`-a chloe/Av.cp.final.gff3`**:
   - This is the genome annotaiton file we have created Monday (`Av.cp.final.gff3`). This contains annotations of genomic features, including coding sequences (CDS), exons, introns, etc., for the chloroplast genome.

3. **`-b Av.hiQ.vcf`**:
   - This specifies the second file, which is in VCF format (`Av.hiQ.vcf`). This file contains high-quality variants (likely SNPs) identified in the chloroplast genome.

4. **`grep "CDS"`**:
   - This filters the output to include only lines that correspond to CDS features, which represent protein-coding regions in the genome.

5. **`sort -nk 10`**:
   - This sorts the filtered output numerically (`-n`) by the 10th column (`-k 10`), which contains the count of intersecting variants.

#### Output Interpretation:

The output will show each CDS feature from the GFF3 file, along with the count of high-quality variants that intersect with each CDS. The sorted output will list these features from the lowest to the highest number of intersecting variants.

#### Example of Output:
```plaintext
...
Av.cp.final     Chloe   CDS     82316   83140   9.001e-03       -       0       ID=rpl2-2.CDS.1;Parent=rpl2-2   0
Av.cp.final     Chloe   CDS     83159   83440   1.554e-03       -       0       ID=rpl23-2.CDS.1;Parent=rpl23-2 0
Av.cp.final     Chloe   CDS     83819   89932   6.275e-03       +       0       ID=ycf2.CDS.1;Parent=ycf2       0
Av.cp.final     Chloe   CDS     91801   92058   2.690e-03       +       0       ID=rps12B.CDS.2;Parent=rps12B   0
Av.cp.final     Chloe   CDS     92117   92584   3.971e-04       +       0       ID=rps7.CDS.1;Parent=rps7       0
Av.cp.final     Chloe   CDS     13145   13888   6.512e-04       -       0       ID=atpI.CDS.1;Parent=atpI       1
Av.cp.final     Chloe   CDS     19377   20993   1.891e-03       -       0       ID=rpoC1.CDS.3;Parent=rpoC1     1
Av.cp.final     Chloe   CDS     22171   25383   2.825e-04       -       0       ID=rpoB.CDS.1;Parent=rpoB       1
Av.cp.final     Chloe   CDS     66322   66702   1.124e-03       +       0       ID=rpl20.CDS.1;Parent=rpl20     1
Av.cp.final     Chloe   CDS     8921    10444   3.111e-04       -       0       ID=atpA.CDS.1;Parent=atpA       1
Av.cp.final     Chloe   CDS     402     1925    3.854e-03       -       0       ID=matK.CDS.1;Parent=matK       2
Av.cp.final     Chloe   CDS     108094  113754  5.212e-02       -       0       ID=ycf1.CDS.1;Parent=ycf1       5
Av.cp.final     Chloe   CDS     38139   40391   7.808e-05       -       0       ID=psaA.CDS.1;Parent=psaA       5
...
```

#### Explanation:
- **Columns 1-9**: These columns represent the standard fields in a GFF3 file, including chromosome, source, feature type (CDS in this case), start and end positions, score, strand, phase, and attribute fields.
  
- **Column 10**: This is the count of high-quality variants from the VCF file that intersect with the corresponding CDS.

- **Sorting**: The sorting by the 10th column allows you to quickly identify which CDS features have the most (or fewest) intersecting variants. 

This output helps to identify regions in the chloroplast genome that are most affected by variants, particularly in protein-coding regions. Regions with higher variant counts might be of particular interest for further functional analysis, as they could indicate areas of high mutation or regions subject to selective pressures.

#### Conclusion from out dataset

1. **No Overlap with High-Quality Variants**:
   - The majority of the CDS regions listed in the output have a count of `0` in the last column. This indicates that these coding sequences do not overlap with any high-quality variants in the `Av.hiQ.vcf` file.
   - This suggests that these regions of the chloroplast genome are either highly conserved or that the variants in the dataset are located outside of these coding regions.

2. **Few CDS Regions with Variants**:
   - A small number of CDS regions show a count of `1` or more, indicating the presence of one or more variants overlapping these coding sequences. These regions might be of particular interest as they could indicate regions with potential functional significance or evolutionary importance.

3. **Specific Genes with Variants**:
   - For example, the gene **`matK`** has 2 variants within its CDS region, and **`ycf1`** has 5 variants. The presence of these variants could suggest potential polymorphisms or evolutionary changes within these genes.

4. **Sorting by Variant Count**:
   - The output is sorted by the number of intersecting variants, with genes that have more variants listed at the bottom. This sorting helps to quickly identify which coding sequences have the most variants, which could be useful for prioritizing further analysis.

### Implications:
- **Conservation of Chloroplast Genes**: The fact that most coding sequences do not overlap with high-quality variants suggests that many chloroplast genes are highly conserved, reflecting their essential roles in the plant's cellular machinery.
- **Potential Functional Impact**: The few coding sequences that do have variants may be of interest for further study to understand the potential impact of these variants on gene function and plant physiology.


🌟 **Woohoo, you've made it through my pracs! Congratulations!** 🎉