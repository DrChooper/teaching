
# Variant Calling

## Content
1. [Mapping Reads to Reference using bbmap](#mapping-reads-to-reference-using-bbmap)
2. [Quality Control and Clean up](#quality-control-and-clean-up)
3. [Analysis and Interpretation](#analysis-and-interpretation)

## Introduction
This workflow details the analysis of chloroplast DNA sequencing data, encompassing read mapping, variant calling, quality control, contamination removal, and annotation to identify genetic variations across the genome.

### Mapping Reads to reference using bbmap

- Map reads to reference, creating a SAM file. SAM stands for *Sequence Alignment/Map*.
- The reads have already been trimmed to remove adapter sequences
- We are using bbmap as the alignment software https://sourceforge.net/projects/bbmap/

```bash
sourcedir=/mnt/s-ws/everyone/SCIEM401/Module_6_Variants
```
- map reads from variants, creating BAM files

```bash
for sample in {1..2}
do
bbmap.sh in1=$sourcedir/$sample.cp.R1.trimmed.fq.gz in2=$sourcedir/$sample.cp.R2.trimmed.fq.gz out=$sample.bam ref=Av.cp.final.fasta mappedonly=t
done
```
- sort and index the BAM files using samtools (https://github.com/samtools/samtools)

```bash
for sample in {1..2}
do
samtools sort -T tmp -o $sample.sorted.bam $sample.bam
samtools index $sample.sorted.bam
done
```
- compare the sizes of the sorted and unsorted BAM files

```bash
ls -lh
```
- are the sorted files generally smaller or larger?

- remove the redundant unsorted BAM files
- note the use of `[0-9]` to indicate any number between 0 and 9
- why can't we simply use `rm *.bam` ?
- (Alternatively, you could use `rm ?.bam` and `rm ??.bam`)

```bash
rm [0-9].bam
```

**Variant Calling with Original Data**
- Use `bcftools` to call variants. This will take several minutes.
- The output file is in VCF format, described here: http://samtools.github.io/hts-specs/
*As for SAM/BAM, there is a binary (compressed) format, BCF if you are generating very large files.*

- In this module, as the VCF files are small, we'll stick to VCF format
- By default, bcftools assumes the samples are from a diploid organism and thus may be heterozygous for variants
- But in this case, the samples are of chloroplast DNA, which like bacterial DNA, can be safely considered to be uniform
- So we can treat these samples as effectively haploid (--ploidy 1)

```bash
bcftools mpileup -Ou -f Av.cp.final.fasta $sourcedir/bams/*.sorted.bam | bcftools call --ploidy 1 -mv -Ov -o Av.vcf
```

- We can see a summary of the results by using bcftools stats
- And just view the per sample counts by passing the output through grep

```bash
bcftools stats -s - Av.vcf | grep "PSC"
```

- As we told bcftools the samples were effectively haploid, all the variants are counted in the nHapRef (same as reference) or nHapAlt (different from reference) columns
- If the samples were diploid, variants would be counted in the nRefHom (homozygous reference), nNonRefHom (homozygous variant) and nHets (heterozygous) columns.

- Which sample shows the least variation to the reference?
- Which sample shows the most?
- One sample shows a surprisingly large number of variations. This might indicate contamination of the sample with sequencing reads from another plant.


### Quality Control and Clean up

**Removal of Arabidopsis Contamination**
- We're going to use `bbmap` to remove reads that map to the Arabidopsis thaliana chloroplast genome.
- Arabidopsis is a likely contaminant for sequencing experiments as it is the most widely used plant in research labs.
- This time we use the Arabidopsis thaliana chloroplast genome as the reference, and we ask bbmap to output unmapped reads (outu=), not mapped ones

```bash
for sample in 6
do
bbmap.sh in1=$sourcedir/$sample.cp.R1.trimmed.fq.gz in2=$sourcedir/$sample.cp.R2.trimmed.fq.gz ref=/mnt/s-ws/everyone/NC_000932.fa outu1=$sample.cp.R1.trimmed.clean.fq.gz outu2=$sample.cp.R2.trimmed.clean.fq.gz minratio=0.9 maxindel=3 bwr=0.16 bw=12 fast minhits=2 qtrim=r trimq=10 untrim idtag printunmappedcount kfilter=25 maxsites=1 k=14
done
```

- What proportion of reads mapped to the Arabidopsis genome?

- I've redone the mapping, sorting, indexing and variant calling steps with the cleaned reads
- map reads from all variants, creating new BAM files in $sourcedir/clean_bams

**Variant Calling with Cleaned Data**
- call variants with the clean bam files

```bash
bcftools mpileup -Ou -f Av.cp.final.fasta $sourcedir/clean_bams/*.sorted.bam | bcftools call --ploidy 1 -mv -Ov -o Av_clean.vcf
```

- has cleaning the reads reduced the number of apparent variants?

```bash
wc -l Av.vcf

wc -l Av_clean.vcf
```

### Variant Querying and Annotation

- querying, filtering and counting variants
- we can filter variants to only retain SNPs or indels, and only 'high-quality' variants

```bash
bcftools filter -i 'TYPE="snp" && QUAL>30 && DP>20' Av_clean.vcf > Av.hiQ.vcf

bcftools stats -s - Av.hiQ.vcf | grep "PSC"
```

- Has the cleaning process fixed the problems with the outlier sample?

- We can query the VCF file to output just the data we want, e.g. a list of SNPs and how they differ from the reference sequence

```bash
bcftools query -i 'TYPE="snp"' -f '%CHROM %POS %REF %ALT\n' Av_clean.vcf
```

**Analysis and Interpretation**
- We can see where the variation is occurring by comparing the VCF file to genome annotations
- We can use the GFF file you generated in the Annotation Module
- Bedtools is a software package that can do many different operations on genomics data, including VCF, BAM and GFF files
- https://bedtools.readthedocs.io/en/latest/
- Here we are using bedtools intersect to count the number of variants in each gene

```bash
bedtools intersect -c -a Av.cp.final.gff3 -b Av.hiQ.vcf | grep "CDS" | sort -nk 10
```

- which genes contain the most variants?
- Most genes contain no variants, and the sum of the variants within genes is much lower than the total number of variants found by bcftools
- As the gene-dense chloroplast genome is about 50% coding sequence, what does this imply about where most sequence variation occurs?

🌟 **Woohoo, you've made it through my pracs! Congratulations!** 🎉
