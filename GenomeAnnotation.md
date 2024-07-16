# Genome Annotation Workshop

## Content

1. [GeSeq Annotation](#geseq-annotation)
2. [Chloe Annotation](#chloe-annotation)
3. [Data Extraction and Preprocessing](#data-extraction-and-preprocessing)
4. [Gene Counting and Comparison](#gene-counting-and-comparison)
5. [Detailed Annotation Comparison](#detailed-annotation-comparison)
6. [Extracting Sequence Features](#extracting-sequence-features)
7. [Aligning and Analyzing Sequences](#aligning-and-analyzing-sequences)


## Introduction
The Genome Annotation Workshop focuses on annotating the Aldrovanda vesiculosa chloroplast genome using GeSeq and Chloë, highlighting their different approaches and output comparisons. Participants will manipulate GFF files and use text-wrangling techniques to prepare data for analysis and phylogenetic studies.

We'll start by annotating the Aldrovanda vesiculosa chloroplast genome with two different software tools in the browsers.
- GeSeq
- Chloe

### GeSeq Annotation
- Geseq uses a database of features to identify equivalent features in the target genome
- It finds protein-coding sequences using tblastx (using a protein sequence to align to translated nucleotide sequence)
- This makes it very sensitive at detecting such features, but it may struggle to accurately place the edges of the features in some cases

1. Upload the the genome sequence(`Av.cp.final.fasta`) to be found in `/mnt/s-ws/everyone/Aldrovanda_vesiculosa/Av.cp.final.fasta` to the [Geseq](https://chlorobox.mpimp-golm.mpg.de/geseq.html) website https://chlorobox.mpimp-golm.mpg.de/geseq.html 
    - You need the fasta file on your computer. You can either download it from LMS or you copy it down from the server.From your computer in a chosen folder in the terminal window you can use the `scp` command and when prompted give you password.

    ```bash
    #please use your account and your server
    scp studentaccount@yourserverip:/mnt/s-ws/everyone/Aldrovanda_vesiculosa/Av.cp.final.fasta .
    ```

2. Use the settings indicated in the Geseq slide in the introductory presentation
3. When Geseq has finished the annotations, download the GFF file to your PC and then upload it to the Nimbus server. 
    - In the process, rename it to something simpler, e.g. Av.geseq.gff3
    - Uploaded to your home folder on the server:
    ```bash
    scp Av.geseq.gff3 studentaccount@yourserverip:~/
    ```

### Chloe Annotation
- Chloë (developed at UWA) uses a different approach, based on whole genome alignments, not feature alignments.
- This allows more accurate placement of features, but at some cost in sensitivity

-  Upload the same genome sequence (/mnt/s-ws/everyone/Aldrovanda_vesiculosa/Av.cp.final.fasta) to the [Chloë](https://chloe.plastid.org/annotate.html) website https://chloe.plastid.org/annotate.html
- Use the settings indicated in the Chloë slide in the introductory presentation
- When Chloë has finished the annotations, download the GFF file to your PC and then upload it to the Nimbus server
- In the process, rename it to something simpler, e.g. Av.chloe.gff3

### Data Extraction and Preprocessing
Take a look at the two GFF files. Although they both follow the GFF3 specification, they are not easy to compare, because the feature names are presented differently.

```bash
less Av.geseq.gff3
```
```bash
less Av.chloe.gff3
```
- Thus to facilitate comparisons, we are going to extract the information we need from each of these two files.
- This sort of 'text-wrangling' is very commonly required in bioinformatics (or in data analysis generally) and thus is a key transferable skill to learn.

- `awk` is a powerful text processor that can be used to filter files
- The following command uses awk to extract lines in a file where the 3rd field is 'gene'
- Then we use cut to extract just columns 4, 5 and 9
```bash
awk '$3 == "gene"' Av.geseq.gff3 | cut -f 4,5,9 | sort -n > Av.geseq.txt
head Av.geseq.txt
```
- We still need to process the problematic column 9 which includes the feature names.
- For this we are going to use a search pattern known as a 'regular expression' or 'regex' (see https://en.wikipedia.org/wiki/Regular_expression)
- The pattern we are going to use is formed of three parts
    1. `[^\t]*;` indicates a run of characters that are not tab characters, followed by a semi-colon
    2. `gene=([^;]*);` indicates the characters gene= followed by a run of characters that are not semi-colons, followed by a semi-colon.
        -The parentheses are not part of the search pattern, but indicate that we want to capture the text that matches this part of the pattern (because this text is the name of the gene)
    3. `.*` indicates a run of any characters

- We are going to use these patterns to direct the program sed
- `sed` is a stream editor used to alter text (e.g. sed 's/find/replace/' input.txt > output.txt will replace all occurrences of 'find' with 'replace' in input.txt
- The command below replaces the text that matches our regex pattern with the captured text (the gene name)

```bash
sed -E -i 's/[^\t]*;gene=([^;]*);.*/\1/' Av.geseq.txt
``` 

- Check that this has worked as expected, i.e. replaced the unintelligible final column with just the gene name
```bash
head Av.geseq.txt
```

🚀 If you have any intentions of a career in bioinformatics or data analysis, practising the use of regular expressions to filter and reformat text files will prove invaluable!🔍✨

### Gene Counting and Comparison
- Now we need to do the same sort of text-wrangling with the Chloe output. Note that the patterns are slightly different because the gene names are presented differently in the Chloe output

```bash
awk '$3 == "gene"' Av.chloe.gff3 | cut -f 4,5,9 | sort -n | sed -E 's/[^\t]*;Name=([^;]*)/\1/' > Av.chloe.txt
head Av.chloe.txt
```
- Now we can count how many distinct genes Chloe found.

```bash
cut -f 3 Av.chloe.txt | uniq | wc -l
```

- We can see which genes are different using `comm`
-This command is using a bash trick called process substitution (https://www.gnu.org/software/bash/manual/html_node/Process-Substitution.html)
- In the output, the first column contains genes found only by Geseq, column two contains genes found only by Chloe
```bash
comm -3 <(cut -f 3 Av.geseq.txt) <(cut -f 3 Av.chloe.txt)
```
- Did Chloe find fewer genes? It generally does, for two reasons; the tblastx search used by GeSeq is more sensitive, so can find genes that Chloe misses (false negative by Chloe).
- However, it is also prone to detecting remnants (pseudogenes) that are no longer functional (false positive by Geseq).
- Distinguishing between 'real' genes (i.e. that are functionally significant) and features that can look like genes but that have little or no functional significance (pseudogenes, transposable elements)
- is a major ongoing issue with genome annotation.

### Detailed Annotation Comparison
- We can look at differences between the annotations in more detail using `diff`
- `diff` catalogues differences between files
- `diff` output is explained here https://www.gnu.org/software/diffutils/manual/diffutils.html#Detailed-Unified
```bash
diff -u Av.geseq.txt Av.chloe.txt
```
- You can compare these annotations visually using IGV
- Download IGV from https://software.broadinstitute.org/software/igv/download
- Load 'Av.cp.final.fasta' into IGV (Genomes->Load Genome from File...)
- Load 'Av.geseq.gff3' into IGV (File->Load from File...)
- Load 'Av.chloe.gff3' into IGV (File->Load from File...)

- You can see that there are quite a few differences in the annotations produced by these two software pipelines, even though both are supposed to be optimised for annotating chloroplast genomes
- Hence the importance, when comparing annotations, of ensuring that all genomes have been annotated consistently with the same tools.
- Otherwise you may find more artefacts due to differences in the annotation pipelines than real biological differences.

**We are now going to compare annotations across several genomes, some closely related to Aldrovanda, others only distantly related. Many of them are also carnivorous plants**
- You can view some basic information on these plants
- column is a useful unix command for tabulating data

- here we are setting a variable 'sourcedir' to the directory that contains all the data for today's lab
```bash
sourcedir=/mnt/s-ws/everyone/annotation

cat $sourcedir/genome_metadata.tsv | column -t -s $'\t' | less -S
```
- You can also download this file to your PC and open it in Excel if you prefer
- Although these are published genomes, downloaded from NCBI GenBank, we can't rely on the annotations provided with these genomes as being correct,
- or even consistent (i.e. annotated in the same way)
- So we're going to re-annotate them with Chloe, using the command line version on the Nimbus server
```bash
mkdir chloe

julia /mnt/s-ws/everyone/chloe_biojulia/chloe.jl annotate -g -o chloe/ $sourcedir/*.fasta
```
- This should have generated 11 GFF files in the chloe directory

- I've provided a file listing the protein-encoding genes in flowering plant chloroplast genomes and some basic information about the proteins they encode

```bash
cat $sourcedir/cp_proteins.csv | column -t -s, | sort | less -S
```

- From this we can extract a list of gene names
```bash
cut -d ',' -f 1 $sourcedir/cp_proteins.csv | tail -n +2 | uniq > gene_names.txt

less gene_names.txt
```
- We can use this list to query all the gff files and count how many of the 11 genomes contain each gene
- This short script reads the gene_names.txt file and for each gene in the file searches all 11 GFF files using grep
- The gene count is saved in gene_numbers.txt

```bash
while read -r gene
do
    numgenes=$(grep "ID=$gene;Name=$gene" chloe/*.gff3 | wc -l)
    echo "$gene $numgenes" >> gene_numbers.txt
done < gene_names.txt
```

- Look at the results
``bash
less gene_numbers.txt
```

- You'll see that nearly all genes are found in all 11 genomes
- Some aren't, and we can sort the list to see which genes are found least often
```bash
sort -nk 2 gene_numbers.txt | head
```
- The ndh genes are each only present in 7 of the 11 genomes
- These genes all encode subunits of the plastoquinone oxido-reductase complex
```bash
grep "ndh" $sourcedir/cp_proteins.csv
```
- Are the losses random, or have the same 4 plants lost all of these genes?
- Check for a few different ndh genes using grep

```bash
grep "ID=ndhA;Name=ndhA" chloe/*.gff3
```


### Extracting Sequence Features
- Which 4 plants have lost all their ndh genes? What do they have in common?

- Very often, we would like to use annotations to analyse or extract different features
- For example, in a later module we will be constructing phylogenetic trees of these genomes
- and for that we need to extract a consistent set of sequences to analyse

- As an example, here's a short script that will calculate the average length of each gene across the 11 genomes

```bash
while read -r gene
do
    sumlength=0
    grep "ID=$gene;Name=$gene" chloe/*.gff3 > $gene.gff
    count=$(wc -l < $gene.gff)
    while read -r a b c begin end f g h i
    do
        sumlength=$((sumlength + end - begin + 1))
    done < $gene.gff
    rm $gene.gff
    echo "$gene $((sumlength/count))"  >> gene_lengths.txt
done < gene_names.txt

sort -nk 2 gene_lengths.txt | less
```

- scroll to the end of the file; you'll see that the genes range in average size from under 100 bp to over 6 kb
- For phylogenetics, the more sequence we can use, the better, provided the sequences are homologous and align well.
- Pick a long (> 1000 bp) gene that is present in all 11 genomes (check gene_numbers.txt)
- Extract the gene coordinates from the GFF files, e.g. here using the example of rpoC2

```bash
grep -h "ID=rpoC2;Name=rpoC2" chloe/*.gff3 | cut -f 1,4,5,7 > rpoC2.tsv
less rpoC2.tsv
```
- Check the lengths of the gene across genomes are consistent (i.e. the difference between the two columns is roughly constant)
- and that the length is divisible by 3 (which it should be if it's full-length)

```bash
while read -r id start end strand
do
length=$((end - start + 1))
modulus=$((length % 3))
echo "$id $length $modulus"
done < rpoC2.tsv
```

- If the lengths are not roughly consistent (Â± 10%), or if the modulus is not zero for any of the genes, pick a different gene and try again

- Now we can extract the gene sequences
- First copy the sequences to your workspace, e.g. in the chloe directory where you have the annotations

```bash
cp $sourcedir/*.fasta chloe/
```
- Add your Av Assembly

```bash
cp Av.cp.final.fasta chloe/
```
- We'll use a homemade Julia script to extract the gene sequences

```bash
julia $sourcedir/extract_from_fasta.jl chloe rpoC2.tsv
```
- This should generate a nucleotide fasta file containing the rpoC2 gene sequence for each genome.
- Ideally we would use the nearly universally available samtools faidx for feature extraction instead of a homemade script.
- But samtools faidx cannot correctly extract features that cross the ends of the genome, which genes in circular genomes often do.

- Check that the extracted sequences start with a start codon and end with a stop codon
```bash
head Av.cp.final.rpoC2.nt.fa

tail Av.cp.final.rpoC2.nt.fa
```

### Aligning and Analyzing Sequences
- We can translate these sequences into protein sequences
```bash
julia /mnt/s-ws/everyone/tools/translatefasta.jl *.rpoC2.nt.fa

less Av.cp.final.rpoC2.nt.protein.fa
```

- We can compare these sequences by doing a multiple sequence alignment
- First concatenate the protein sequence files then align with MAFFT https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3603318/
```bash
cat *rpoC2.nt.protein.fa > allrpoC2.protein.fa

mafft --maxiterate 1000 --globalpair allrpoC2.protein.fa > rpoC2.protein.msa
```
- You can upload the .msa file to online viewers, e.g.
- https://www.ncbi.nlm.nih.gov/projects/msaviewer/
- https://www.ebi.ac.uk/Tools/msa/mview/

**This ends the annotation lab; you will need several of the files you created for the phylogenetics lab, so don't delete them yet.**


