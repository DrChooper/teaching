## To compare the SNPs between the cleaned and filtered VCF files, you can follow these steps:

### **1. Extract SNP Data from Both VCF Files**

First, extract the SNP data from both the cleaned and filtered VCF files:

```bash
bcftools query -i 'TYPE="snp"' -f '%CHROM %POS %REF %ALT\n' Av_clean.vcf > snps_cleaned.txt
bcftools query -i 'TYPE="snp"' -f '%CHROM %POS %REF %ALT\n' Av.hiQ.vcf > snps_filtered.txt
```

This will generate two text files (`snps_cleaned.txt` and `snps_filtered.txt`) that list the SNPs in each dataset.

### **2. Compare the SNP Data**

You can now compare these two lists using `diff`, `comm`, or any other comparison tool:


To see the SNPs that are common in both `snps_cleaned.txt` and `snps_filtered.txt`, as well as those that are unique to the cleaned file, you can use the `comm` command with the following options:

### Step 1: Find Common SNPs in Both Files
To see the SNPs that are common in both the cleaned and filtered files:

```bash
comm -12 <(sort snps_cleaned.txt) <(sort snps_filtered.txt)
```

### Step 2: Find SNPs Only in the Cleaned File
To see the SNPs that are present only in the cleaned file but not in the filtered file:

```bash
comm -23 <(sort snps_cleaned.txt) <(sort snps_filtered.txt)
```

### Explanation:
- **`comm -12`**: This option shows only the lines that are common between the two sorted files.
- **`comm -23`**: This option shows the lines that are unique to the cleaned file (i.e., present in `snps_cleaned.txt` but not in `snps_filtered.txt`).

### Example Usage:
1. **Common SNPs**:
   ```bash
   comm -12 <(sort snps_cleaned.txt) <(sort snps_filtered.txt)
   ```

   This will output the SNPs that are found in both the cleaned and filtered datasets.

2. **Unique SNPs in Cleaned**:
   ```bash
   comm -23 <(sort snps_cleaned.txt) <(sort snps_filtered.txt)
   ```

   This will output the SNPs that are only in the cleaned dataset and were filtered out in the subsequent processing.

### Result Interpretation:
- **Common SNPs**: These are the variants that passed through both cleaning and filtering, suggesting they are higher confidence calls.
- **Unique SNPs in Cleaned**: These variants did not meet the filtering criteria, which could be due to lower quality, insufficient read depth, or other factors that were used as thresholds in the filtering process. 

This approach will give you a clear view of which SNPs were retained through filtering and which were excluded.

### **3. Analyze the Differences**

- **Identify Missing/Additional SNPs**: Review the output to see if certain SNPs were removed or added during filtering.
- **Quality of SNPs**: Filtering usually retains high-quality SNPs, so any discrepancies might indicate lower-quality SNPs that were removed or new high-quality SNPs that were identified.

### **4. Summary**
- **Increased Consistency**: If the filtered SNPs are more consistent across samples, it suggests that filtering has improved the data quality.
- **Removed Ambiguities**: Differences might also show which SNPs were ambiguous or low-quality, leading to their removal in the filtered data.