Here's an annotation for the command:

```bash
cut -d ',' -f 1 $sourcedir/cp_proteins.csv | tail -n +2 | uniq > gene_names.txt
```

### Explanation:

1. **`cut -d ',' -f 1 $sourcedir/cp_proteins.csv`**:
   - **`cut`**: Extracts sections from each line of files.
   - **`-d ','`**: Specifies the delimiter for fields in the input file, which is a comma (`,`). This indicates that the file is in CSV (Comma-Separated Values) format.
   - **`-f 1`**: Indicates that only the first field (column) should be extracted from each line of the file.
   - **`$sourcedir/cp_proteins.csv`**: Specifies the input file path using the `sourcedir` variable, which refers to the directory where `cp_proteins.csv` is located.

2. **`| tail -n +2`**:
   - **`|`**: Pipes the output of the `cut` command into the `tail` command.
   - **`tail -n +2`**: Outputs the file starting from the second line onward. The `+2` argument tells `tail` to skip the first line, which is typically the header line in CSV files.

3. **`| uniq`**:
   - **`|`**: Pipes the output of the `tail` command into the `uniq` command.
   - **`uniq`**: Removes duplicate lines from the input. It ensures that only unique values are retained in the output.

4. **`> gene_names.txt`**:
   - **`>`**: Redirects the output to a file.
   - **`gene_names.txt`**: The file where the output will be saved. It will contain the unique gene names extracted from the first column of the CSV file, excluding the header line.

### Summary:
This command sequence extracts the first column from a CSV file located in the specified directory, removes the header line, eliminates duplicate values, and saves the result (unique gene names) to a file named `gene_names.txt`.