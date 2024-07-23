Here’s a detailed annotation of the command:

```bash
awk '$3 == "gene"' Av.chloe.gff3 | cut -f 4,5,9 | sort -n | sed -E 's/[^\t]*;Name=([^;]*)/\1/' > Av.chloe.txt
```

### Breakdown

1. **`awk '$3 == "gene"' Av.chloe.gff3`**:
   - **`awk`**: A powerful text processing tool.
   - **`'$3 == "gene"'`**: An `awk` pattern that selects lines where the third field (`$3`) is equal to `"gene"`.
   - **`Av.chloe.gff3`**: The input file for `awk`.

   **Result**: This command filters lines in `Av.chloe.gff3` to include only those where the third field is `"gene"`.

2. **`cut -f 4,5,9`**:
   - **`cut`**: A command-line utility for cutting sections from each line of files.
   - **`-f 4,5,9`**: Specifies fields to be extracted (fields 4, 5, and 9) from the input.

   **Result**: After filtering by `awk`, this command extracts columns 4, 5, and 9 from the filtered lines.

3. **`sort -n`**:
   - **`sort`**: A command-line utility for sorting lines of text files.
   - **`-n`**: Sorts numerically.

   **Result**: Sorts the output numerically based on the values in the lines.

4. **`sed -E 's/[^\t]*;Name=([^;]*)/\1/'`**:
   - **`sed`**: Stream editor for parsing and transforming text.
   - **`-E`**: Enables extended regular expressions (ERE).
   - **`'s/[^\t]*;Name=([^;]*)/\1/'`**:
     - **`s/.../.../`**: Substitution command in `sed`.
     - **`[^\t]*`**: Matches any sequence of characters except for a tab character.
     - **`;Name=`**: Matches the literal string `;Name=`.
     - **`([^;]*)`**: 
       - **`(` and `)`**: Define a capture group.
       - **`[^;]*`**: Matches any sequence of characters except for a semicolon.
     - **`\1`**: Refers to the first capture group, which contains the text matched by `([^;]*)`.

   **Result**: This command transforms each line by replacing the portion before and including `;Name=` with the content captured after `;Name=`.

5. **`> Av.chloe.txt`**:
   - **`>`**: Redirects the output to a file.
   - **`Av.chloe.txt`**: The output file where the final result is saved.

   **Result**: The transformed data is written to `Av.chloe.txt`.

### Summary

The overall command performs the following steps:

1. Filters lines from `Av.chloe.gff3` where the third field is `"gene"`.
2. Extracts columns 4, 5, and 9 from these lines.
3. Sorts the extracted lines numerically.
4. Transforms each line to extract the part after `;Name=`, discarding the rest.
5. Saves the final output to `Av.chloe.txt`.

**Example**:

Given an input line:

```
something\t\tgene\t123\t456\t;Name=GeneName;otherdata
```

After processing, the output might be:

```
GeneName
```