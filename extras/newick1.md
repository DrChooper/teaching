## Script Annotation

### Overview
This script processes a Newick tree string by replacing accession numbers with corresponding tax IDs from a TSV file.

### Steps

1. **Set Directory Path**
   ```bash
   sourcedir="/mnt/s-ws/everyone/annotation/"
   ```
   - Specifies the directory where the `genome_metadata.tsv` file is located.

2. **Initialize Newick Tree String**
   ```bash
   tree="((((((NC_041258:0.020608038493454474,NC_051455:0.01917829831140714)..."
   ```
   - Defines the Newick formatted tree with accession numbers to be replaced.

3. **Read TSV and Update Tree**
   ```bash
   while IFS=$'\t' read -r latin name taxid accession carnivorous australian
   do
       tree=$(sed "s/$accession/$taxid/g" <<< "$tree")
   done < <(tail -n +2 "$sourcedir/genome_metadata.tsv")
   ```
   - Reads each line from the TSV file, excluding the header.
   - Replaces each accession number in the tree string with its corresponding tax ID using `sed`.

4. **Output the Modified Tree**
   ```bash
   echo "$tree"
   ```
   - Prints the updated Newick tree with tax IDs.

### Notes
- The TSV file should be tab-separated and contain the columns: `latin`, `name`, `taxid`, `accession`, `carnivorous`, `australian`.
- Ensure the header row is excluded by using `tail -n +2`.

This script efficiently updates a Newick tree by substituting accession numbers with tax IDs from a metadata file.

---

#### Mashtree
HERE the actual example that we ran
```bash
# Directory containing the input file
sourcedir="/mnt/s-ws/everyone/annotation/"

# Newick tree string placeholder
tree="((((((NC_041258:0.020608038493454474,NC_051455:0.01917829831140714)Node_5:0.003407144767460306,(((NC_035415:0.0029979458688291823,triffid0104:0.031308218707373964)Node_1:0.00989078007801644,NC_035417:0.03143793544625087)Node_2:0.0029300500394906642,Av.cp.final:0.02453883519364107)Node_3:0.028925162654889328)Node_6:0.0022608164161328145,NC_010776:0.04037827266534525)Node_7:0.007501864123288876,NC_008326:0.035398386435476964)Node_8:0.001897550306508941,NC_022402:0.029023260143435453)Node_9:0.0009020351129601575,(NC_042597:0.026622618724706764,NC_051971:0.030041951163609974)Node_4:0.00686151403585894,NC_026134:0.039746905986761394)Node_10:1.0;"

# Process the TSV file to replace accession numbers with tax IDs
while IFS=$'\t' read -r latin name taxid accession carnivorous australian
do
    # Replace accession with tax ID in the tree string
    tree=$(sed "s/$accession/$taxid/g" <<< "$tree")
done < <(tail -n +2 "$sourcedir/genome_metadata.tsv")

# Print the modified Newick tree
echo "$tree"
```

Or you can use the species names and that will look eye friendly
```bash
# Directory containing the input file
sourcedir="/mnt/s-ws/everyone/annotation/"

# Newick tree string placeholder
tree="((((((NC_041258:0.020608038493454474,NC_051455:0.01917829831140714)Node_5:0.003407144767460306,(((NC_035415:0.0029979458688291823,triffid0104:0.031308218707373964)Node_1:0.00989078007801644,NC_035417:0.03143793544625087)Node_2:0.0029300500394906642,Av.cp.final:0.02453883519364107)Node_3:0.028925162654889328)Node_6:0.0022608164161328145,NC_010776:0.04037827266534525)Node_7:0.007501864123288876,NC_008326:0.035398386435476964)Node_8:0.001897550306508941,NC_022402:0.029023260143435453)Node_9:0.0009020351129601575,(NC_042597:0.026622618724706764,NC_051971:0.030041951163609974)Node_4:0.00686151403585894,NC_026134:0.039746905986761394)Node_10:1.0;"

# Replace accession with species name in the tree string
while IFS=$'\t' read -r species _ _ accession _ _
do
    tree=$(sed "s/$accession/$species/g" <<< "$tree")
done < <(tail -n +2 "$sourcedir/genome_metadata.tsv")

# Print the modified Newick tree
echo "$tree"
```

#### IQ-TREE

Here the actual tree:

```bash
# Directory containing the input file
sourcedir="/mnt/s-ws/everyone/annotation/"

# Newick tree string placeholder
tree="((((((Av.cp.final.rpoC2:0.1052793967,NC_035417.rpoC2:0.1587527461)100/100:0.0161788686,(NC_035415.rpoC2:0.0007903547,triffid0104.rpoC2:0.0838105936)100/100:0.0466704676)100/100:0.1431645840,NC_041258.rpoC2:0.0766205550)100/100:0.0073415864,NC_051455.rpoC2:0.1008749587)100/100:0.0456123691,NC_010776.rpoC2:0.2761098574)100/100:0.0947961713,NC_008326.rpoC2:0.1905579329,(NC_022402.rpoC2:0.1873226950,(NC_026134.rpoC2:0.2368698654,(NC_042597.rpoC2:0.1370663007,NC_051971.rpoC2:0.1932912155)100/100:0.0535984419)100/100:0.0256457427)100/100:0.0310506528);"

# Replace accession with species name in the tree string
while IFS=$'\t' read -r species _ _ accession _ _
do
    tree=$(sed "s/$accession/$species/g" <<< "$tree")
done < <(tail -n +2 "$sourcedir/genome_metadata.tsv")

# Print the modified Newick tree
echo "$tree"
```

#### MrBayes

Here the actual tree:

```bash
# Nexus tree not yet annotated
tree="[&U] (1[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:3.545467e-02[&length_mean=3.50257832e-02,length_median=3.54546700e-02,length_95%HPD={2.90613900e-02,3.97492700e-02}],7[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:5.315716e-02[&length_mean=5.28309676e-02,length_median=5.31571600e-02,length_95%HPD={4.72365000e-02,5.85820900e-02}],(((((2[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:6.604020e-02[&length_mean=6.51734753e-02,length_median=6.60402000e-02,length_95%HPD={5.26125000e-02,7.31287800e-02}],(4[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:6.220791e-02[&length_mean=6.22263394e-02,length_median=6.22079100e-02,length_95%HPD={5.58731400e-02,6.75169500e-02}],(5[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:7.821340e-02[&length_mean=7.98002797e-02,length_median=7.82134000e-02,length_95%HPD={7.23440700e-02,8.94922500e-02}],(9[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:4.610593e-02[&length_mean=4.63465888e-02,length_median=4.61059300e-02,length_95%HPD={3.80167600e-02,5.43496800e-02}],11[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:6.523289e-02[&length_mean=6.56442368e-02,length_median=6.52328900e-02,length_95%HPD={5.73263000e-02,7.31237900e-02}])[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:1.894015e-02[&length_mean=1.89293994e-02,length_median=1.89401500e-02,length_95%HPD={1.36761900e-02,2.27847900e-02}])[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:6.711367e-03[&length_mean=7.26324518e-03,length_median=6.71136700e-03,length_95%HPD={4.02061500e-03,1.04333800e-02}])[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:1.126808e-02[&length_mean=1.10690108e-02,length_median=1.12680800e-02,length_95%HPD={7.77257900e-03,1.52072500e-02}])[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:3.109043e-02[&length_mean=3.16918276e-02,length_median=3.10904300e-02,length_95%HPD={2.46074300e-02,3.84590000e-02}],3[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:9.394743e-02[&length_mean=9.28133738e-02,length_median=9.39474300e-02,length_95%HPD={8.29263100e-02,1.02087100e-01}])[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:1.480231e-02[&length_mean=1.47342397e-02,length_median=1.48023100e-02,length_95%HPD={1.14382400e-02,1.87087700e-02}],10[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:3.465900e-02[&length_mean=3.44781638e-02,length_median=3.46590000e-02,length_95%HPD={2.77199000e-02,3.90480900e-02}])[&prob=9.41176471e-01,prob_stddev=8.31890331e-02,prob_range={8.82352941e-01,1.00000000e+00},prob(percent)="94",prob+-sd="94+-8"]:1.920188e-03[&length_mean=2.30552003e-03,length_median=1.92018800e-03,length_95%HPD={7.27882000e-04,3.68785700e-03}],8[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:2.568806e-02[&length_mean=2.59187574e-02,length_median=2.56880600e-02,length_95%HPD={2.15827000e-02,2.94196700e-02}])[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:4.891049e-02[&length_mean=4.81608729e-02,length_median=4.89104900e-02,length_95%HPD={4.04868600e-02,5.55982800e-02}],(6[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:5.403450e-04[&length_mean=6.25692220e-04,length_median=5.40345000e-04,length_95%HPD={3.65281900e-05,1.47964200e-03}],12[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:2.831855e-02[&length_mean=2.87012465e-02,length_median=2.83185500e-02,length_95%HPD={2.52717800e-02,3.36914100e-02}])[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:1.587571e-02[&length_mean=1.60931885e-02,length_median=1.58757100e-02,length_95%HPD={1.19296500e-02,1.99041900e-02}])[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:5.287911e-03[&length_mean=5.51374997e-03,length_median=5.28791100e-03,length_95%HPD={2.95650500e-03,7.71555700e-03}]);"
```
#### carnivore tree

Here the actual data tree:

```bash
(((Fagopyrum_esculentum:99.295536,((Nepenthes_khasiana:78.324329,Ancistrocladus_tectorius:78.324329):4.66512,(Drosera_regia:63.040033,(Aldrovanda_vesiculosa:54.959634,Dionaea_muscipula:54.959634):8.0804):19.949415):16.306087):24.438701,((Acacia_ligulata:115.785565,(Cephalotus_follicularis:91.346395,Oxalis_corniculata:91.346395):24.43917)mrcaott2ott371:2.793039,Eucalyptus_diversicolor:118.578604)mrcaott2ott96:5.155633)Pentapetalae:12.17795,Liriodendron_tulipifera:135.912187)Mesangiospermae;
```