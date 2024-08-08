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
The Newick string NCBI
```bash
(((NC_042597:0.026590858524101986,NC_051971:0.030032337254327683)Node_3:0.006881205739188441,NC_026134:0.03977029703373102)Node_8:0.0008734026701547552,((((NC_041258:0.02057999146509946,NC_051455:0.01903169980006732)Node_4:0.0032818007568811494,((NC_035415:0.014364544390109256,NC_035417:0.03342188333322965)Node_1:0.0012238676593345453,Av.cp.final:0.024342589470773012)Node_2:0.028789123730448153)Node_5:0.0024728919891094762,NC_010776:0.04009233544831099)Node_6:0.007594193081470097,NC_008326:0.03532286573462001)Node_7:0.0018929424416855804,NC_022402:0.029059106129023043)Node_9:1.0;
```

HERE the actual example that we ran
```bash
# Directory containing the input file
sourcedir="/mnt/s-ws/everyone/annotation/"

# Newick tree string placeholder
tree="(((NC_042597:0.026590858524101986,NC_051971:0.030032337254327683)Node_3:0.006881205739188441,NC_026134:0.03977029703373102)Node_8:0.0008734026701547552,((((NC_041258:0.02057999146509946,NC_051455:0.01903169980006732)Node_4:0.0032818007568811494,((NC_035415:0.014364544390109256,NC_035417:0.03342188333322965)Node_1:0.0012238676593345453,Av.cp.final:0.024342589470773012)Node_2:0.028789123730448153)Node_5:0.0024728919891094762,NC_010776:0.04009233544831099)Node_6:0.007594193081470097,NC_008326:0.03532286573462001)Node_7:0.0018929424416855804,NC_022402:0.029059106129023043)Node_9:1.0;"

# Process the TSV file to replace accession numbers with tax IDs
while IFS=$'\t' read -r latin name taxid accession carnivorous australian
do
    # Replace accession with tax ID in the tree string
    tree=$(sed "s/$accession/$taxid/g" <<< "$tree")
done < <(tail -n +2 "$sourcedir/genome_metadata.tsv")

# Print the modified Newick tree
echo "$tree"
```


##### Mash Newick with taxaID
```bash
(((3775:0.026590858524101986,212256:0.030032337254327683)Node_3:0.006881205739188441,138025:0.03977029703373102)Node_8:0.0008734026701547552,((((714108:0.02057999146509946,122310:0.01903169980006732)Node_4:0.0032818007568811494,((4371:0.014364544390109256,4362:0.03342188333322965)Node_1:0.0012238676593345453,Av.cp.final:0.024342589470773012)Node_2:0.028789123730448153)Node_5:0.0024728919891094762,180217:0.04009233544831099)Node_6:0.007594193081470097,3415:0.03532286573462001)Node_7:0.0018929424416855804,166933:0.029059106129023043)Node_9:1.0;
```


#### IQ-TREE

Here the actual tree:

```bash
# Directory containing the input file
sourcedir="/mnt/s-ws/everyone/annotation/"

# Newick tree string placeholder
tree="((((((Av.cp.final.rpoC2:0.1059226071,NC_035417.rpoC2:0.1593292621)100/100:0.0161163570,NC_035415.rpoC2:0.047504210
0)100/100:0.1441383004,NC_041258.rpoC2:0.0768405933)100/100:0.0073042781,NC_051455.rpoC2:0.1013574707)100/100:0.04
57779817,NC_010776.rpoC2:0.2783070161)100/100:0.0952923125,NC_008326.rpoC2:0.1916193482,(NC_022402.rpoC2:0.1884444
229,(NC_026134.rpoC2:0.2385038359,(NC_042597.rpoC2:0.1377749320,NC_051971.rpoC2:0.1945001615)100/100:0.0537714440)
100/100:0.0255907905)100/100:0.0309114382);"

# Process the TSV file to replace accession numbers with tax IDs
while IFS=$'\t' read -r latin name taxid accession carnivorous australian
do
    # Replace accession with tax ID in the tree string
    tree=$(sed "s/$accession/$taxid/g" <<< "$tree")
done < <(tail -n +2 "$sourcedir/genome_metadata.tsv")

# Print the modified Newick tree
echo "$tree"
```

##### IQtree with taxaID
```bash
((((((Av.cp.final.rpoC2:0.1059226071,4362.rpoC2:0.1593292621)100/100:0.0161163570,4371.rpoC2:0.047504210
0)100/100:0.1441383004,714108.rpoC2:0.0768405933)100/100:0.0073042781,122310.rpoC2:0.1013574707)100/100:0.04
57779817,180217.rpoC2:0.2783070161)100/100:0.0952923125,3415.rpoC2:0.1916193482,(166933.rpoC2:0.1884444
229,(138025.rpoC2:0.2385038359,(3775.rpoC2:0.1377749320,212256.rpoC2:0.1945001615)100/100:0.0537714440)
100/100:0.0255907905)100/100:0.0309114382);
```


#### MrBayes

Here the Bayes tree:

```bash
# Nexus tree not yet annotated
```bash    
  #NEXUS
[ID: 5280131283]
begin taxa;
        dimensions ntax=11;
        taxlabels
                Av.cp.final.rpoC2
                NC_008326.rpoC2
                NC_010776.rpoC2
                NC_022402.rpoC2
                NC_026134.rpoC2
                NC_035415.rpoC2
                NC_035417.rpoC2
                NC_041258.rpoC2
                NC_042597.rpoC2
                NC_051455.rpoC2
                NC_051971.rpoC2
                ;
end;
begin trees;
        translate
                1       Av.cp.final.rpoC2,
                2       NC_008326.rpoC2,
                3       NC_010776.rpoC2,
                4       NC_022402.rpoC2,
                5       NC_026134.rpoC2,
                6       NC_035415.rpoC2,
                7       NC_035417.rpoC2,
                8       NC_041258.rpoC2,
                9       NC_042597.rpoC2,
                10      NC_051455.rpoC2,
                11      NC_051971.rpoC2
                ;
   tree con_50_majrule = [&U] (1[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:3.452907e-02[&length_mean=3.45212778e-02,length_median=3.45290700e-02,length_95%HPD={2.96225100e-02,3.83524500e-02}],7[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:5.288674e-02[&length_mean=5.34552884e-02,length_median=5.28867400e-02,length_95%HPD={4.78573900e-02,5.82505700e-02}],(((((2[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:6.524875e-02[&length_mean=6.46920862e-02,length_median=6.52487500e-02,length_95%HPD={5.65023700e-02,7.24066900e-02}],(4[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:6.467878e-02[&length_mean=6.42238491e-02,length_median=6.46787800e-02,length_95%HPD={5.61046000e-02,6.88557300e-02}],(5[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:8.107971e-02[&length_mean=8.13590269e-02,length_median=8.10797100e-02,length_95%HPD={7.12122100e-02,8.83254200e-02}],(9[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:4.562096e-02[&length_mean=4.50742869e-02,length_median=4.56209600e-02,length_95%HPD={3.86517800e-02,5.12176800e-02}],11[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:6.589138e-02[&length_mean=6.58006497e-02,length_median=6.58913800e-02,length_95%HPD={5.81845000e-02,7.33511300e-02}])[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:1.900464e-02[&length_mean=1.90827156e-02,length_median=1.90046400e-02,length_95%HPD={1.50565100e-02,2.27622300e-02}])[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:7.761626e-03[&length_mean=7.81708409e-03,length_median=7.76162600e-03,length_95%HPD={4.66618300e-03,1.07821900e-02}])[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:1.075370e-02[&length_mean=1.10413548e-02,length_median=1.07537000e-02,length_95%HPD={6.97324600e-03,1.46796500e-02}])[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:3.158395e-02[&length_mean=3.15678328e-02,length_median=3.15839500e-02,length_95%HPD={2.37178600e-02,3.77140000e-02}],3[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:9.407518e-02[&length_mean=9.37498978e-02,length_median=9.40751800e-02,length_95%HPD={8.15330000e-02,1.01788400e-01}])[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:1.474397e-02[&length_mean=1.49450287e-02,length_median=1.47439700e-02,length_95%HPD={1.06209100e-02,1.92830000e-02}],10[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:3.331415e-02[&length_mean=3.36189069e-02,length_median=3.33141500e-02,length_95%HPD={2.90139400e-02,3.88283300e-02}])[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:1.904403e-03[&length_mean=2.35553543e-03,length_median=1.90440300e-03,length_95%HPD={8.54679500e-04,3.93488100e-03}],8[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:2.543118e-02[&length_mean=2.60639181e-02,length_median=2.54311800e-02,length_95%HPD={2.22574900e-02,2.96783500e-02}])[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:4.832307e-02[&length_mean=4.85161122e-02,length_median=4.83230700e-02,length_95%HPD={3.93723300e-02,5.45909000e-02}],6[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:1.561311e-02[&length_mean=1.55678619e-02,length_median=1.56131100e-02,length_95%HPD={1.15494300e-02,1.82530200e-02}])[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:5.952369e-03[&length_mean=5.74763612e-03,length_median=5.95236900e-03,length_95%HPD={3.23652600e-03,8.28774200e-03}]);
end;length_95%HPD={3.23652600e-03,8.28774200e-03}]);
```   

Translated with species in iTOL friendly format
```bash
  tree con_50_majrule = [&U] (Av.cp.final.rpoC2:3.452907e-02,
Dionaea_muscipula:5.288674e-02,
(((((Liriodendron_tulipifera:6.524875e-02,
(Eucalyptus_diversicolor:6.467878e-02,
(Acacia_ligulata:8.107971e-02,
(Cephalotus_follicularis:4.562096e-02,
Oxalis_corniculata:6.589138e-02):1.900464e-02):7.761626e-03):1.075370e-02):3.158395e-02,
Fagopyrum_esculentum:9.407518e-02):1.474397e-02,
Nepenthes_khasiana:3.331415e-02):1.904403e-03,
Ancistrocladus_tectorius:2.543118e-02):4.832307e-02,
Drosera_regia:1.561311e-02):5.952369e-03);
```

#### Here the actual data tree:

```bash
(((Fagopyrum_esculentum:99.295536,((Nepenthes_khasiana:78.324329,Ancistrocladus_tectorius:78.324329):4.66512,(Drosera_regia:63.040033,(Aldrovanda_vesiculosa:54.959634,Dionaea_muscipula:54.959634):8.0804):19.949415):16.306087):24.438701,((Acacia_ligulata:115.785565,(Cephalotus_follicularis:91.346395,Oxalis_corniculata:91.346395):24.43917)mrcaott2ott371:2.793039,Eucalyptus_diversicolor:118.578604)mrcaott2ott96:5.155633)Pentapetalae:12.17795,Liriodendron_tulipifera:135.912187)Mesangiospermae;
```
