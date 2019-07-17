Commands used:
-----------
-  raxmlHPC-AVX was used for creating the parameters of the model tha epa needs (RAxML_info.file):

```
raxmlHPC-AVX -f e -s all-genes.fasta.phy -t speciestree.tre  -m GTRGAMMAX
```

- epa-ng was run on reference alignments and tree and the query sequence as follows:

```
epa-ng --ref-msa pruned-all-genes.fasta --tree RAxML_result.pruned --query query.fas --outdir OUTPUT_DIR --model RAxML_info.file --redo -T 1 &> epa.err
```
query.fas contains the query sequence. RAxML\_info.file is the info file created by previous command with raxmlHPC\-AVX. The RAxML_result.pruned is the optimized tree with one leaf pruned.

- guppy was used to create trees from placement outputs:

```
guppy tog -o tog_epa.tre epa_result.jplace 
```

- INSTRAL was run using the following command:

```
java -Xmx2G -jar instral.5.13.4.jar -i genetrees.tree -f pruned.speciestree.tre -o instral.out --placement pruned.speciesname --no-scoring -C -T 1 2> instral.err
```

