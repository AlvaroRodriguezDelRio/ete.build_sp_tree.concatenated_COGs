# ete.build_sp_tree.concatenated_COGs

Script for creating a phylogenetic tree from species proteomes.

For running the script, a miniconda environment with ete3 installed has to be installed. Directory with the paths to the COGs has to be
specified in the script.

###  Usage

```
python ete_build.COG_concat.py dir_input_faa [-n dir_input_fna] [-e evalue_threshold]
```

- `dir_input_faa`             Directory with all protein prediction fasta files of the different species (one file per proteome, "faa" extension)
- `-n dir_input_fna`          Directory with all gene prediction files of the different species (one file per proteome, "fna" extension)
- `-e evalue_threshold`       Minimum e-value for considering a hmmsearch vs the COGs as significant (default = 0.001)

TO-DO:
- only one proteome with a hit for a COG > ete throws error > avoid
- if hit to COG0086 > COG0085 > don't include COG0085
- tab out for parsing hmm search out


DUDAS:
- More than one hit for one species in a COG? (taking the best now)
- One gene hit to more than one COG? (now renaming in ete command)
