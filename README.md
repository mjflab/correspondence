# Introduction
This repository contains code to reproduce the results of the correspondence.

# Files
## Jupyter notebooks
- `jeme-ripple.ipynb`: To run JEME and RIPPLE using RIPPLE data with various parameter settings.
- `jeme.ipynb`: To run JEME and TargetFinder using JEME "random targets" datasets with various parameter settings.
- `motifs prediction targetfinder-results.ipynb`: To obtain the prediction performances of using Wang et al motifs and shuffled motifs.
- `reproducing TF paper-eep.ipynb`: To obtain the cross-validation with shuffling and chromosome-split results of TargetFinder on its EE/P datasets.
- `reproducing TF paper.ipynb`: To obtain the cross-validation with shuffling and chromsome-split results of TargetFinder on its E/P/W datasets.
- `shuffle_motif.ipynb`: To generate shuffled motifs.

## HTML files
- `tomtom-shuffled-motifs.html`: TomTom results using shuffled motifs against JASPAR 2014 Vertebrates motif database.
- `tomtom-wang-motifs.html`: TomTom results using Wang et al motifs against JASPAR 2014 Vertebrates motif database.

## meme files
- `wang_etal_motifs.meme`: The motifs from Wang et al in MEME format.
- `wang_motifs_all_shuffled.meme`: The shuffled Wang et al motifs in MEME format.

## targetfinder-master folder
### Bash scripts
- `count_overlap_eep_same_promoter.sh`: Count the number of extended enhancers overlap with each other at different cutoff thresholds and share the same promoter.
- `count_overlap_eep.sh`: Count the number of extended enhancers overlap with each other at different cutoff thresholds.
- `count_overlap.same_promoter.sh`: Count the number of windows overlap with each other at different cutoff thresholds and share the same promoter.
- `count_overlap.sh`: Count the number of windows overlap with each other at different cutoff thresholds.
- `generate_pairs_index_eep.sh`: Produce the name tuples of samples with overlapping extended-enhancers at different overlapping cutoffs.
- `generate_pairs_index.sh`: Produce the name tuples of sampels with overlapping windows at different overlapping cutoffs.

### Python scripts
- `generate_window3.py`: Generate window regions as Bed3 format.
