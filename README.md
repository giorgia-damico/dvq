# DVQ - DNA Visualisations and Quick comparisons 
Abstract
We introduce DVQ (DNA Visualisations and Quick comparisons), an open-source Python library for exploring nucleotide sequences using a variety of methods. Understanding DNA sequences intuitively is important for a variety of tasks in biology. DVQ aims to be a one-stop comprehensive library that makes explainable DNA easy for geneticists, researchers, and practitioners who need explanations. For practitioners, the library provides an easy-to-use interface to generate visualisations for their sequences by only writing a few lines of code. In this report, we demonstrate several example use cases across different types of sequences as well as visualisations. 

Simple early version preprint: http://dx.doi.org/10.13140/RG.2.2.19227.89125

## Methods Implemented:
### Visual
- [x] [Persistant Homological Representations](https://american-cse.org/csci2022-ieee/pdfs/CSCI2022-2lPzsUSRQukMlxf8K2x89I/202800b599/202800b599.pdf)
- [x] [ColorSquare](https://match.pmf.kg.ac.rs/electronic_versions/Match68/n2/match68n2_621-637.pdf)
- [x] [C-Curve](https://pubmed.ncbi.nlm.nih.gov/23246806/) - Removed from development plan due it being redundant vs a 2D Line
- [x] [Spider Representation](https://www.researchgate.net/publication/260971259_Spider_Representation_of_DNA_Sequences) - Removed from development plan due to it being found to be difficult to use for large sequences
- [x] [2D Line](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC162336/)
- [x] [Chaos Game Representation](Fwww.sciencedirect.com%2Fscience%2Farticle%2Fpii%2FS2001037021004736&usg=AOvVaw38odDudWfUCAqbc626rD2e&opi=89978449)

### Statistical
- [x] Deng entropy
- [x] KL Divergence plain, computes distros of chunks similar to deng entropy.
- [x] KL Divergence KL div like measure using deng entropy.
- [ ] [KL Divergence](https://pubmed.ncbi.nlm.nih.gov/31981184/)
- [ ] [Perpelxity](https://arxiv.org/pdf/1202.2518.pdf)
- [x] [Entropy](https://pubmed.ncbi.nlm.nih.gov/9344742/)
- [x] K-Mer overlap
- [ ] *fast* K-Mer overlap, (current implementation is too slow or stuck)
- [x] [Wen's Method](https://pubmed.ncbi.nlm.nih.gov/29765099/)
- [ ] JS Divergence
- [ ] Wasserstein Distance

## Overview:

* How to use dvq
```python
from dvq import visual

visual.plot_2d_comparison([seqs_1, seqs_2], ['seq_1', 'seq_2'])

```

![An example graphic comparing dna sequences for the same virus ](Untitled.png "2D Comparison - Same virus")

```python
from dvq import statistical

statistical.similarity_wen([seqs_1, seqs_2])
# 0.99

```




# Docs
## Funtional

This module provides utilities for detecting repeat patterns in DNA sequences, including fixed-length n-gram repeats and tandem repeats of both fixed and variable length.

### Features

- Detect earliest repeated n-grams
- Detect fixed-length tandem repeats
- Detect variable-length tandem repeats (motif discovery)
- Batch/multiprocessing support

### Classes & Functions

### `DNARepeatDetector`

```python
detector = DNARepeatDetector(n=4, min_repeats=2)
```

#### `.detect_repeats(sequence, return_earliest=False)`
- Returns list of repeated n-grams and their positions.
- If `return_earliest=True`, returns the first repeated n-gram and its first two positions.

#### `.detect_tandem_repeats(sequence)`
- Detects adjacent (tandem) repeats of fixed-length `n`.

### `detect_variable_tandem_repeats(sequence, min_n=2, max_n=6)`
- Finds tandem repeats with motif lengths between `min_n` and `max_n`.

---

#### Multiprocessing Tools

All tools use `multiprocessing.Pool` for performance.

##### `repeat_onset_positions_multiprocess(seqs, n, min_repeats, num_cores)`
- Returns earliest repeat onset positions per sequence.

##### `tandem_repeats_multiprocess(seqs, n, min_repeats, num_cores)`
- Returns list of tandem repeats for each sequence.

##### `variable_tandem_repeats_multiprocess(seqs, min_n, max_n, min_repeats, num_cores)`
- Returns variable-length tandem repeats for each sequence.

---

### Example

```python
from functional.repeat_analysis import (
    repeat_onset_positions_multiprocess,
    tandem_repeats_multiprocess,
    variable_tandem_repeats_multiprocess
)

seqs = ["ACGTACGTACGT", "ATCGATCGATCG", "AGCTTTCGAAGCTTTCGAA"]

onsets = repeat_onset_positions_multiprocess(seqs, n=4, min_repeats=2)
tandems = tandem_repeats_multiprocess(seqs, n=4, min_repeats=2)
variable = variable_tandem_repeats_multiprocess(seqs, min_n=2, max_n=6, min_repeats=2)

print(onsets)
print(tandems)
print(variable)
```

---

### Notes

- Designed for biological sequence analysis and model evaluation.
- Place this module in `src/functional/` and import as needed in larger packages.

