## TICHR: Transcriptional regulation analysis by integrating epigenome, 3D genome and transcriptome


## Introduction

<img src="images/logo.png" width="600">

Deciphering transcriptional regulation across multiple omics layers is essential for understanding cellular processes and disease mechanisms, yet remains challenging due to the limited direct connections among epigenomic, 3D genomic, and transcriptomic data. Here, we present TICHR, a scalable multiomics integration framework that quantifies both site-to-gene and gene-level regulation through diverse weighting strategies, enabling genome-wide characterization of regulatory programs. 


## Functions

<img src="images/overall.png" width="600">

TICHR offered sophisticated downstream functions to address diverse research questions for transcriptional regulation, including enhancer prediction, attribution of transcriptional changes, assessment of regulation–transcription concordance, identification of context-specific regulations, cross-sample analysis of large-scale and single-cell data, and characterization of time-series transcriptional dynamics. These multi-task capabilities establish TICHR as an efficient framework for studying complex transcriptional mechanisms.

## Installation

You can install the latest version of TICHR from PyPI using pip:
``` shell
pip install tichr
```

## Usage

There are two ways to use TICHR

1. **Command Line Interface (CLI)** — This is the most straightforward method for most users. After installation, you can check the available commands with:

``` shell
tichr --help
```

2. **Python Module** — Advanced users can import TICHR as a Python package and use its functions programmatically. For example:

``` python
import sys
from tichr import *
```

## Document and tutorial 

You can find a complete tutorial in https://tichr.readthedocs.io/


## Citations
Jiankang Wang, Puxuan Sun, Bo Zhou, Yaqing Liu, Ke-wei Zheng, Jun Wu, Haoping Chen. TICHR: Transcriptional regulation analysis by integrating epigenome, 3D genome and transcriptome. In preparation, 2025.

