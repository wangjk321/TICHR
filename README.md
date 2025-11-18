## TICHR: Transcriptional regulation analysis by integrating epigenome, 3D genome and transcriptome


<p align="center">
  <a href="https://pypi.org/project/tichr/">
    <img src="https://img.shields.io/pypi/v/tichr?color=4C8BF5&label=PyPI%20Version">
  </a>
  <a href="https://tichr.readthedocs.io/en/latest/">
    <img src="https://readthedocs.org/projects/tichr/badge/?version=latest&color=4C8BF5">
  </a>
  <a href="https://pepy.tech/project/tichr">
    <img src="https://static.pepy.tech/badge/tichr" alt="Downloads">
  </a>
  <img src="https://img.shields.io/badge/Python-3.8%2B-3A6EDB">
  <a href="https://doi.org/10.5281/zenodo.17629590">
    <img src="https://zenodo.org/badge/DOI/10.5281/zenodo.17629590.svg">
  </a>
  <img src="https://img.shields.io/pypi/l/tichr?color=39CC8F&label=License">
  <img src="https://img.shields.io/github/last-commit/wangjk321/tichr?color=4C8BF5">
  <img src="https://img.shields.io/github/repo-size/wangjk321/tichr?color=3A6EDB">
  <a href="https://github.com/wangjk321/tichr/stargazers">
    <img src="https://img.shields.io/github/stars/wangjk321/tichr?style=social">
  </a>
  <!--
  <img src="https://img.shields.io/github/watchers/wangjk321/tichr?style=social">
  <img src="https://img.shields.io/github/followers/wangjk321?style=social">
  -->
</p>




## Introduction

<p align="center">
  <img src="image/logo.png" width="400">
</p>

Deciphering transcriptional regulation across multiple omics layers is essential for understanding cellular processes and disease mechanisms, yet remains challenging due to the limited direct connections among epigenomic, 3D genomic, and transcriptomic data. Here, we present TICHR, a scalable multiomics integration framework that quantifies both site-to-gene and gene-level regulation through diverse weighting strategies, enabling genome-wide characterization of regulatory programs. 


## Functions

TICHR offered sophisticated downstream functions to address diverse research questions for transcriptional regulation, including enhancer prediction, attribution of transcriptional changes, assessment of regulation–transcription concordance, identification of context-specific regulations, cross-sample analysis of large-scale and single-cell data, and characterization of time-series transcriptional dynamics. These multi-task capabilities establish TICHR as an efficient framework for studying complex transcriptional mechanisms.

<p align="center">
  <img src="image/overall.png" width="800">
</p>


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
from tichr import *
```

## Document and tutorial 

You can find a complete tutorial in https://tichr.readthedocs.io/

Although the individual functions of TICHR were described separately, they could be streamlined into a typical workflow:

<p align="center">
  <img src="image/workflow.png" width="500">
</p>

Deposited Data for the TICHR Project: https://zenodo.org/records/17629590

## Contact information
You can open an issue in this repository, or contact the developer (wangjk321@gmail.com)

## Citations
In preparation, 2025.


