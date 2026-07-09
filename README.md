![Test Status](https://github.com/griffithlab/VAtools/actions/workflows/tests.yml/badge.svg)
[![Docs](https://readthedocs.org/projects/pvactools/badge/?version=latest)](http://vatools.readthedocs.io/en/latest/?badge=latest)
<a href="https://pypi.python.org/pypi/vatools/">
    <img src="https://img.shields.io/pypi/v/vatools.svg?maxAge=1000" alt="PyPI" />
</a>


# VCF Annotation Tools (VAtools)

VAtools (VCF Annotation Tools) is a python package for easy manipulation of genomic data stored in the common VCF format

**vcf-readcount-annotator**

A tool that will add the data from bam-readcount files to the VCF sample column.

**vcf-expression-annotator**

A tool that will add gene expression data to the VCF FORMAT column on a per-sample basis. StringTie, Kallisto, and Cufflinks outputs are supported natively, or there is a `custom` option to annotate with data from any tab-delimited file.

**vcf-info-annotator**

A general-purpose tool that will add data from a tab-delimited file to any user-specified field in the VCF INFO column.

**vcf-genotype-annotator**

A tool to add a new sample to an existing VCF file.

**vep-annotation-reporter**

A tool to extract human-readable variant annotation data (in TSV format) from a VEP-annotated VCF.

**ref-transcript-mismatch-reporter**

A tool to identify variants in a VCF where the reference genome (used to
align and call variants) doesn't match the Ensembl reference transcript
(used by VEP for variant consequence annotations).

**transform-split-values**

A tool that extracts and manipulates values from existing sample fields and outputs the results to a TSV file.

## Documentation

Please see [vatools.org](http://vatools.org) for the full documentation.

## Install

`pip install vatools`

## Container images

VAtools is available as a Docker Image at <a href="https://hub.docker.com/r/griffithlab/vatools/">DockerHub griffithlab/vatools</a>.

## Stable release with DOI

[![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.21079855-blue)](https://doi.org/10.5281/zenodo.21079855)
