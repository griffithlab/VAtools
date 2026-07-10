---
title: 'VAtools: a python package for customizable annotation of VCF files'
tags:
  - Python
  - genomics
  - annotation
  - immunotherapy
  - vcf
authors:
  - name: Susanna Kiwala
    orcid: 0000-0003-4378-7328
    equal-contrib: true
    affiliation: "1" 
  - name: Christopher A. Miller
    orcid: 0000-0003-4266-6700
    equal-contrib: true 
    affiliation: "1, 2"
  - name: John E. Garza
    orcid: 0009-0002-2565-0774
    affiliation: "1" 
  - name: Alexander J. Paul
    affiliation: "1" 
  - name:  Jason Walker
    affiliation: "1" 
  - name: Obi L. Griffith
    orcid: 0000-0002-0843-4271
    affiliation: "1, 2" 
  - name: Malachi Griffith
    orcid: 0000-0002-6388-446X
    affiliation: "1, 2" 
affiliations:
 - name: Washington University in St Louis, USA
   index: 1
 - name: Siteman Cancer Center, St Louis, MO
   index: 2
date: 4 July 2026
bibliography: pubs.bib
---

# Summary:

VAtools (VCF Annotation Tools) is a package designed to ease manipulation of genomic data stored in the common VCF format, so that annotations can be added or extracted in a user-friendly manner.

# Statement of need:

The VCF format is the de facto standard for reporting and sharing genomic variant data, but it exists within an ecosystem of bioinformatics tools that output a wide array of often ill-defined formats. In order to store data in a consistent manner and simplify downstream analyses, computational biologists and genomics researchers often need to integrate these results into a VCF file. Doing so with existing toolchains often requires complex multi-step machinations. We built VAtools to simplify some of the most common operations.

# State of the field:

Several python libraries exist that allow programmatic access to a VCF, such as pysam, VCFPy, and PyVCF [@Heger2026pysam; @Holtgrewe2016-nw; @Vermaat2023jamescasbon]. Using these generally involves writing custom wrappers and requires detailed knowledge of both the structure of a VCF and of the input data. The VCF format is complex (with a 50-page specifications document as of version 4.5: https://samtools.github.io/hts-specs/VCFv4.5.pdf), it can be difficult to parse and manipulate, and this represents a substantial hurdle for many users. Other tools, like bcftools, support general command line manipulation, but rely on the same deep knowledge and have a complicated syntax[@Danecek2021-ta]. VAtools was designed to provide a simple interface for annotating VCFs with some of the most common types of genomic data.

# Software design:

VAtools was built with ease of use and modularity as the primary goals. VAtools was motivated by our own difficulty in preparing inputs for variant analyses, and provides tools that simplify a variety of common annotation tasks. Some are designed specifically around common genomic programs: `vcf-readcount-annotator` adds detailed metrics from bam-readcount, while `vcf-expression-annotator` can link variants to specific gene or transcript expression information from Kallisto or Stringtie, among others[@Khanna2022; @Bray2016-fg; @Pertea2015-ka]. In addition, the package provides more general tools that can take arbitrary data with genomic coordinates and incorporate them into a VCF using `vcf-info-annotator` (e.g. dbSNP allele frequencies or ClinVar pathogenicity classifications). It also includes utilities for extracting variant effect annotation data from VCFs into tabular format via `vep-annotation-reporter`, fixing common issues with malformed VCFs with `vcf-genotype-annotator` and for sanity checking of results with `ref-transcript-mismatch-reporter`. 

Creating this set of tools installable with a single call to `pip` fits well with the paradigms of genomic analysis, which favors rapid installation and the ability to prototype and iterate via the command line. These prototypes often mature into more complex genomic pipelines, for which modularity and containerization are critical. Distributing via pip and Bioconda enables easy installation and version management, as well as containerization that is critical to modern bioinformatics workflow engines like Toil, Cromwell, or Nextflow[@Vivian2017-bo; @Voss2017-za; @Di-Tommaso2017-nx].

VAtools was built on top of the `VCFPy` library, selected over two alternatives in part because of its intuitive syntax that eases development and code review. `PyVCF` (or its new `PyVCF3` fork) was not used because it lacks convenient accessors for modifying VCF headers and per-sample genotype fields. Another alternative, `pysam`, has faster I/O, but requires strict checking of `##contig` lines in the header. As many genomic pipelines do not produce such fields, we ultimately determined that using pysam would limit the reach of this tool.  Though not built on the fastest library, VAtools is nonetheless performant, annotating 1 million records with vcf-info-annotator in 47.8 seconds on a Macbook Pro M2 with 16Gb of RAM and using Python 3.14.0. Further optimizations are possible, but unnecessary for the typical use case of this tool, which operates on the scale of individual exomes (thousands of variants) to whole genomes (~3 million variants).

# Research impact statement:

VAtools was originally designed to work in concert with pVACtools [@Hundal2020-qk; @hoang2026], an immunotherapy tool that incorporates many different sources of data (VEP annotations, gene expression from multiple sources, allele specific expression, readcounts). Custom pipelines using VAtools provided the underlying data for neoantigen and immunotherapy studies in pediatric ependymoma, lung cancer, and high-grade serous ovarian cancer [@Miller2018-od; @Nicholas2026-pt; @Garsed2022-eu]. More recently, VAtools components have been incorporated into a pipeline called ImmunoNX [@Singhal2025-ec] that provides an end-to-end workflow for cancer vaccine design. It has also been used to create neoantigen vaccines in clinical trials for triple-negative breast cancer, follicular lymphoma, and glioblastoma[@Zhang2024-wm; @Ramirez2024-dl; @Garfinkle2026-br], as well as at least 10 other clinical trials currently in progress.

Outside of cancer immunotherapy, VAtools has found use in a variety of other contexts, including importing sequencing metrics to enable filtering of genomic variants in studies of leukemia and malignant peripheral nerve sheath tumors [@Abel2023; @Dehner2021-ei] or in  mouse models of skin cancer [@Shea2023-sz]. It was also used in evolutionary studies of the amoeba Dictyostelium discoideum[@Walker2024].

We expect that VAtools will continue to be generally useful to the genomics community, enabling applications ranging from large-scale automated pipelines to small ad-hoc workflows. VAtools is released under an MIT license to allow for broad reuse. The code is available at https://github.com/griffithlab/vatools, with documentation at https://vatools.org. Versioned releases are made available on PyPi (under the vatools package name), GitHub (https://github.com/griffithlab/VAtools/releases), DockerHub (griffithlab/vatools), BioConda, and Zenodo (https://doi.org/10.5281/zenodo.21079855).

# AI usage disclosure:

No AI was used in the first ~8 years of VAtools development, but starting in 2026, we have used Claude Sonnet 4.6 for occasional assistance in adding new features or refactoring. All code (whether generated by humans or AI) has undergone rigorous code review by the authors of the package and each change must pass a comprehensive test suite integrated into the repository.

# Acknowledgements

MG was supported by the National Human Genome Research Institute (NHGRI) of the National Institutes of Health (NIH) under Award Number R00HG007940. SK, MG and OLG were supported by the NIH National Cancer Institute (NCI) under Award Numbers U01CA209936, U01CA231844, U24CA237719, U24CA275783, U24CA305456, and R01CA301054. SK, MG and OLG were supported by NIH National Center for Advancing Translational Sciences (NCATS) under Award Number UL1TR002345. MG was supported by the V Foundation for Cancer Research under Award Number V2018-007. MG and SK were supported by a gift from the Goldberg Family Foundation. MG, OLG, SK and CAM were supported by the NCI under Award Number U01CA248235. CAM was supported by the NCI under Award Number R50CA211782.
