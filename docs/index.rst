.. VAtools  documentation master file, created by
   sphinx-quickstart on Wed Aug  8 09:54:20 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

VCF Annotation Tools (VAtools)
==============================

VAtools (VCF Annotation Tools) is a python package for easy manipulation of genomic data stored in the common VCF format

:doc:`vcf-readcount-annotator <vcf_readcount_annotator>`
    A tool that will add the data from `bam-readcount <https://github.com/genome/bam-readcount>`_ files to the VCF sample
    column. Writes depth, allele counts, and VAFs; optionally also writes
    detailed per-read quality metrics (mapping quality, base quality, strand
    counts, and more) as additional FORMAT fields.

:doc:`vcf-expression-annotator <vcf_expression_annotator>`
    A tool that will add the data from several expression tools' output files
    to the VCF FORMAT column, on a per-sample basis (use ``-s`` to select the
    sample for multi-sample VCFs). Directly supports outputs from `StringTie <https://github.com/gpertea/stringtie>`_,
    `Kallisto <https://kallisto.readthedocs.io/en/latest/>`_, and
    `Cufflinks <https://github.com/cole-trapnell-lab/cufflinks>`_.
    There also is a ``custom`` option to annotate with data
    from any tab-delimited file.

:doc:`vcf-info-annotator <vcf_info_annotator>`
    A general-purpose tool that will add data from a tab-delimited file into VCF INFO fields.
    Supports mapping multiple TSV columns to multiple INFO fields in a single
    pass.

:doc:`vcf-genotype-annotator <vcf_genotype_annotator>`
    A tool to add a new sample to an existing VCF file or fill in the GT field
    for an existing sample in a VCF. Fills a need for genotype manipulation in
    VCFs that don't contain one, which can cause errors in downstream tools.

:doc:`vep-annotation-reporter <vep_annotation_reporter>`
    A tool to parse the complex `VEP <https://www.ensembl.org/vep>`_-added CSQ field from a VCF and create a tab-delimited (TSV) file of variants and their VEP annotations.

:doc:`ref-transcript-mismatch-reporter <ref_transcript_mismatch_reporter>`
    A tool to identify variants in a VCF where the reference genome used to
    align and call variants doesn't match the Ensembl reference transcript
    used by VEP for variant consequence annotations.

:doc:`transform-split-values <transform_split_values>`
    A tool that extracts and manipulates values from existing sample fields
    and outputs the results to a TSV file.

.. toctree::
   :maxdepth: 1

   vcf_readcount_annotator
   vcf_expression_annotator
   vcf_info_annotator
   vcf_genotype_annotator
   vep_annotation_reporter
   ref_transcript_mismatch_reporter
   transform_split_values
   install
   license
   contact
