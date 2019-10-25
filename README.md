# VCF Annotation Tools (VAtools)

VCF Annotation Tools is a python package that includes several tools to annotate VCF files with data from other tools.

**vcf-readcount-annotator**

A tool that will add the data from bam-readcount files to the VCF sample column.

**vcf-expression-annotator**

A tool that will add the data from several expression tools’ output files to the VCF INFO column. Supported tools are StringTie, Kallisto, and Cufflinks. There also is a `custom` option to annotate with data from any tab-delimited file.

**vcf-info-annotator**

A tool that will add data from a tab-delimited file to any user-specified field in the VCF INFO column.

**vcf-genotype-annotator**

A tool to add a new sample to an existing VCF file.

**vep-annotation-reporter**

A tool to create a tab-delimited (TSV) file of variants in a VCF and their VEP annotations.

**transform-split-values**

A tool that extracts and manipulates values from existing sample fields and outputs the results to a TSV file.

## Documentation

Please see [vatools.org](http://vatools.org) for the full documentation.

## Install

`pip install vatools`
