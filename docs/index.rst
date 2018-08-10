.. VCF Annotation Tools documentation master file, created by
   sphinx-quickstart on Wed Aug  8 09:54:20 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

VCF Annotation Tools
====================

VCF Annotation Tools is a python package that includes several tools to
annotate VCF files with data from other tools.

**vcf-readcount-annotator**
    A tool that will add the data from bam-readcount files to the VCF sample
    column.

**vcf-expression-annotator**
    A tool that will add the data from several expression tools' output files
    to the VCF INFO column. Supported tools are StringTie, Kallisto,
    and Cufflinks. There also is a ``custom`` option to annotate with data
    from any tab-delimited file.

**vcf-info-annotator**
    A tool that will add data from a tab-delimited file to any user-specified
    field in the VCF INFO column.

.. toctree::
   :maxdepth: 1

   vcf_readcount_annotator
   vcf_expression_annotator
   vcf_info_annotator
   install
   license
   contact
