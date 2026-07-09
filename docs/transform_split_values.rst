Transform Split Values
======================

The Transform Split Values tool extracts and manipulates values from existing
sample fields in a VCF and outputs the results to a TSV file. 

Usage
-----

.. program-output:: transform-split-values -h

Details
-------

The field in the TSV to manipulate is chosen via the second positional argument.

Supported operations are the following:

- ``ref``: Extract the first value in a R-number field (the reference value).
- ``alt``: Extract the second value in a R-number field (the alt value).
- ``sum``: Calculate the sum of all the numbers in the field.
- ``min``: Calculate the minimum of all the numbers in the field.
- ``max``: Calculate the maximum of all the numbers in the field.
- ``mean``: Calculate the mean of all the numbers in the field.
- ``median``: Calculate the median of all the numbers in the field.
- ``stdev``: Calculate the standard deviation of all the numbers in the field.
- ``ref_ratio``: The first value in a R-number field divided by the sum of all the numbers (the reference ratio).
- ``alt_ratio``: The second value in a R-number field divided by the sum of all the numbers (the alt ratio).

If your VCF is a multi-sample VCF, you have to pick one of the samples in
your VCF by setting the ``--sample-name`` option. This is the sample whose
field values will be extracted.

By default the output TSV will be written to a ``.tsv`` file next to
your input VCF file. You can set a different output file using the
``--output-tsv`` parameter.


Examples
--------

Download the example data used below:

.. code-block:: none

   curl -LO https://vatools.readthedocs.io/en/latest/_static/vatools-examples.tar.gz
   tar xzf vatools-examples.tar.gz
   cd vatools-examples

First, annotate the VCF with depth information from the readcount data:

.. code-block:: none

   vcf-readcount-annotator sample.vcf sample.snv.bam_readcount DNA \
     -o sample.readcount.vcf

This adds the ``AD`` field (allelic depths, ``Number=R``: ref count, alt count), which the examples below extract and summarize.

**1. Basic, single operation**

.. code-block:: none

   transform-split-values sample.readcount.vcf AD mean -o sample.transform.tsv

``AD`` holds the ref and alt depths as a list (``102,5``); ``mean`` averages them, giving ``53.5``.

**2. Multiple operations in one run**

.. code-block:: none

   transform-split-values sample.readcount.vcf AD mean median min max \
     -o sample.transform.tsv

Produces one output column per operation: ``H_NJ-HCC1395-HCC1395-AD-mean=53.5``, ``...-median=53.5``, ``...-min=5``, ``...-max=102``.

**3. Multi-sample VCF, selecting a sample**

.. code-block:: none

   vcf-readcount-annotator sample.multi_sample.vcf sample.snv.bam_readcount DNA \
     -s H_NJ-HCC1395-HCC1395 -o sample.multi.readcount.vcf

   transform-split-values sample.multi.readcount.vcf AD mean \
     -s H_NJ-HCC1395-HCC1395 -o sample.transform.tsv

Demonstrates ``-s``/``--sample-name`` for both the readcount annotator and ``transform-split-values`` on a multi-sample VCF — only the selected sample's depth values are extracted.

**4. Extracting the alt ratio (aka VAF)**

.. code-block:: none

   transform-split-values sample.readcount.vcf AD alt_ratio -o sample.transform.tsv

``alt_ratio`` divides the alt depth by the total depth (``5 / 107 = 0.0467``) — the variant allele fraction. 

**5.** ``-t`` **merge into an existing TSV**

.. code-block:: none

   transform-split-values sample.readcount.vcf AD ref alt \
     -t sample.variant_report.tsv -o sample.merged.tsv

Same merge behavior as ``vep-annotation-reporter -t``, appending computed columns onto the existing rows instead of starting fresh.
