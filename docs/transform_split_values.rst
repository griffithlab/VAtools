Transform Split Values
======================

The Transform Split Values tool extracts and manipulates values from existing
sample fields in a VCF and outputs the results to a TSV file. The field to
manipulate is chosen via the second positional argument.

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

If your VCF is a multi-sample VCF, you have to pick one of the sample in
your VCF by setting the ``--sample-name`` option. This is the sample that the
readcounts will be written for.

By default the output TSV will be written to a ``.tsv`` file next to
your input VCF file. You can set a different output file using the
``--output-tsv`` parameter.

Usage
-----

.. program-output:: transform-split-values -h
