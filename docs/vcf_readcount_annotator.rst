VCF Readcount Annotator
=======================

The VCF Readcount Annotator will take an output file from
`bam-readcount <https://github.com/genome/bam-readcount>`_
and add its data to your VCF. It supports both DNA and RNA readcounts.

DNA readcounts are identified by specifying ``DNA`` in the list of
positional arguments. Depth, allele counts, and VAFs are then written to the
DP, AD, and AF fields, respectively.

RNA readcounts are identified by specifying ``RNA`` in the list of positional
arguments. Depth, allele counts, and VAFs are then written tot he RDP, RAD,
and RAF fields, respectively.

If your VCF is a multi-sample VCF, you have to pick one of the sample in
your VCF by setting the ``--sample-name`` option. This is the sample that the
readcounts will be written for.

By default the output VCF will be written to a ``.readcount.vcf`` file next to
your input VCF file. You can set a different output file using the
``--output-vcf`` parameter.

Usage
-----

.. program-output:: vcf-readcount-annotator -h
