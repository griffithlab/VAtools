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

Snvs and indels are usually run separately through bam-readcount because indels
require to be run in insertion-centric mode (``-i`` option). When using the
``-vcf-readcount-annotator``, the
``--variant-type`` option can then be used to annotate your VCF with those two
files separately. For example, you could run the ``vcf-readcount`` annotator
once with the ``--variant-type snv`` option to run in snv-only mode using the snv
bam-readcount output file and then annotate the output file from that step with
indel information by using the ``--variant-type indel`` option and the
indel bam-readcount output file. This is generally recommended because the
``all`` option in conjunction with a concatenated
bam-readcount output file (containing both snvs and indels) will not be able to handle
cases with a snv and indel at the same position. This situation results in
duplicated bam-readcount entries in the concatenated file, one from the snv
and one from the indel, that might contain conflicting information that can't
be resolved by the ``vcf-readcount-anntator``.

**Example commands for running the vcf-readcount-annotator with snvs and indels
separately**

.. code-block:: none

   vcf-readcount-annotator <input_vcf> <snv_bam_readcount_file> <DNA|RNA> \
   -s <sample_name> -t snv -o <snv_annotated_vcf>

   vcf-readcount-annotator <snv_annotated_vcf> <indel_bam_readcount_file> <DNA|RNA> \
   -s <sample_name> -t indel -o <annotated_vcf>

Usage
-----

.. program-output:: vcf-readcount-annotator -h
