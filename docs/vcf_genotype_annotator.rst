VCF Genotype Annotator
======================

The VCF Genotype Annotator will take an existing VCF file and add a new sample
to it. The sample name is set via the second positional argument.
The sample's GT field is pre-populated with a default value given by
the third positional argument. Options are ``0/1``, ``1/1``, ``0/0``, or
``.``.

It can also be used to add a GT field to an existing sample, e.g. for VCFs
created by Strelka which does not output a GT field for its calls.

By default the output VCF will be written to a ``.genotype.vcf`` file next to
your input VCF file. You can set a different output file using the
``--output-vcf`` parameter.

Usage
-----

.. program-output:: vcf-genotype-annotator -h
