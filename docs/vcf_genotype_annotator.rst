VCF Genotype Annotator
======================

The VCF Genotype Annotator will take an existing VCF file and add a new sample
to it. 

Usage
-----

.. program-output:: vcf-genotype-annotator -h

Details
-------

The sample name is set via the second positional argument.
The sample's GT field is pre-populated with a default value given by
the third positional argument. Options are ``0/1``, ``1/1``, ``0/0``, or
``.``.  This is useful because some tools do not generate sample/genotype fields, and some downstream tools require them. 

It can also be used to add a GT field to an existing sample, e.g. for VCFs
created by Strelka which does not output a GT field for its calls.

By default the output VCF will be written to a ``.genotype.vcf`` file next to
your input VCF file. You can set a different output file using the
``--output-vcf`` parameter.



Examples
--------

Download the example data used below:

.. code-block:: none

   curl -LO https://vatools.readthedocs.io/en/latest/_static/vatools-examples.tar.gz
   tar xzf vatools-examples.tar.gz
   cd vatools-examples

**1. Add a brand-new sample**

.. code-block:: none

   vcf-genotype-annotator sample.vcf NORMAL 0/0 -o sample.genotype.vcf

``sample.vcf`` only has sample ``H_NJ-HCC1395-HCC1395``; this adds a second sample column ``NORMAL`` with ``GT=0/0``.

**2. Fill in GT for an existing sample missing it** (e.g. Strelka output)

.. code-block:: none

   vcf-genotype-annotator sample.no_gt.vcf H_NJ-HCC1395-HCC1395 0/1 \
     -o sample.genotype.vcf

``sample.no_gt.vcf``'s FORMAT column has no ``GT``; this inserts one for the existing sample rather than adding a new sample column.
