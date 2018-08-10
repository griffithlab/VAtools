VCF Info Annotator
==================

The VCF Info Annotator will add data form a tab-delimited (TSV) file to a
VCF's INFO column.

The TSV file needs to contain three columns in the following order:

- chromsome
- position
- the value for your field at that position

To define the new INFO field you need to specify a info field name in the
positional parameters. This term will be used as the ID field in the INFO
header. You will also need to specify a description in quotes that will be
used as the Description field in the INFO header. Lastly, you will need to
specify the format of your data. This can be either be ``Integer``, ``Float``,
``Flag``, ``Character``, or ``String``.

Optional, you can also set the Source and Version fields of the INFO header
using the ``--source`` and ``--version`` parameters, respectively.

By default the output VCF will be written to a ``.info.vcf`` file next to
your input VCF file. You can set a different output file using the
``--output-vcf`` parameter.

Usage
-----

.. program-output:: vcf-info-annotator -h
