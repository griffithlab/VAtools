Reference Transcript Mismatch Reporter
======================================

This tool can be used to identify variants where the reference genome build
doesn't match the Ensembl reference transcript used by VEP for variant consequence
annotations. In these cases, the REF nucleotide(s) at a variant position will differ
from the Ensembl transcript nucleotide(s) at the corresponding mutation
position. Any resulting amino acid change predictions found in the ``Amino_acids``
field of the VEP CSQ annotation will then be different from the translated Ensembl transcript
amino acids at that position.

This will lead to errors in some downstream tools, e.g. pVACseq, which rely on
the ``Amino_acids`` field as well as the translated Ensembl transcript peptide
sequence - as reported by the ``Wildtype`` plugin to make predictions about
the impact of the mutation on the transcript peptide sequence. Since those two
fields will be in disagreement in such cases as described above, pVACseq
cannot make predictions on such variants.

Such errors might occur in a small number of variants if there are only minor
differences between the reference used and the Ensembl transcripts but they
might also be more widespread, for example, if users aligned to GRCh37 but
used a GRCh38 VEP cachce.

The input VCF needs to be annotated by VEP and requires annotation with the
``Wildtype`` VEP plugin available as part of `pVACtools <https://pvactools.readthedocs.io/en/latest/pvacseq/input_file_prep/vep.html#installing-vep>`_.

This tool will report on the number of variants and transcripts in a VCF that
are affected by this issuei and output this information to stdout. It will
write a ``.mismatch.tsv`` file next to the VCF that provides further details
on the problematic variants.

This tool also allows the user to either soft-filter or hard-filter the VCF
using the ``--filter [soft|hard]`` parameter. Soft-filtering will tag the
problematic variants with a custom VCF FILTER ``CSQ_MISMACH`` while hard-filtering
will produce a new VCF that has these variants removed. When using a filter,
the output VCF will be written to a ``filtered.vcf`` file next to
your input VCF file. You can set a different output file using the
``--output-vcf`` parameter.

Usage
-----

.. program-output:: ref-transcript-mismatch-reporter -h

Example Command
---------------

.. code-block:: none

   csq-mismatch-report input.vcf --filter soft
