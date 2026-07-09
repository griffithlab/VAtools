VCF Readcount Annotator
=======================

The VCF Readcount Annotator takes an output file from
`bam-readcount <https://github.com/genome/bam-readcount>`_
and adds its data to your VCF. It supports both DNA and RNA readcounts.

Usage
-----

.. program-output:: vcf-readcount-annotator -h

## Details

DNA readcounts are identified by specifying ``DNA`` in the list of
positional arguments. Depth, allele counts, and VAFs are then written to the
DP, AD, and AF fields, respectively. Forward and reverse strand allele counts
are written to the ADF and ADR fields, respectively.

RNA readcounts are identified by specifying ``RNA`` in the list of positional
arguments. Depth, allele counts, and VAFs are then written to the RDP, RAD,
and RAF fields, respectively. Forward and reverse strand allele counts
are written to the RADF and RADR fields, respectively.

If your VCF is a multi-sample VCF, you must specify a sample using the
``--sample-name`` option. The readcounts will be written for that sample only.

By default the output VCF will be written to a ``.readcount.vcf`` file next to
your input VCF file. You can set a different output file using the
``--output-vcf`` parameter.

Gzip-compressed bam-readcount files are also accepted.

### SNVs and Indels

SNVs and indels are usually run separately through bam-readcount because indels
require insertion-centric mode (the ``-i`` option in bam-readcount). The
``--variant-type`` option can then be used to annotate your VCF with each
output file separately. For example, run the annotator once with
``--variant-type snv`` using the SNV bam-readcount file, then annotate the
output of that step with ``--variant-type indel`` using the indel file. This
two-pass approach is generally recommended because the ``all`` option used
with a concatenated bam-readcount file cannot handle the case where a SNV and
an indel exist at the same position — the duplicate entries cannot be
resolved cleanly.

### Extra bam-readcount fields

bam-readcount records per-base quality statistics beyond simple allele
counts. These additional metrics can be written to the VCF as optional
FORMAT fields using the flags below.

Use ``--all-fields`` (``-a``) as a convenience flag to enable all extra
fields at once, or select individual fields using the flags in the table
below.

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Flag
     - FORMAT tag
     - Description
   * - ``-a`` / ``--all-fields``
     - *(all below)*
     - Enable all extra bam-readcount fields
   * - ``-q`` / ``--avg-mapping-quality``
     - ``VAMQ``
     - Avg mapping quality of variant-supporting reads
   * - ``-b`` / ``--avg-basequality``
     - ``VABQ``
     - Avg base quality of variant-supporting reads
   * - ``-e`` / ``--avg-se-mapping-quality``
     - ``VASEMQ``
     - Avg SE mapping quality of variant-supporting reads
   * - ``-r`` / ``--strand-counts``
     - ``ADF``, ``ADR``
     - Forward/reverse strand read counts. In DNA mode, ADF and ADR are
       already written by default; this flag is a no-op with a warning.
   * - ``-f`` / ``--avg-pos-fraction``
     - ``VAPF``
     - Avg position of variant reads as a fraction of read length
   * - ``-m`` / ``--avg-mismatches``
     - ``VAMF``
     - Avg mismatches per variant-supporting read (as a fraction)
   * - ``-k`` / ``--sum-mismatch-qual``
     - ``VAMQS``
     - Avg sum of mismatch base qualities for variant reads
   * - ``-2`` / ``--num-q2-reads``
     - ``VAQ2``
     - Number of variant-supporting reads containing a Q2 base
   * - ``-d`` / ``--avg-q2-distance``
     - ``VAQD``
     - Avg distance to Q2 start in Q2-containing reads
   * - ``-c`` / ``--avg-clipped-length``
     - ``VACL``
     - Avg clipped read length for variant-supporting reads
   * - ``-3`` / ``--avg-3p-distance``
     - ``VA3P``
     - Avg distance to effective 3' end for variant reads


Examples
--------

Download the example data used below:

.. code-block:: none

   curl -LO https://vatools.readthedocs.io/en/latest/_static/vatools-examples.tar.gz
   tar xzf vatools-examples.tar.gz
   cd vatools-examples

**1. Basic DNA annotation**

.. code-block:: none

   vcf-readcount-annotator sample.vcf sample.snv.bam_readcount DNA \
     -o sample.dna.readcount.vcf

Adds ``DP``, ``AD``, ``ADF``, ``ADR``, ``AF`` to the FORMAT column for the single sample.

**2. RNA annotation** (same VCF, same readcounts, different data type)

.. code-block:: none

   vcf-readcount-annotator sample.vcf sample.snv.bam_readcount RNA \
     -o sample.rna.readcount.vcf

Adds ``RDP``, ``RAD``, ``RADF``, ``RADR``, ``RAF`` instead (R for RNA), so DNA and RNA readcounts can coexist on one VCF.

**3. Multi-sample VCF, selecting a sample**

.. code-block:: none

   vcf-readcount-annotator sample.multi_sample.vcf sample.snv.bam_readcount DNA \
     -s H_NJ-HCC1395-HCC1395 -o sample.multi.readcount.vcf

Demonstrates ``-s``/``--sample-name``, which is required because the input has 2 samples. Only the selected sample's FORMAT values are filled in.

**4. SNV/indel two-pass workflow with** ``--variant-type``

.. code-block:: none

   vcf-readcount-annotator sample.snvs_and_indels.vcf sample.snv.bam_readcount DNA \
     -t snv -o snv_annotated.vcf

   vcf-readcount-annotator snv_annotated.vcf sample.indel.bam_readcount DNA \
     -t indel -o annotated.vcf

``sample.snvs_and_indels.vcf`` has both SNVs and indels; the first pass annotates only the SNVs. The second command adds the indel readcounts. This two-pass approach is recommended because bam-readcount needs ``-i`` (insertion-centric mode) for indels and is normally run separately for SNVs vs indels.

**5. Extra per-read quality fields**

.. code-block:: none

   vcf-readcount-annotator sample.vcf sample.snv.bam_readcount DNA \
     -a -o sample.all_fields.readcount.vcf

``--all-fields`` adds ``VAMQ``, ``VABQ``, ``VASEMQ``, ``VAPF``, ``VAMF``, ``VAMQS``, ``VAQ2``, ``VAQD``, ``VACL``, ``VA3P``. These fields can be individually added by using combinations of the short flags (``-q -b -e -f -m -k -2 -d -c``).
