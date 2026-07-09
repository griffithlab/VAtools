VCF Info Annotator
==================

The VCF Info Annotator adds data from a tab-delimited (TSV) file to a
VCF's INFO column. It supports annotating one or more INFO fields in a
single run, and can handle any VCF-spec data type.

Usage
-----

.. program-output:: vcf-info-annotator -h



Details
----------------

**Input TSV format**

The TSV file must have a header row. The first two columns must be
chromosome and position (one-based coordinates, used to match rows to VCF records). Any
additional columns can be mapped to VCF INFO fields by name.

Example TSV with two data columns::

   chrom   pos         freq    classification
   chr1    168192360   0.042   benign
   chr1    230456789   0.187   pathogenic

Gzip-compressed TSV files (``.tsv.gz``) are also accepted.

**Defining column mappings**

Use the ``-m`` / ``--column-mappings`` flag to specify how TSV columns
map to VCF INFO fields. Each mapping is a colon-delimited string with
four required fields and two optional fields:

.. code-block:: none

   source_col:info_field:type:description[:source[:version]]

- **source_col** — the column name in the TSV header
- **info_field** — the ID to use for the INFO field in the VCF
- **type** — the VCF data type: ``Integer``, ``Float``, ``Flag``, ``Character``, or ``String``
- **description** — free-text description written to the INFO header line
- **source** *(optional)* — the source database or tool name
- **version** *(optional)* — the source version (requires source to be set)

To annotate multiple INFO fields in one run, separate mappings with a
comma:

.. code-block:: none

   -m "col1:FIELD1:type:description,col2:FIELD2:type:description"

**Overwriting existing fields**

By default, the tool raises an error if the VCF already contains an
INFO field with the same ID as a mapped field. Use ``--overwrite``
(``-w``) to allow writing to existing fields.

``--clear-existing`` extends this behavior: when set, the existing value
is removed from **every** record before annotation, so records that have
no matching TSV entry will have no value for that field rather than
retaining the old one. ``--clear-existing`` requires ``--overwrite``.

**Output**

By default the output VCF is written to a ``.info.vcf`` file next to
your input VCF. Use ``--output-vcf`` to specify a different path.


Examples
--------

Download the example data used below:

.. code-block:: none

   curl -LO https://vatools.readthedocs.io/en/latest/_static/vatools-examples.tar.gz
   tar xzf vatools-examples.tar.gz
   cd vatools-examples

**1. Annotate a single Integer field**

.. code-block:: none

   vcf-info-annotator sample.vcf sample.info_single.tsv \
     -m "value:DEPTH_SCORE:Integer:Placeholder depth-based score" \
     -o sample.info.vcf

Adds ``DEPTH_SCORE=2`` to the 22:18644673 record.

**2. Annotate two fields in one run**

.. code-block:: none

   vcf-info-annotator sample.vcf sample.info_multi.tsv \
     -m "mapping_quality:MQ0:Integer:Mapping quality score,clinvar_classification:CVCLASS:String:ClinVar classification" \
     -o sample.info.vcf

Adds ``MQ0=30;CVCLASS=Pathogenic`` to the same record.

**3. Include source and version in the INFO header**

.. code-block:: none

   vcf-info-annotator sample.vcf sample.info_multi.tsv \
     -m "mapping_quality:MQ0:Integer:Mapping quality score:ClinVar database:2.1" \
     -o sample.info.vcf

Adds ``Source=ClinVar database,Version=2.1`` to the ``MQ0`` header line.

**4. Overwrite an existing field, clearing it from records not in the TSV**

.. code-block:: none

   vcf-info-annotator sample.info.vcf sample.info_multi.tsv \
     -m "mapping_quality:MQ0:Integer:Mapping quality score" \
     -w --clear-existing -o sample.info.overwritten.vcf

Re-runs against the already-annotated output from example 2, which requires ``-w``/``--overwrite``. Adding ``--clear-existing`` means any record *not* in the TSV would lose its existing ``MQ0`` value rather than keep it (not visible here with only 1 record in the file).
