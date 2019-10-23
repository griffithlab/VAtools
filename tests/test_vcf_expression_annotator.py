import unittest
import sys
import os
import py_compile
from vatools import vcf_expression_annotator
import tempfile
from filecmp import cmp
import logging
from testfixtures import LogCapture, StringComparison as S

class VcfExpressionEncoderTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        base_dir          = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
        cls.executable    = os.path.join(base_dir, 'vatools', 'vcf_expression_annotator.py')
        cls.test_data_dir = os.path.join(base_dir, 'tests', 'test_data')

    def test_source_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_error_custom_format_id_column_not_set(self):
        with self.assertRaises(Exception) as context:
            command = [
                os.path.join(self.test_data_dir, 'input.vcf'),
                os.path.join(self.test_data_dir, 'genes.fpkm_tracking'),
                'custom',
                'gene',
                '-e', 'FPKM',
            ]
            vcf_expression_annotator.main(command)
        self.assertTrue('--id-column is not set. This is required when using the `custom` format.' == str(context.exception))

    def test_error_custom_format_expression_column_not_set(self):
        with self.assertRaises(Exception) as context:
            command = [
                os.path.join(self.test_data_dir, 'input.vcf'),
                os.path.join(self.test_data_dir, 'genes.fpkm_tracking'),
                'custom',
                'gene',
                '-i', 'tracking_id',
            ]
            vcf_expression_annotator.main(command)
        self.assertTrue('--expression-column is not set. This is required when using the `custom` format.' == str(context.exception))

    def test_error_more_than_one_sample_without_sample_name(self):
        with self.assertRaises(Exception) as context:
            command = [
                os.path.join(self.test_data_dir, 'multiple_samples.vcf'),
                os.path.join(self.test_data_dir, 'genes.fpkm_tracking'),
                'cufflinks',
                'gene',
            ]
            vcf_expression_annotator.main(command)
        self.assertTrue('contains more than one sample. Please use the -s option to specify which sample to annotate.' in str(context.exception))

    def test_error_more_than_one_sample_with_wrong_sample_name(self):
        with self.assertRaises(Exception) as context:
            command = [
                os.path.join(self.test_data_dir, 'multiple_samples.vcf'),
                os.path.join(self.test_data_dir, 'genes.fpkm_tracking'),
                'cufflinks',
                'gene',
                '-s', 'nonexistent_sample',
            ]
            vcf_expression_annotator.main(command)
        self.assertTrue('does not contain a sample column for sample nonexistent_sample.' in str(context.exception))

    def test_error_no_csq(self):
        with self.assertRaises(Exception) as context:
            command = [
                os.path.join(self.test_data_dir, 'no_csq.vcf'),
                os.path.join(self.test_data_dir, 'genes.fpkm_tracking'),
                'cufflinks',
                'gene',
            ]
            vcf_expression_annotator.main(command)
        self.assertTrue('is not VEP-annotated. Please annotate the VCF with VEP before running this tool.' in str(context.exception))

    def test_error_already_TX_annotated(self):
        with self.assertRaises(Exception) as context:
            command = [
                os.path.join(self.test_data_dir, 'input.tx.vcf'),
                os.path.join(self.test_data_dir, 'isoforms.fpkm_tracking'),
                'cufflinks',
                'transcript',
            ]
            vcf_expression_annotator.main(command)
        self.assertTrue('is already transcript expression annotated. TX format header already exists.' in str(context.exception))

    def test_error_already_GX_annotated(self):
        with self.assertRaises(Exception) as context:
            command = [
                os.path.join(self.test_data_dir, 'input.gx.vcf'),
                os.path.join(self.test_data_dir, 'genes.fpkm_tracking'),
                'cufflinks',
                'gene',
            ]
            vcf_expression_annotator.main(command)
        self.assertTrue('is already gene expression annotated. GX format header already exists.' in str(context.exception))

    def test_error_id_column_nonexistent_in_file(self):
        with self.assertRaises(Exception) as context:
            command = [
                os.path.join(self.test_data_dir, 'input.vcf'),
                os.path.join(self.test_data_dir, 'genes.fpkm_tracking'),
                'custom',
                'gene',
                '-i', 'nonexistent_column',
                '-e', 'FPKM',
            ]
            vcf_expression_annotator.main(command)
        self.assertTrue('ERROR: id_column header nonexistent_column does not exist in expression_file' in str(context.exception))

    def test_error_expression_column_nonexistent_in_file(self):
        with self.assertRaises(Exception) as context:
            command = [
                os.path.join(self.test_data_dir, 'input.vcf'),
                os.path.join(self.test_data_dir, 'genes.fpkm_tracking'),
                'custom',
                'gene',
                '-i', 'tracking_id',
                '-e', 'nonexistent_column',
            ]
            vcf_expression_annotator.main(command)
        self.assertTrue('ERROR: expression_column header nonexistent_column does not exist in expression_file' in str(context.exception))

    def test_warning_no_csq_for_variants(self):
        logging.disable(logging.NOTSET)
        with LogCapture() as l:
            temp_path = tempfile.TemporaryDirectory()
            os.symlink(os.path.join(self.test_data_dir, 'input.no_csq.vcf'), os.path.join(temp_path.name, 'input.vcf'))
            command = [
                os.path.join(temp_path.name, 'input.vcf'),
                os.path.join(self.test_data_dir, 'genes.fpkm_tracking'),
                'cufflinks',
                'gene',
            ]
            vcf_expression_annotator.main(command)
            temp_path.cleanup()
            l.check_present(('root', 'WARNING', S("Variant is missing VEP annotation. INFO column doesn't contain CSQ field for variant")))

    def test_skip_variant_without_gene_in_csq(self):
        temp_path = tempfile.TemporaryDirectory()
        os.symlink(os.path.join(self.test_data_dir, 'input.ensr_transcript.vcf'), os.path.join(temp_path.name, 'input.vcf'))
        command = [
            os.path.join(temp_path.name, 'input.vcf'),
            os.path.join(self.test_data_dir, 'genes.fpkm_tracking'),
            'cufflinks',
            'gene',
        ]
        vcf_expression_annotator.main(command)
        self.assertTrue(cmp(os.path.join(self.test_data_dir, 'input.ensr_transcript.gx.vcf'), os.path.join(temp_path.name, 'input.gx.vcf')))
        temp_path.cleanup()

    def test_warning_mutation_without_matching_expression_value(self):
        logging.disable(logging.NOTSET)
        with LogCapture() as l:
            temp_path = tempfile.TemporaryDirectory()
            os.symlink(os.path.join(self.test_data_dir, 'no_matching_expression.vcf'), os.path.join(temp_path.name, 'input.vcf'))
            command = [
                os.path.join(temp_path.name, 'input.vcf'),
                os.path.join(self.test_data_dir, 'genes.fpkm_tracking'),
                'cufflinks',
                'gene',
            ]
            vcf_expression_annotator.main(command)
            temp_path.cleanup()
            l.check_present(('root', 'WARNING', "1 of 1 transcripts did not have an expression entry for their gene id."))

    def test_multi_sample_vcf(self):
        temp_path = tempfile.TemporaryDirectory()
        os.symlink(os.path.join(self.test_data_dir, 'multiple_samples.vcf'), os.path.join(temp_path.name, 'input.vcf'))
        command = [
            os.path.join(temp_path.name, 'input.vcf'),
            os.path.join(self.test_data_dir, 'genes.fpkm_tracking'),
            'cufflinks',
            'gene',
            '-s', 'H_NJ-HCC1395-HCC1395',
        ]
        vcf_expression_annotator.main(command)
        self.assertTrue(cmp(os.path.join(self.test_data_dir, 'multiple_samples.gx.vcf'), os.path.join(temp_path.name, 'input.gx.vcf')))
        temp_path.cleanup()

    def test_multiple_transcripts(self):
        temp_path = tempfile.TemporaryDirectory()
        os.symlink(os.path.join(self.test_data_dir, 'multiple_transcripts.vcf'), os.path.join(temp_path.name, 'input.vcf'))
        command = [
            os.path.join(temp_path.name, 'input.vcf'),
            os.path.join(self.test_data_dir, 'genes.fpkm_tracking'),
            'cufflinks',
            'gene',
        ]
        vcf_expression_annotator.main(command)
        self.assertTrue(cmp(os.path.join(self.test_data_dir, 'multiple_transcripts.gx.vcf'), os.path.join(temp_path.name, 'input.gx.vcf')))
        temp_path.cleanup()

    def test_multiple_transcripts_with_one_missing_gene(self):
        temp_path = tempfile.TemporaryDirectory()
        os.symlink(os.path.join(self.test_data_dir, 'multiple_transcripts_one_missing_gene.vcf'), os.path.join(temp_path.name, 'input.vcf'))
        command = [
            os.path.join(temp_path.name, 'input.vcf'),
            os.path.join(self.test_data_dir, 'genes.fpkm_tracking'),
            'cufflinks',
            'gene',
        ]
        vcf_expression_annotator.main(command)
        self.assertTrue(cmp(os.path.join(self.test_data_dir, 'multiple_transcripts_one_missing_gene.gx.vcf'), os.path.join(temp_path.name, 'input.gx.vcf')))
        temp_path.cleanup()

    def test_cufflinks_format_gene_mode(self):
        temp_path = tempfile.TemporaryDirectory()
        os.symlink(os.path.join(self.test_data_dir, 'input.vcf'), os.path.join(temp_path.name, 'input.vcf'))
        command = [
            os.path.join(temp_path.name, 'input.vcf'),
            os.path.join(self.test_data_dir, 'genes.fpkm_tracking'),
            'cufflinks',
            'gene',
        ]
        vcf_expression_annotator.main(command)
        self.assertTrue(cmp(os.path.join(self.test_data_dir, 'input.cufflinks.gx.vcf'), os.path.join(temp_path.name, 'input.gx.vcf')))
        temp_path.cleanup()

    def test_cufflinks_format_transcript_mode(self):
        temp_path = tempfile.TemporaryDirectory()
        os.symlink(os.path.join(self.test_data_dir, 'input.vcf'), os.path.join(temp_path.name, 'input.vcf'))
        command = [
            os.path.join(temp_path.name, 'input.vcf'),
            os.path.join(self.test_data_dir, 'isoforms.fpkm_tracking'),
            'cufflinks',
            'transcript',
        ]
        vcf_expression_annotator.main(command)
        self.assertTrue(cmp(os.path.join(self.test_data_dir, 'input.cufflinks.tx.vcf'), os.path.join(temp_path.name, 'input.tx.vcf')))
        temp_path.cleanup()

    def test_kallisto_format_gene_mode(self):
        temp_path = tempfile.TemporaryDirectory()
        os.symlink(os.path.join(self.test_data_dir, 'input.vcf'), os.path.join(temp_path.name, 'input.vcf'))
        command = [
            os.path.join(temp_path.name, 'input.vcf'),
            os.path.join(self.test_data_dir, 'kallisto.genes'),
            'kallisto',
            'gene',
        ]
        vcf_expression_annotator.main(command)
        self.assertTrue(cmp(os.path.join(self.test_data_dir, 'input.kallisto.gx.vcf'), os.path.join(temp_path.name, 'input.gx.vcf')))
        temp_path.cleanup()

    def test_kallisto_format_transcript_mode(self):
        temp_path = tempfile.TemporaryDirectory()
        os.symlink(os.path.join(self.test_data_dir, 'input.vcf'), os.path.join(temp_path.name, 'input.vcf'))
        command = [
            os.path.join(temp_path.name, 'input.vcf'),
            os.path.join(self.test_data_dir, 'kallisto.transcripts'),
            'kallisto',
            'transcript',
        ]
        vcf_expression_annotator.main(command)
        self.assertTrue(cmp(os.path.join(self.test_data_dir, 'input.kallisto.tx.vcf'), os.path.join(temp_path.name, 'input.tx.vcf')))
        temp_path.cleanup()

    def test_stringtie_format_gene_mode(self):
        temp_path = tempfile.TemporaryDirectory()
        os.symlink(os.path.join(self.test_data_dir, 'input.vcf'), os.path.join(temp_path.name, 'input.vcf'))
        command = [
            os.path.join(temp_path.name, 'input.vcf'),
            os.path.join(self.test_data_dir, 'genes.tsv'),
            'stringtie',
            'gene',
        ]
        vcf_expression_annotator.main(command)
        self.assertTrue(cmp(os.path.join(self.test_data_dir, 'input.stringtie.gx.vcf'), os.path.join(temp_path.name, 'input.gx.vcf')))
        temp_path.cleanup()

    def test_stringtie_format_transcript_mode(self):
        temp_path = tempfile.TemporaryDirectory()
        os.symlink(os.path.join(self.test_data_dir, 'input.vcf'), os.path.join(temp_path.name, 'input.vcf'))
        command = [
            os.path.join(temp_path.name, 'input.vcf'),
            os.path.join(self.test_data_dir, 'transcripts.gtf'),
            'stringtie',
            'transcript',
        ]
        vcf_expression_annotator.main(command)
        self.assertTrue(cmp(os.path.join(self.test_data_dir, 'input.stringtie.tx.vcf'), os.path.join(temp_path.name, 'input.tx.vcf')))
        temp_path.cleanup()

    def test_kallisto_with_transcript_version_in_vcf(self):
        temp_path = tempfile.TemporaryDirectory()
        os.symlink(os.path.join(self.test_data_dir, 'input.transcript_version.vcf'), os.path.join(temp_path.name, 'input.vcf'))
        command = [
            os.path.join(temp_path.name, 'input.vcf'),
            os.path.join(self.test_data_dir, 'kallisto.transcripts'),
            'kallisto',
            'transcript',
            "--ignore-ensembl-id-version",
        ]
        vcf_expression_annotator.main(command)
        self.assertTrue(cmp(os.path.join(self.test_data_dir, 'input.kallisto.with_version.tx.vcf'), os.path.join(temp_path.name, 'input.tx.vcf')))
        temp_path.cleanup()

    def test_kallisto_with_transcript_version_in_expression_file(self):
        temp_path = tempfile.TemporaryDirectory()
        os.symlink(os.path.join(self.test_data_dir, 'input.vcf'), os.path.join(temp_path.name, 'input.vcf'))
        command = [
            os.path.join(temp_path.name, 'input.vcf'),
            os.path.join(self.test_data_dir, 'kallisto.transcript_version.transcripts'),
            'kallisto',
            'transcript',
            "--ignore-ensembl-id-version",
        ]
        vcf_expression_annotator.main(command)
        self.assertTrue(cmp(os.path.join(self.test_data_dir, 'input.kallisto.tx.vcf'), os.path.join(temp_path.name, 'input.tx.vcf')))
        temp_path.cleanup()

    def test_warning_kallisto_with_transcript_version_in_vcf(self):
        logging.disable(logging.NOTSET)
        with LogCapture() as l:
            temp_path = tempfile.TemporaryDirectory()
            os.symlink(os.path.join(self.test_data_dir, 'input.transcript_version.vcf'), os.path.join(temp_path.name, 'input.vcf'))
            command = [
                os.path.join(temp_path.name, 'input.vcf'),
                os.path.join(self.test_data_dir, 'kallisto.transcripts'),
                'kallisto',
                'transcript',
            ]
            vcf_expression_annotator.main(command)
            temp_path.cleanup()
            l.check_present(('root', 'WARNING', "1 of 1 transcripts did not have an expression entry for their transcript id."))

    def test_warning_kallisto_with_transcript_version_in_expression_file(self):
        logging.disable(logging.NOTSET)
        with LogCapture() as l:
            temp_path = tempfile.TemporaryDirectory()
            os.symlink(os.path.join(self.test_data_dir, 'input.vcf'), os.path.join(temp_path.name, 'input.vcf'))
            command = [
                os.path.join(temp_path.name, 'input.vcf'),
                os.path.join(self.test_data_dir, 'kallisto.transcript_version.transcripts'),
                'kallisto',
                'transcript',
            ]
            vcf_expression_annotator.main(command)
            temp_path.cleanup()
            l.check_present(('root', 'WARNING', "1 of 1 transcripts did not have an expression entry for their transcript id."))

    def test_skip_ENSR_transcript(self):
        temp_path = tempfile.TemporaryDirectory()
        os.symlink(os.path.join(self.test_data_dir, 'input.ensr_transcript.vcf'), os.path.join(temp_path.name, 'input.vcf'))
        command = [
            os.path.join(temp_path.name, 'input.vcf'),
            os.path.join(self.test_data_dir, 'kallisto.transcripts'),
            'kallisto',
            'transcript',
        ]
        vcf_expression_annotator.main(command)
        self.assertTrue(cmp(os.path.join(self.test_data_dir, 'input.ensr_transcript.tx.vcf'), os.path.join(temp_path.name, 'input.tx.vcf')))
        temp_path.cleanup()
