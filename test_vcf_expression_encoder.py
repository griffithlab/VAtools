import unittest
import sys
import os
import py_compile
import vcf_expression_encoder

class VcfExpressionEncoderTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        base_dir          = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
        cls.executable    = os.path.join(base_dir, 'vcf_expression_encoder.py')
        cls.test_data_dir = os.path.join(base_dir, 'test_data')

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
            vcf_expression_encoder.main(command)
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
            vcf_expression_encoder.main(command)
        self.assertTrue('--expression-column is not set. This is required when using the `custom` format.' == str(context.exception))

    def test_error_more_than_one_sample_without_sample_name(self):
        with self.assertRaises(Exception) as context:
            command = [
                os.path.join(self.test_data_dir, 'multiple_samples.vcf'),
                os.path.join(self.test_data_dir, 'genes.fpkm_tracking'),
                'cufflinks',
                'gene',
            ]
            vcf_expression_encoder.main(command)
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
            vcf_expression_encoder.main(command)
        self.assertTrue('does not contain a sample column for sample nonexistent_sample.' in str(context.exception))

    def test_error_no_csq(self):
        with self.assertRaises(Exception) as context:
            command = [
                os.path.join(self.test_data_dir, 'no_csq.vcf'),
                os.path.join(self.test_data_dir, 'genes.fpkm_tracking'),
                'cufflinks',
                'gene',
            ]
            vcf_expression_encoder.main(command)
        self.assertTrue('is not VEP-annotated. Please annotate the VCF with VEP before running this tool.' in str(context.exception))

    def test_error_already_TX_annotated(self):
        with self.assertRaises(Exception) as context:
            command = [
                os.path.join(self.test_data_dir, 'input.tx.vcf'),
                os.path.join(self.test_data_dir, 'isoforms.fpkm_tracking'),
                'cufflinks',
                'transcript',
            ]
            vcf_expression_encoder.main(command)
        self.assertTrue('is already transcript expression annotated. TX format header already exists.' in str(context.exception))

    def test_error_already_GX_annotated(self):
        with self.assertRaises(Exception) as context:
            command = [
                os.path.join(self.test_data_dir, 'input.gx.vcf'),
                os.path.join(self.test_data_dir, 'genes.fpkm_tracking'),
                'cufflinks',
                'gene',
            ]
            vcf_expression_encoder.main(command)
        self.assertTrue('is already gene expression annotated. GX format header already exists.' in str(context.exception))
