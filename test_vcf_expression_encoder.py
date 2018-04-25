import unittest
import sys
import os
import py_compile

class VcfExpressionEncoderTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        base_dir       = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
        cls.executable = os.path.join(base_dir, 'vcf_expression_encoder.py')

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
