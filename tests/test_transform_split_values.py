import unittest
import sys
import os
import py_compile
from vatools import transform_split_values
import tempfile
from filecmp import cmp
import io
import logging
from testfixtures import LogCapture, StringComparison as S

class TransformSplitValuesTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        base_dir          = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
        cls.executable    = os.path.join(base_dir, 'vatools', 'transform_split_values.py')
        cls.test_data_dir = os.path.join(base_dir, 'tests', 'test_data', 'transform_split_values')

    def test_source_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_error_more_than_one_sample_without_sample_name(self):
        with self.assertRaises(Exception) as context:
            command = [
                os.path.join(self.test_data_dir, '..', 'multiple_samples.vcf'),
                'AD',
                'ref',
            ]
            transform_split_values.main(command)
        self.assertTrue('contains more than one sample. Please use the -s option to specify which sample to annotate.' in str(context.exception))

    def test_error_more_than_one_sample_with_wrong_sample_name(self):
        with self.assertRaises(Exception) as context:
            command = [
                os.path.join(self.test_data_dir, '..', 'multiple_samples.vcf'),
                'AD',
                'ref',
                '-s', 'nonexistent_sample',
            ]
            transform_split_values.main(command)
        self.assertTrue('does not contain a sample column for sample nonexistent_sample.' in str(context.exception))

    def test_error_nonexistent_format_field(self):
        with self.assertRaises(Exception) as context:
            command = [
                os.path.join(self.test_data_dir, '..', 'input.vcf'),
                'nonexistent_field',
                'ref',
            ]
            transform_split_values.main(command)
        self.assertTrue('does not contain a format field nonexistent_field.' in str(context.exception))

    def test_error_field_not_of_number_R(self):
        with self.assertRaises(Exception) as context:
            command = [
                os.path.join(self.test_data_dir, '..', 'input.readcount.vcf'),
                'AD',
                'ref',
            ]
            transform_split_values.main(command)
        self.assertTrue('format field AD incompatible with operation ref. Not of Number R.' in str(context.exception))

    def test_error_field_not_a_number(self):
        with self.assertRaises(Exception) as context:
            command = [
                os.path.join(self.test_data_dir, '..', 'input.readcount.vcf'),
                'GT',
                'sum',
            ]
            transform_split_values.main(command)
        self.assertTrue('format field GT incompatible with operation sum. Not of Type Integer or Float.' in str(context.exception))

    def test_error_missing_column_in_tsv(self):
        with self.assertRaises(Exception) as context:
            command = [
                os.path.join(self.test_data_dir, 'input.empty_ad.vcf'),
                'AD',
                'ref',
                '-t', os.path.join(self.test_data_dir, '..', 'vep_annotation_reporter', 'variants.no_ALT.tsv'),
            ]
            transform_split_values.main(command)
        self.assertTrue("doesn't contain required column 'ALT'." in str(context.exception))

    def test_all_operations(self):
        output_file = tempfile.NamedTemporaryFile()
        command = [
            os.path.join(self.test_data_dir, '..', 'single_sample.dna.readcount.vcf'),
            'AD',
            'ref',
            'alt',
            'sum',
            'min',
            'max',
            'mean',
            'median',
            'stdev',
            'ref_ratio',
            'alt_ratio',
            '-o', output_file.name
        ]
        transform_split_values.main(command)
        self.assertTrue(cmp(os.path.join(self.test_data_dir, 'output.all_operations.tsv'), output_file.name))

    def test_multiallelic_site(self):
        output_file = tempfile.NamedTemporaryFile()
        command = [
            os.path.join(self.test_data_dir, '..', 'input.multiallelic.vcf.gz'),
            'AD',
            'ref',
            '-o', output_file.name
        ]
        transform_split_values.main(command)
        self.assertTrue(cmp(os.path.join(self.test_data_dir, 'output.multiallelic.tsv'), output_file.name))

    def test_variant_with_empty_ad(self):
        output_file = tempfile.NamedTemporaryFile()
        command = [
            os.path.join(self.test_data_dir, 'input.empty_ad.vcf'),
            'AD',
            'ref',
            '-o', output_file.name
        ]
        transform_split_values.main(command)
        self.assertTrue(cmp(os.path.join(self.test_data_dir, 'output.empty_ad.tsv'), output_file.name))

    def test_variant_without_ad(self):
        output_file = tempfile.NamedTemporaryFile()
        command = [
            os.path.join(self.test_data_dir, 'input.no_ad.vcf'),
            'AD',
            'ref',
            '-o', output_file.name
        ]
        transform_split_values.main(command)
        self.assertTrue(cmp(os.path.join(self.test_data_dir, 'output.no_ad.tsv'), output_file.name))

    def test_input_tsv(self):
        output_file = tempfile.NamedTemporaryFile()
        command = [
            os.path.join(self.test_data_dir, '..', 'single_sample.dna.readcount.vcf'),
            'AD',
            'ref',
            '-t', os.path.join(self.test_data_dir, 'variants.tsv'),
            '-o', output_file.name
        ]
        transform_split_values.main(command)
        self.assertTrue(cmp(os.path.join(self.test_data_dir, 'output.with_input_tsv.tsv'), output_file.name))
