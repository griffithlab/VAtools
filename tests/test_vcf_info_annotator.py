import unittest
import sys
import os
import py_compile
from vatools import vcf_info_annotator
from vatools.vcf_info_annotator import coerce_value
import tempfile
from filecmp import cmp

class VcfInfoEncoderTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        base_dir          = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
        cls.executable    = os.path.join(base_dir, 'vatools', 'vcf_info_annotator.py')
        cls.test_data_dir = os.path.join(base_dir, 'tests', 'test_data')

    def test_source_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_error_already_INFO_annotated(self):
        with self.assertRaises(Exception) as context:
            command = [
                os.path.join(self.test_data_dir, 'input.vcf'),
                os.path.join(self.test_data_dir, 'info.tsv'),
                '-m', 'value:CSQ:String:CSQ annotation',
                '-o', 'ztest.vcf'
            ]
            vcf_info_annotator.main(command)
        self.assertTrue('INFO already contains a CSQ field. Choose a different label, or use the --overwrite flag to retain this field and overwrite values' in str(context.exception))

    def test_error_invalid_column_mapping_format(self):
        with self.assertRaises(Exception) as context:
            command = [
                os.path.join(self.test_data_dir, 'input.vcf'),
                os.path.join(self.test_data_dir, 'info.tsv'),
                '-m', 'value:TEST',
                '-o', 'ztest.vcf'
            ]
            vcf_info_annotator.main(command)
        self.assertTrue("Invalid column mapping 'value:TEST'. Expected format: source_col:info_field:type:description[:source[:version]]" in str(context.exception))

    def test_simple_caseq(self):
        temp_path = tempfile.TemporaryDirectory()
        command = [
            os.path.join(self.test_data_dir, 'input.vcf'),
            os.path.join(self.test_data_dir, 'info.tsv'),
            '-m', 'value:TEST:Integer:test',
            '-o', os.path.join(temp_path.name, 'info_annotation.vcf')
        ]
        vcf_info_annotator.main(command)
        self.assertTrue(cmp(os.path.join(self.test_data_dir, 'info_annotation.vcf'), os.path.join(temp_path.name, 'info_annotation.vcf')))
        temp_path.cleanup()

    def test_gzipped_values_file(self):
        temp_path = tempfile.TemporaryDirectory()
        command = [
            os.path.join(self.test_data_dir, 'input.vcf'),
            os.path.join(self.test_data_dir, 'info.tsv.gz'),
            'TEST',
            '-d', "test",
            '-f', 'Integer',
            '-o', os.path.join(temp_path.name, 'info_annotation.vcf')
        ]
        vcf_info_annotator.main(command)
        self.assertTrue(cmp(os.path.join(self.test_data_dir, 'info_annotation.vcf'), os.path.join(temp_path.name, 'info_annotation.vcf')))
        temp_path.cleanup()

    def test_simple_string(self):
        temp_path = tempfile.TemporaryDirectory()
        command = [
            os.path.join(self.test_data_dir, 'input.vcf'),
            os.path.join(self.test_data_dir, 'info3.tsv'),
            '-m', 'value:TEST:String:test',
            '-o', os.path.join(temp_path.name, 'info3_output.vcf')
        ]
        vcf_info_annotator.main(command)
        self.assertTrue(cmp(os.path.join(self.test_data_dir, 'info3_output.vcf'), os.path.join(temp_path.name, 'info3_output.vcf')))
        temp_path.cleanup()

    def test_addwhile_overwriteset(self):
        temp_path = tempfile.TemporaryDirectory()
        command = [
            os.path.join(self.test_data_dir, 'info2_input.vcf'),
            os.path.join(self.test_data_dir, 'info2.tsv'),
            '-m', 'value:MQ0:Integer:Mapping quality',
            '-w',
            '-o', os.path.join(temp_path.name, 'info2_output.vcf')
        ]
        vcf_info_annotator.main(command)
        self.assertTrue(cmp(os.path.join(self.test_data_dir, 'info2_output.vcf'), os.path.join(temp_path.name, 'info2_output.vcf')))
        temp_path.cleanup()

    def test_overwrite_existing_field(self):
        temp_path = tempfile.TemporaryDirectory()
        command = [
            os.path.join(self.test_data_dir, 'input.vcf'),
            os.path.join(self.test_data_dir, 'info.tsv'),
            '-m', 'value:CSQ:String:CSQ annotation',
            '-w',
            '-o', os.path.join(temp_path.name, 'info_annotation.vcf')
        ]
        vcf_info_annotator.main(command)
        self.assertTrue(cmp(os.path.join(self.test_data_dir, 'info_overwrite.vcf'), os.path.join(temp_path.name, 'info_annotation.vcf')))
        temp_path.cleanup()

    def test_source_and_version_in_mapping(self):
        temp_path = tempfile.TemporaryDirectory()
        command = [
            os.path.join(self.test_data_dir, 'input.vcf'),
            os.path.join(self.test_data_dir, 'multi_col.tsv'),
            '-m', 'mapping_quality:MQ0:Integer:Mapping quality:TestSource:1.0',
            '-o', os.path.join(temp_path.name, 'source_version.info.vcf'),
        ]
        vcf_info_annotator.main(command)
        self.assertTrue(cmp(os.path.join(self.test_data_dir, 'source_version.info.vcf'), os.path.join(temp_path.name, 'source_version.info.vcf')))
        temp_path.cleanup()

    def test_coerce_value_special_values(self):
        import math
        # VCF universal missing marker → "." for all types
        for sentinel in ('.', ''):
            self.assertEqual(coerce_value(sentinel, 'Integer'), '.')
            self.assertEqual(coerce_value(sentinel, 'Float'), '.')
            self.assertEqual(coerce_value(sentinel, 'String'), '.')
        # Float special values
        for nan_str in ('NA', 'na', 'NaN', 'nan'):
            self.assertTrue(math.isnan(coerce_value(nan_str, 'Float')))
        for inf_str in ('Inf', 'inf', '+Inf', '+inf'):
            self.assertEqual(coerce_value(inf_str, 'Float'), float('inf'))
        for neginf_str in ('-Inf', '-inf'):
            self.assertEqual(coerce_value(neginf_str, 'Float'), float('-inf'))
        # Integer: NA/NaN → None; Inf → error
        for nan_str in ('NA', 'na', 'NaN', 'nan'):
            self.assertIsNone(coerce_value(nan_str, 'Integer'))
        for inf_str in ('Inf', '+Inf', '-Inf'):
            with self.assertRaises(ValueError):
                coerce_value(inf_str, 'Integer')
        # normal values coerce correctly
        self.assertEqual(coerce_value('42', 'Integer'), 42)
        self.assertAlmostEqual(coerce_value('3.14', 'Float'), 3.14)
        self.assertEqual(coerce_value('hello', 'String'), 'hello')

    def test_clear_existing_error_without_overwrite(self):
        with self.assertRaises(SystemExit):
            vcf_info_annotator.main([
                os.path.join(self.test_data_dir, 'info2_input.vcf'),
                os.path.join(self.test_data_dir, 'no_match.tsv'),
                '-m', 'value:MQ0:Integer:Mapping quality',
                '--clear-existing',
            ])

    def test_clear_existing_removes_unmatched_values(self):
        temp_path = tempfile.TemporaryDirectory()
        command = [
            os.path.join(self.test_data_dir, 'info2_input.vcf'),
            os.path.join(self.test_data_dir, 'no_match.tsv'),
            '-m', 'value:MQ0:Integer:Mapping quality',
            '-w', '--clear-existing',
            '-o', os.path.join(temp_path.name, 'info2_cleared.vcf'),
        ]
        vcf_info_annotator.main(command)
        self.assertTrue(cmp(os.path.join(self.test_data_dir, 'info2_cleared.vcf'), os.path.join(temp_path.name, 'info2_cleared.vcf')))
        temp_path.cleanup()

    def test_multi_column_mappings(self):
        temp_path = tempfile.TemporaryDirectory()
        command = [
            os.path.join(self.test_data_dir, 'input.vcf'),
            os.path.join(self.test_data_dir, 'multi_col.tsv'),
            '-m', 'mapping_quality:MQ0:Integer:Mapping quality,clinvar_classification:CVCLASS:String:ClinVar variant classification',
            '-o', os.path.join(temp_path.name, 'multi_col.info.vcf'),
        ]
        vcf_info_annotator.main(command)
        self.assertTrue(cmp(os.path.join(self.test_data_dir, 'multi_col.info.vcf'), os.path.join(temp_path.name, 'multi_col.info.vcf')))
        temp_path.cleanup()
