import unittest
import sys
import os
import py_compile
from vatools import vcf_info_annotator
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
                'CSQ', 
                '-d', "test",
                '-f', 'Integer',
                '-o', 'ztest.vcf'
            ]
            vcf_info_annotator.main(command)
        self.assertTrue('INFO already contains a CSQ field. Choose a different label, or use the --overwrite flag to retain this field and overwrite values' in str(context.exception))

    def test_error_new_field_no_description(self):
        with self.assertRaises(Exception) as context:
            command = [
                os.path.join(self.test_data_dir, 'input.vcf'),
                os.path.join(self.test_data_dir, 'info.tsv'),
                'TEST', 
                '-o', 'ztest.vcf'
            ]
            vcf_info_annotator.main(command)
        self.assertTrue("the --description and --value_format arguments are required unless updating/overwriting an existing field (with flag --overwrite)" in str(context.exception))


    def test_overwrite_when_field_doesnt_exist(self):
        with self.assertRaises(Exception) as context:
            command = [
                os.path.join(self.test_data_dir, 'input.vcf'),
                os.path.join(self.test_data_dir, 'info.tsv'),
                'TEST', 
                '-o', 'ztest.vcf',
                '-w'
            ]
            vcf_info_annotator.main(command)
        self.assertTrue("INFO field TEST does not exist and thus cannot be overwritten!" in str(context.exception))

    def test_simple_caseq(self):
        temp_path = tempfile.TemporaryDirectory()
        print(temp_path)
        command = [
            os.path.join(self.test_data_dir, 'input.vcf'),
            os.path.join(self.test_data_dir, 'info.tsv'),
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
        print(temp_path)
        command = [
            os.path.join(self.test_data_dir, 'input.vcf'),
            os.path.join(self.test_data_dir, 'info3.tsv'),
            'TEST', 
            '-d', "test",
            '-f', 'String',
            '-o', os.path.join(temp_path.name, 'info3_output.vcf')
        ]
        vcf_info_annotator.main(command)
        self.assertTrue(cmp(os.path.join(self.test_data_dir, 'info3_output.vcf'), os.path.join(temp_path.name, 'info3_output.vcf')))
        temp_path.cleanup()

    def test_addwhile_overwriteset(self):
        temp_path = tempfile.TemporaryDirectory()
        print(temp_path)
        command = [
            os.path.join(self.test_data_dir, 'info2_input.vcf'),
            os.path.join(self.test_data_dir, 'info2.tsv'),
            'MQ0', 
            '-w',
            '-o', os.path.join(temp_path.name, 'info2_output.vcf')
        ]
        vcf_info_annotator.main(command)
        self.assertTrue(cmp(os.path.join(self.test_data_dir, 'info2_output.vcf'), os.path.join(temp_path.name, 'info2_output.vcf')))
        temp_path.cleanup()

    def test_overwrite_existing_field(self):
        temp_path = tempfile.TemporaryDirectory()
        print(temp_path)
        command = [
            os.path.join(self.test_data_dir, 'input.vcf'),
            os.path.join(self.test_data_dir, 'info.tsv'),
            'CSQ', 
            '-d', "test",
            '-f', 'Integer',
            '-w',
            '-o', os.path.join(temp_path.name, 'info_annotation.vcf')
        ]
        vcf_info_annotator.main(command)
        self.assertTrue(cmp(os.path.join(self.test_data_dir, 'info_overwrite.vcf'), os.path.join(temp_path.name, 'info_annotation.vcf')))
        temp_path.cleanup()
