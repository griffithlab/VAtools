import unittest
import sys
import os
import py_compile
from vcf_annotation_tools import vcf_info_annotator
import tempfile
from filecmp import cmp

class VcfInfoEncoderTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        base_dir          = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
        cls.executable    = os.path.join(base_dir, 'vcf_annotation_tools', 'vcf_info_annotator.py')
        cls.test_data_dir = os.path.join(base_dir, 'tests', 'test_data')

    def test_source_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_error_already_INFO_annotated(self):
        with self.assertRaises(Exception) as context:
            command = [
                os.path.join(self.test_data_dir, 'input.vcf'),
                os.path.join(self.test_data_dir, 'info.tsv'),
                'QSI', 
                "test",
                'Float',
                'ztest.vcf'
            ]
            vcf_info_annotator.main(command)
        self.assertTrue('INFO already contains a QSI field. Choose a different label' in str(context.exception))

