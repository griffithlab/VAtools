import unittest
import sys
import os
import py_compile
from vatools import vcf_genotype_annotator
import tempfile
from filecmp import cmp

class VcfExpressionEncoderTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        base_dir          = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
        cls.executable    = os.path.join(base_dir, 'vatools', 'vcf_genotype_annotator.py')
        cls.test_data_dir = os.path.join(base_dir, 'tests', 'test_data')

    def test_source_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_error_sample_name_already_exists_with_GT_field(self):
        with self.assertRaises(Exception) as context:
            command = [
                os.path.join(self.test_data_dir, 'input.vcf'),
                'H_NJ-HCC1395-HCC1395',
                '0/1',
            ]
            vcf_genotype_annotator.main(command)
        self.assertTrue('VCF already contains a sample column for sample H_NJ-HCC1395-HCC1395 with a GT field.' in str(context.exception))

    def test_no_sample_vcf(self):
        temp_path = tempfile.TemporaryDirectory()
        os.symlink(os.path.join(self.test_data_dir, 'input.no_sample.vcf'), os.path.join(temp_path.name, 'input.vcf'))
        command = [
            os.path.join(temp_path.name, 'input.vcf'),
            'TUMOR',
            '0/1',
        ]
        vcf_genotype_annotator.main(command)
        self.assertTrue(cmp(os.path.join(self.test_data_dir, 'no_sample.genotype.vcf'), os.path.join(temp_path.name, 'input.genotype.vcf')))
        temp_path.cleanup()

    def test_single_sample_vcf(self):
        temp_path = tempfile.TemporaryDirectory()
        os.symlink(os.path.join(self.test_data_dir, 'input.vcf'), os.path.join(temp_path.name, 'input.vcf'))
        command = [
            os.path.join(temp_path.name, 'input.vcf'),
            'TUMOR',
            '0/1',
        ]
        vcf_genotype_annotator.main(command)
        self.assertTrue(cmp(os.path.join(self.test_data_dir, 'single_sample.genotype.vcf'), os.path.join(temp_path.name, 'input.genotype.vcf')))
        temp_path.cleanup()

    def test_no_gt_in_format(self):
        temp_path = tempfile.TemporaryDirectory()
        os.symlink(os.path.join(self.test_data_dir, 'input.no_gt_in_format.vcf'), os.path.join(temp_path.name, 'input.vcf'))
        command = [
            os.path.join(temp_path.name, 'input.vcf'),
            'TUMOR',
            '0/1',
        ]
        vcf_genotype_annotator.main(command)
        self.assertTrue(cmp(os.path.join(self.test_data_dir, 'no_gt_in_format.genotype.vcf'), os.path.join(temp_path.name, 'input.genotype.vcf')))
        temp_path.cleanup()

    def test_adding_gt_in_existing_sample(self):
        temp_path = tempfile.TemporaryDirectory()
        os.symlink(os.path.join(self.test_data_dir, 'input.no_gt_in_format.vcf'), os.path.join(temp_path.name, 'input.vcf'))
        command = [
            os.path.join(temp_path.name, 'input.vcf'),
            'H_NJ-HCC1395-HCC1395',
            '0/1',
        ]
        vcf_genotype_annotator.main(command)
        self.assertTrue(cmp(os.path.join(self.test_data_dir, 'existing_sample.genotype.vcf'), os.path.join(temp_path.name, 'input.genotype.vcf')))
        temp_path.cleanup()
