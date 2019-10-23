import unittest
import sys
import os
import py_compile
from vatools import vep_annotation_reporter
import tempfile
from filecmp import cmp
import io
import logging
from testfixtures import LogCapture, StringComparison as S

class VcfExpressionEncoderTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        base_dir          = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
        cls.executable    = os.path.join(base_dir, 'vatools', 'vep_annotation_reporter.py')
        cls.test_data_dir = os.path.join(base_dir, 'tests', 'test_data', 'vep_annotation_reporter')

    def test_source_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_error_no_vep_annotation(self):
        with self.assertRaises(Exception) as context:
            command = [
                os.path.join(self.test_data_dir, '..', 'no_csq.vcf'),
                'Consequence',
                '-t', os.path.join(self.test_data_dir, 'variants.tsv'),
            ]
            vep_annotation_reporter.main(command)
        self.assertTrue('is not VEP-annotated. Please annotate the VCF with VEP before running this tool.' in str(context.exception))

    def test_error_missing_column_in_tsv(self):
        with self.assertRaises(Exception) as context:
            command = [
                os.path.join(self.test_data_dir, 'input.vcf.gz'),
                'Consequence',
                '-t', os.path.join(self.test_data_dir, 'variants.no_ALT.tsv'),
            ]
            vep_annotation_reporter.main(command)
        self.assertTrue("doesn't contain required column 'ALT'." in str(context.exception))

    def test_single_vep_field(self):
        temp_path = tempfile.TemporaryDirectory()
        os.symlink(os.path.join(self.test_data_dir, 'input.vcf.gz'), os.path.join(temp_path.name, 'input.vcf.gz'))
        command = [
            os.path.join(temp_path.name, 'input.vcf.gz'),
            'Consequence',
            '-t', os.path.join(self.test_data_dir, 'variants.tsv'),
        ]
        vep_annotation_reporter.main(command)
        self.assertTrue(cmp(os.path.join(self.test_data_dir, 'output.single_field.tsv'), os.path.join(temp_path.name, 'input.tsv')))
        temp_path.cleanup()

    def test_multiple_vep_fields(self):
        temp_path = tempfile.TemporaryDirectory()
        os.symlink(os.path.join(self.test_data_dir, 'input.vcf.gz'), os.path.join(temp_path.name, 'input.vcf.gz'))
        command = [
            os.path.join(temp_path.name, 'input.vcf.gz'),
            'Consequence',
            'Gene',
            '-t', os.path.join(self.test_data_dir, 'variants.tsv'),
        ]
        vep_annotation_reporter.main(command)
        self.assertTrue(cmp(os.path.join(self.test_data_dir, 'output.multiple_fields.tsv'), os.path.join(temp_path.name, 'input.tsv')))
        temp_path.cleanup()

    def test_nonexistent_vep_field(self):
        temp_path = tempfile.TemporaryDirectory()
        os.symlink(os.path.join(self.test_data_dir, 'input.vcf.gz'), os.path.join(temp_path.name, 'input.vcf.gz'))
        command = [
            os.path.join(temp_path.name, 'input.vcf.gz'),
            'Nonexistent_Field',
            '-t', os.path.join(self.test_data_dir, 'variants.tsv'),
        ]
        vep_annotation_reporter.main(command)
        self.assertTrue(cmp(os.path.join(self.test_data_dir, 'output.nonexistent_field.tsv'), os.path.join(temp_path.name, 'input.tsv')))
        temp_path.cleanup()

    def test_multiple_multiallelic_site(self):
        temp_path = tempfile.TemporaryDirectory()
        os.symlink(os.path.join(self.test_data_dir, '..', 'input.multiallelic.vcf.gz'), os.path.join(temp_path.name, 'input.vcf.gz'))
        command = [
            os.path.join(temp_path.name, 'input.vcf.gz'),
            'Consequence',
            '-t', os.path.join(self.test_data_dir, 'variants.multiallelic.tsv'),
        ]
        vep_annotation_reporter.main(command)
        self.assertTrue(cmp(os.path.join(self.test_data_dir, 'output.multiallelic.tsv'), os.path.join(temp_path.name, 'input.tsv')))
        temp_path.cleanup()

    def test_no_input_tsv(self):
        temp_path = tempfile.TemporaryDirectory()
        os.symlink(os.path.join(self.test_data_dir, 'input.vcf.gz'), os.path.join(temp_path.name, 'input.vcf.gz'))
        command = [
            os.path.join(temp_path.name, 'input.vcf.gz'),
            'Consequence',
        ]
        vep_annotation_reporter.main(command)
        self.assertTrue(cmp(os.path.join(self.test_data_dir, 'output.no_input_tsv.tsv'), os.path.join(temp_path.name, 'input.tsv')))
        temp_path.cleanup()

    def test_vcf_entries_with_no_csq_runs_successfully(self):
        temp_path = tempfile.TemporaryDirectory()
        os.symlink(os.path.join(self.test_data_dir, 'input.variants_without_csq.vcf.gz'), os.path.join(temp_path.name, 'input.vcf.gz'))
        command = [
            os.path.join(temp_path.name, 'input.vcf.gz'),
            'Consequence',
        ]
        vep_annotation_reporter.main(command)
        self.assertTrue(cmp(os.path.join(self.test_data_dir, 'output.variants_without_csq.tsv'), os.path.join(temp_path.name, 'input.tsv')))
        temp_path.cleanup()

    def test_vcf_with_multiple_transcripts_and_no_pick(self):
        temp_path = tempfile.TemporaryDirectory()
        os.symlink(os.path.join(self.test_data_dir, 'input.merge_multiple_transcripts.vcf.gz'), os.path.join(temp_path.name, 'input.vcf.gz'))
        command = [
            os.path.join(temp_path.name, 'input.vcf.gz'),
            'SYMBOL',
        ]
        vep_annotation_reporter.main(command)
        self.assertTrue(cmp(os.path.join(self.test_data_dir, 'output.merge_multiple_transcripts.tsv'), os.path.join(temp_path.name, 'input.tsv')))
        temp_path.cleanup()
