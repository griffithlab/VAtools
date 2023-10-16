import unittest
import sys
import os
import py_compile
from vatools import ref_transcript_mismatch_reporter
import tempfile
from filecmp import cmp
import io
import logging
from testfixtures import LogCapture, StringComparison as S
from io import StringIO
from unittest.mock import patch

class RefTranscriptMismatchReporterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        base_dir          = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
        cls.executable    = os.path.join(base_dir, 'vatools', 'ref_transcript_mismatch_reporter.py')
        cls.test_data_dir = os.path.join(base_dir, 'tests', 'test_data', 'ref_transcript_mismatch_reporter')

    def test_source_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))

    def test_error_no_vep_annotation(self):
        with self.assertRaises(Exception) as context:
            command = [
                os.path.join(self.test_data_dir, '..', 'no_csq.vcf'),
            ]
            ref_transcript_mismatch_reporter.main(command)
        self.assertTrue('is not VEP-annotated. Please annotate the VCF with VEP before running this tool.' in str(context.exception))

    #TODO: add test for no Wildtype plugin

    def test_no_filter(self):
        logging.disable(logging.NOTSET)
        with LogCapture() as l:
            temp_path = tempfile.TemporaryDirectory()
            os.symlink(os.path.join(self.test_data_dir, 'csq_mismatch.vcf'), os.path.join(temp_path.name, 'input.vcf'))
            command = [
                os.path.join(temp_path.name, 'input.vcf'),
            ]
            ref_transcript_mismatch_reporter.main(command)
            self.assertTrue(cmp(os.path.join(self.test_data_dir, 'output.mismatches.tsv'), os.path.join(temp_path.name, 'input.mismatches.tsv')))
            temp_path.cleanup()
            self.assertTrue("Total number of variants: 2" in str(l))
            self.assertTrue("Total number of processable variants (at least one missense, inframe indels, or frameshift transcript): 1" in str(l))

    def test_soft_filter(self):
        temp_path = tempfile.TemporaryDirectory()
        os.symlink(os.path.join(self.test_data_dir, 'csq_mismatch.vcf'), os.path.join(temp_path.name, 'input.vcf'))
        command = [
            os.path.join(temp_path.name, 'input.vcf'),
            "-f",
            "soft"
        ]
        ref_transcript_mismatch_reporter.main(command)
        self.assertTrue(cmp(os.path.join(self.test_data_dir, 'output.soft.filtered.vcf'), os.path.join(temp_path.name, 'input.filtered.vcf')))
        temp_path.cleanup()

    def test_hard_filter(self):
        temp_path = tempfile.TemporaryDirectory()
        os.symlink(os.path.join(self.test_data_dir, 'csq_mismatch.vcf'), os.path.join(temp_path.name, 'input.vcf'))
        command = [
            os.path.join(temp_path.name, 'input.vcf'),
            "-f",
            "hard"
        ]
        ref_transcript_mismatch_reporter.main(command)
        self.assertTrue(cmp(os.path.join(self.test_data_dir, 'output.hard.filtered.vcf'), os.path.join(temp_path.name, 'input.filtered.vcf')))
        temp_path.cleanup()

    def test_no_protein_position_filter(self):
        temp_path = tempfile.TemporaryDirectory()
        os.symlink(os.path.join(self.test_data_dir, 'input.no_protein_pos.vcf'), os.path.join(temp_path.name, 'input.vcf'))
        command = [
            os.path.join(temp_path.name, 'input.vcf'),
            "-f",
            "hard"
        ]
        ref_transcript_mismatch_reporter.main(command)
        self.assertTrue(cmp(os.path.join(self.test_data_dir, 'output.no_protein_pos.filtered.vcf'), os.path.join(temp_path.name, 'input.filtered.vcf')))
        temp_path.cleanup()
