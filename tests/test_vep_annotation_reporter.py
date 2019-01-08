import unittest
import sys
import os
import py_compile
from vcf_annotation_tools import vep_annotation_reporter
import tempfile
from filecmp import cmp
import io
import logging
from testfixtures import LogCapture, StringComparison as S

class VcfExpressionEncoderTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        base_dir          = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
        cls.executable    = os.path.join(base_dir, 'vcf_annotation_tools', 'vep_annotation_reporter.py')
        cls.test_data_dir = os.path.join(base_dir, 'tests', 'test_data', 'vep_annotation_reporter')

    def test_source_compiles(self):
        self.assertTrue(py_compile.compile(self.executable))
