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
