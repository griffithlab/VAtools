import contextlib
import csv
import gzip
import io
from pathlib import Path
import tempfile
import unittest

import vcfpy

from vatools import vcf_facets_maf_annotator as annotator


class VcfFacetsMafAnnotatorTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.vcf = self.root / 'input.vcf'
        self.maf = self.root / 'input.maf'
        self.output = self.root / 'output.vcf'
        self.mapping = 'ccf_Mcopies:CCF:Float:FACETS CCF'

    def write_vcf(self, records, extra_headers=''):
        self.vcf.write_text(
            '##fileformat=VCFv4.2\n'
            '##contig=<ID=chr1>\n'
            '##contig=<ID=chrM>\n'
            '##FILTER=<ID=LowQual,Description="Low quality">\n'
            '##INFO=<ID=KEEP,Number=1,Type=String,Description="Existing annotation">\n'
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            + extra_headers +
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\tNORMAL\n'
            + '\n'.join(records) + '\n')

    def record(self, pos=100, ref='G', alt='A', chrom='chr1', info='KEEP=original'):
        return '{}\t{}\t.\t{}\t{}\t60\tLowQual\t{}\tGT\t0/1\t0/0'.format(chrom, pos, ref, alt, info)

    def row(self, pos=100, start=100, ref='G', alt='A', value='0.9', chrom='1', sample='TUMOR'):
        return dict(Chromosome=chrom, vcf_pos=str(pos), Start_Position=str(start),
                    End_Position=str(start + 1 if ref == '-' else start + len(ref) - 1),
                    Reference_Allele=ref, Tumor_Seq_Allele2=alt, Tumor_Sample_Barcode=sample,
                    ccf_Mcopies=value)

    def write_maf(self, rows):
        with self.maf.open('w') as stream:
            stream.write('#version 2.4\n# FACETS MAF\n')
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)

    def run_annotation(self, *extra, defaults=False):
        args = [str(self.vcf), str(self.maf), '-o', str(self.output)]
        if not defaults:
            args += ['-m', self.mapping]
        with contextlib.redirect_stderr(io.StringIO()) as stderr:
            result = annotator.main(args + list(extra))
        self.assertIsNone(result)  # console_scripts must exit successfully
        with vcfpy.Reader.from_path(str(self.output)) as reader:
            return reader.header, list(reader), stderr.getvalue()

    def test_snp_and_mnp_at_same_coordinate_remain_distinct(self):
        self.write_vcf([self.record(), self.record(ref='GG', alt='AA')])
        self.write_maf([self.row(), self.row(ref='GG', alt='AA', value='0.4')])
        header, records, summary = self.run_annotation('--require-all-maf-matched')
        self.assertEqual([r.INFO['CCF'] for r in records], [[0.9], [0.4]])
        self.assertEqual(header.get_info_field_info('CCF').number, 'A')
        for record in records:
            self.assertEqual(record.FILTER, ['LowQual'])
            self.assertEqual(record.INFO['KEEP'], 'original')
            self.assertEqual(record.calls[0].data['GT'], '0/1')
        self.assertIn('matched_alleles=2', summary)

    def test_indels_and_prefix_trimmed_substitutions(self):
        self.write_vcf([self.record(100, 'ATG', 'A'), self.record(200, 'A', 'ATG'),
                        self.record(300, 'ACG', 'AT'), self.record(400, 'AG', 'AT')])
        self.write_maf([self.row(100, 101, 'TG', '-', '0.1'),
                        self.row(200, 200, '-', 'TG', '0.2'),
                        self.row(300, 301, 'CG', 'T', '0.3'),
                        self.row(400, 401, 'G', 'T', '0.4')])
        _, records, _ = self.run_annotation('--require-all-maf-matched')
        self.assertEqual([r.INFO['CCF'] for r in records], [[0.1], [0.2], [0.3], [0.4]])

    def test_second_alt_only_and_unmatched_record(self):
        self.write_vcf([self.record(alt='A,T'), self.record(200)])
        self.write_maf([self.row(alt='T')])
        _, records, summary = self.run_annotation()
        self.assertEqual(records[0].INFO['CCF'], [None, 0.9])
        self.assertNotIn('CCF', records[1].INFO)
        self.assertIn('unmatched_alleles=2', summary)

    def test_all_alt_values_keep_original_order(self):
        self.write_vcf([self.record(alt='A,T')])
        self.write_maf([self.row(alt='T', value='0.2'), self.row(alt='A', value='0.7')])
        _, records, _ = self.run_annotation()
        self.assertEqual(records[0].INFO['CCF'], [0.7, 0.2])

    def test_original_vcf_position_is_part_of_key(self):
        self.write_vcf([self.record(101, 'G', 'A')])
        self.write_maf([self.row(100, 101, 'G', 'A')])
        _, records, summary = self.run_annotation()
        self.assertNotIn('CCF', records[0].INFO)
        self.assertIn('unmatched_maf_rows=1', summary)

    def test_mitochondrial_alias_and_gzip_io(self):
        self.write_vcf([self.record(chrom='chrM')])
        self.write_maf([self.row(chrom='MT')])
        for field in ('vcf', 'maf'):
            path = getattr(self, field)
            compressed = path.with_suffix(path.suffix + '.gz')
            with gzip.open(compressed, 'wt') as stream:
                stream.write(path.read_text())
            setattr(self, field, compressed)
        self.output = self.output.with_suffix('.vcf.gz')
        _, records, _ = self.run_annotation()
        self.assertEqual(records[0].INFO['CCF'], [0.9])

    def test_duplicate_maf_keys_fail_before_output(self):
        self.write_vcf([self.record()])
        self.write_maf([self.row(), self.row(value='0.2')])
        with self.assertRaisesRegex(ValueError, 'Duplicate MAF'):
            self.run_annotation()
        self.assertFalse(self.output.exists())

    def test_duplicate_vcf_matches_fail_without_replacing_output(self):
        self.write_vcf([self.record(), self.record()])
        self.write_maf([self.row()])
        self.output.write_text('previous output')
        with self.assertRaisesRegex(ValueError, 'Multiple VCF alleles'):
            self.run_annotation()
        self.assertEqual(self.output.read_text(), 'previous output')
        self.assertFalse(list(self.root.glob('.vcf-maf-*')))

    def test_multiple_tumors_require_selection(self):
        self.write_vcf([self.record()])
        self.write_maf([self.row(), self.row(sample='OTHER', value='0.1')])
        with self.assertRaisesRegex(ValueError, 'multiple tumors'):
            self.run_annotation()
        _, records, _ = self.run_annotation('-s', 'TUMOR')
        self.assertEqual(records[0].INFO['CCF'], [0.9])
        with self.assertRaisesRegex(ValueError, 'not a sample in the VCF'):
            self.run_annotation('-s', 'OTHER')

    def test_header_collision_and_overwrite_clear_stale_unmatched_values(self):
        self.write_vcf([self.record(info='KEEP=original;CCF=0.1'),
                        self.record(200, info='KEEP=original;CCF=0.2')],
                       '##INFO=<ID=CCF,Number=1,Type=Float,Description="Old annotation">\n')
        self.write_maf([self.row()])
        with self.assertRaisesRegex(ValueError, 'already exist'):
            self.run_annotation()
        header, records, _ = self.run_annotation('--overwrite')
        self.assertEqual(header.get_info_field_info('CCF').number, 'A')
        self.assertEqual(records[0].INFO['CCF'], [0.9])
        self.assertNotIn('CCF', records[1].INFO)

    def test_missing_and_nonfinite_values(self):
        sentinels = ['NA', 'NaN', '.', '', 'Inf', '-Inf', '1e999', '0']
        self.write_vcf([self.record(100 + i) for i in range(len(sentinels))])
        self.write_maf([self.row(100+i, 100+i, value=value) for i, value in enumerate(sentinels)])
        _, records, _ = self.run_annotation()
        self.assertTrue(all('CCF' not in r.INFO for r in records[:-1]))
        self.assertEqual(records[-1].INFO['CCF'], [0.0])

    def test_all_22_defaults_include_cf_and_clonality(self):
        self.write_vcf([self.record()])
        row = self.row()
        for mapping in annotator.facets_mappings():
            row[mapping['source_col']] = ('CLONAL' if mapping['type'] == 'String' else
                                           '2' if mapping['type'] == 'Integer' else '0.8')
        self.write_maf([row])
        header, records, _ = self.run_annotation(defaults=True)
        self.assertEqual(len(annotator.facets_mappings()), 22)
        self.assertEqual(records[0].INFO['cf'], [0.8])
        self.assertEqual(records[0].INFO['clonality'], ['CLONAL'])
        self.assertEqual(records[0].INFO['tcn'], [2])
        self.assertEqual(header.get_info_field_info('tcn').type, 'Integer')

    def test_strict_matching_leaves_no_partial_output(self):
        self.write_vcf([self.record()])
        self.write_maf([self.row(), self.row(200, 200)])
        with self.assertRaisesRegex(ValueError, '1 MAF rows have no VCF allele match'):
            self.run_annotation('--require-all-maf-matched')
        self.assertFalse(self.output.exists())

    def test_missing_required_column_and_bad_coordinates(self):
        self.write_vcf([self.record()])
        row = self.row()
        del row['vcf_pos']
        self.write_maf([row])
        with self.assertRaisesRegex(ValueError, 'missing columns: vcf_pos'):
            self.run_annotation()
        row = self.row()
        row['End_Position'] = '99'
        self.write_maf([row])
        with self.assertRaisesRegex(ValueError, 'Inconsistent vcf2maf coordinates'):
            self.run_annotation()

    def test_symbolic_vcf_alleles_preserved_unannotated(self):
        self.write_vcf([self.record(), self.record(200, alt='<DEL>'), self.record(300, alt='*')])
        self.write_maf([self.row()])
        _, records, _ = self.run_annotation()
        self.assertEqual(records[1].ALT[0].serialize(), '<DEL>')
        self.assertNotIn('CCF', records[1].INFO)
        self.assertNotIn('CCF', records[2].INFO)

    def test_output_cannot_replace_input(self):
        self.write_vcf([self.record()])
        self.write_maf([self.row()])
        original = self.vcf.read_bytes()
        self.output = self.vcf
        with self.assertRaises(SystemExit):
            self.run_annotation()
        self.assertEqual(self.vcf.read_bytes(), original)


if __name__ == '__main__':
    unittest.main()
