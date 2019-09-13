
import unittest
from collections import namedtuple

from vcfsyncer.variant import MultiVariantRecord

Var = namedtuple('Var', ['chrom', 'pos', 'id', 'ref', 'alts', 'qual', 'filter',
    'info', 'format', 'samples'])

class TestVCFVariant(unittest.TestCase):
    def test_missing_alt(self):
        x = Var('1', 100, '.', 'A', ['G'], 100, 'PASS', {}, 'AD', {'A': '10,10'})
        y = Var('1', 100, '.', 'A', ['G'], 100, 'PASS', {}, 'AD', {'B': '10,10'})
        z = Var('1', 100, '.', 'A', ['C'], 100, 'PASS', {}, 'AD', {'C': '10,10'})
        m = Var('1', 100, '.', 'A', None, 100, 'PASS', {}, None, {})
        multivar = MultiVariantRecord(x, y, z, m)
        self.assertEqual(multivar.alts, ['G', 'C'])
