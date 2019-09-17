
from pkg_resources import resource_filename
import unittest

from vcfsyncer.vcf import VCFSyncer

class TestVCFSyncer(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        paths = ['a.vcf.gz', 'b.vcf.gz', 'c.vcf.gz', 'no_samples.vcf.gz',
            'unknown_chrom.vcf.gz', 'no_alts.vcf.gz']
        labels = [x.split('.')[0] for x in paths]
        paths = [resource_filename(__name__, 'data/' + x) for x in paths]
        self.paths = dict(zip(labels, paths))
    
    def test__match(self):
        ''' test that VCFSyncer._match works correctly
        '''
        
        class Var:
            def __init__(self, chrom, pos, ref):
                self.chrom = chrom
                self.pos = pos
                self.ref = ref
        
        a = Var('X', 1000, 'A')
        
        vcf = VCFSyncer(self.paths['a'])
        self.assertFalse(vcf._match(a, Var('X', 1001, 'A')))
        
        # variant at the same site with same ref should pass
        self.assertTrue(vcf._match(a, Var('X', 1000, 'A')))
        
        # variant at the same site with different ref doesn't match
        self.assertFalse(vcf._match(a, Var('X', 1000, 'T')))
        
        # variant at the same site with different ref doesn't match
        self.assertFalse(vcf._match(a, Var('1', 1000, 'A')))
    
    def test__int(self):
        ''' test that VCFSyncer._int works correctly
        '''
        vcf = VCFSyncer(self.paths['a'])
        self.assertEqual(vcf._int('1'), 0)
        self.assertEqual(vcf._int('5'), 4)
        self.assertEqual(vcf._int('X'), 22)
        self.assertEqual(vcf._int('chrX'), 22)
        
        # unknown chromosomes get a much higher sort order
        self.assertTrue(vcf._int('ZZZ') > 22)
        
        # check user defined chromosome order is respected
        vcf = VCFSyncer(self.paths['a'], chromosomes=['5', '4'])
        self.assertEqual(vcf._int('5'), 0)
        self.assertEqual(vcf._int('4'), 1)
    
    def test__first(self):
        ''' test that VCFSyncer._first works correctly
        '''
        vcf = VCFSyncer(self.paths['a'], self.paths['c'])
        self.assertEqual(vcf._first(), [True, False])
        
        # now load three files where two files are at the leading edge
        vcf = VCFSyncer(self.paths['a'], self.paths['b'], self.paths['c'])
        self.assertEqual(vcf._first(), [True, True, False])
    
    def test___next__(self):
        ''' test that VCFSyncer.__next__ works correctly
        '''
        vcf = VCFSyncer(self.paths['a'], self.paths['c'])
        
        # only the first VCF should have a sample at the first site. Note that
        # this doesn't extensively test the variant output, those tests are in
        # the MultiVariantRecord tests.
        var = next(vcf)
        self.assertEqual(var.chrom, '1')
        self.assertEqual(var.pos, 25000000)
        self.assertEqual(var.ref, 'A')
        self.assertEqual(list(var.samples), ['sampleA'])
        
        # jump forward two variant, where both VCFs have data for a variant
        var = next(vcf)
        var = next(vcf)
        self.assertEqual(var.chrom, '2')
        self.assertEqual(var.pos, 179415988)
        self.assertEqual(var.ref, 'C')
        self.assertEqual(list(var.samples), ['sampleA', 'sampleC', 'sampleD'])
    
    def test___next___unknown_chrom(self):
        ''' test that VCFSyncer.__next__ works correctly with unknown chroms
        '''
        vcf = VCFSyncer(self.paths['a'], self.paths['unknown_chrom'])
        var = next(vcf)
        self.assertEqual(var.chrom, '1')
        self.assertEqual(var.pos, 25000000)
        self.assertEqual(var.ref, 'A')
        self.assertEqual(list(var.samples), ['sampleA', 'sampleU'])
    
    def test_fetch(self):
        ''' test that VCFSyncer.fetch works correctly
        '''
        vcf = VCFSyncer(self.paths['a'], self.paths['c'])
        vars = vcf.fetch('2', 179415987, 179446000)
        
        self.assertEqual([x.pos for x in vars], [179415988])
        
        # check across a wider position range
        vars = vcf.fetch('2', 179415987, 179447000)
        self.assertEqual([x.pos for x in vars], [179415988, 179446218])
    
    def test_no_alt_allele(self):
        ''' test that we can load VCFs with sites without ALT alleles
        
        These lines can occur when we want to indicate that a genomic position
        (or range) matches the reference genome.
        '''
        vcf = VCFSyncer(self.paths['no_alts'], self.paths['a'])
        var = next(vcf)
        self.assertEqual(var.pos, 10337)
