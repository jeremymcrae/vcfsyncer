
import pysam

from vcfsyncer.variant import MultiVariantRecord

CHROMS = list(map(str, range(1, 23))) + ['X', 'Y', 'MT']

class VCFSyncer:
    """ this synchronises the coordinates of multiple VCFs
    
    I want to be able to take multiple VCFs and iterate through the unique set
    of variants in all VCFs, but keep the VCFs synchronised, despite the VCFs
    potentially containing different variants.
    """
    def __init__(self, *args, chromosomes=None):
        """ initialize class with VCFs and defined chromosomes
        
        Args:
            args: list of VCF paths, or pysam.VariantFile objects
            chromosomes: list of chromosomes (sorted correctly) that we want to
                consider variants from. Without this the VCFs can desynchronize
                if they contain calls on non-standard chromosomes.
        """
        if chromosomes is None:
            chromosomes = CHROMS
        
        # ensure VCFs sort by chrom correctly. Strip any 'chr' prefix, then
        # convert to dict, chrom map to list position, and finally add another
        # set of chromosomes with chr prefixes, but with same sort order.
        chromosomes = [x.lstrip('chr') for x in chromosomes]
        self.chroms = dict(zip(chromosomes, range(len(chromosomes))))
        for k in list(self.chroms):
            self.chroms['chr' + k] = self.chroms[k]
        
        self.vcfs = list(args)
        if type(self.vcfs[0]) == str:
            self.vcfs = [pysam.VariantFile(str(x)) for x in self.vcfs]
        
        # add header field, so we can store alts within INFO dicts from each VCF
        for x in self.vcfs:
            x.header.add_line('##INFO=<ID=alts,Number=A,Type=String,Description="Alternate alleles">')
        
        self.variants = [ next(x) for x in self.vcfs ]
        self.leading = self._first()
        self.cache = {'chrom': None, 'pos': None, 'refs': {}}
    
    def __repr__(self):
        vcfs = ', '.join(map(str, self.vcfs))
        return 'VCFSyncer({})'.format(self.vcfs)
    
    def _match(self, a, b):
        """ check if two variants are at the same site and share a ref allele
        """
        if b is None:
            return False
        
        return a.chrom == b.chrom and a.pos == b.pos and a.ref == b.ref
    
    def _int(self, chrom):
        ''' convert a chromosome to an integer for correct sorting
        '''
        try:
            return self.chroms[chrom]
        except KeyError:
            # allow for non-standard contigs
            return abs(hash(chrom))
    
    def _first(self):
        """ find which of the vcfs have a variant at the lowest genome position
        
        I'm assuming the VCFs are sorted, and we want to use variants which
        come first in the genome (ie are on an earlier chromosome and have a
        lower nucleotide position).
        
        Returns:
            list of booleans showing if each variant is at the lowest position
        """
        vars = ( x for x in self.variants if x is not None )
        var = min(vars, key=lambda x: (self._int(x.chrom), x.pos))
        return [ self._match(var, x) for x in self.variants ]
    
    def __iter__(self):
        return self
    
    def __next__(self):
        """ iterate through vcfs, returning lowest site variants at each step
        """
        if all( x is None for x in self.vcfs ):
            raise StopIteration
        
        refs = sorted(self.cache['refs'])
        if len(refs) == 0:
            self._fill_cache()
            refs = sorted(self.cache['refs'])
        
        ref = refs[0]
        vars = self.cache['refs'][ref]
        del self.cache['refs'][ref]
        
        return MultiVariantRecord(*vars)
    
    def _increment(self):
        """ increment the vcfs/variants ranked first (so iteration works)
        """
        for i, in_lead in enumerate(self.leading):
            if in_lead and self.vcfs[i] is not None:
                try:
                    self.variants[i] = next(self.vcfs[i])
                except StopIteration:
                    self.vcfs[i] = None
                    self.variants[i] = None
            
            # avoid variants on non-standard chromosomes, otherwise we easily
            # get out of sync
            while self.variants[i] is not None and self.variants[i].chrom not in self.chroms:
                try:
                    self.variants[i] = next(self.vcfs[i])
                except StopIteration:
                    self.vcfs[i] = None
                    self.variants[i] = None
        
        try:
            self.leading = self._first()
        except ValueError:
            self.leading = []
    
    def fetch(self, contig=None, start=None, stop=None, region=None, reopen=False, end=None, reference=None):
        ''' imitate pysam fetch behavior
        '''
        
        iterators = (x.fetch(contig, start, stop, region, reopen, end, reference) for x in self.vcfs)
        self.variants = [ next(x) for x in iterators ]
        self.leading = self._first()
        
        for x in self:
            if x.chrom != contig or stop is not None and x.pos > stop:
                break
            yield x
    
    def _fill_cache(self):
        """ fill cache with variants at the same site
        """
        
        if len(self.cache['refs']) == 0:
            var = [ x for x, y in zip(self.variants, self.leading) if y ][0]
            self.cache['chrom'] = var.chrom
            self.cache['pos'] = var.pos
        
        while True:
            vars = [ x for x, y in zip(self.variants, self.leading) if y ]
            if not any( x.pos == self.cache['pos'] for x in vars):
                break
            
            for x in vars:
                if x.ref not in self.cache['refs']:
                    self.cache['refs'][x.ref] = []
                
                self.cache['refs'][x.ref].append(x)
            
            self._increment()
