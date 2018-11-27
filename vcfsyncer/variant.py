
from vcfsyncer.sample import MultiSamples

class MultiVariantRecord(object):
    ''' represent variants from different VCFs, at same site, as single object
    
    You'll need to be careful about accessing the attributes afterwards, since
    the chrom, pos and ref attributes are valid, whereas alts is now a
    deduplicated list of alts from the variant records, info is a list of info
    dictionaries, and samples is a specialised class which mimics the pysam API.
    
    Possibly could be better integrated?
    '''
    def __init__(self, *variants):
        self.chrom = variants[0].chrom
        self.contig = self.chrom
        self.pos = variants[0].pos
        self.ref = variants[0].ref
        self.alts = variants
        self.filter = [x.filter for x in variants]
        self.format = [x.format for x in variants]
        self.info = [x.info for x in variants]
        self.samples = MultiSamples(*[x.samples for x in variants])
        self.id = [x.id for x in variants]
        self.qual = [x.qual for x in variants]
    
    @property
    def alts(self):
        return self._alts
    @alts.setter
    def alts(self, variants):
        self._alts = (alt for var in variants for alt in var.alts)
        seen = set()
        self._alts = [x for x in self._alts if not (x in seen or seen.add(x))]
    
    @property
    def alleles(self):
        return [self.ref] + self.alts
