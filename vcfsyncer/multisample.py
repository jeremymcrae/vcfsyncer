
from warnings import warn

class MultiSamples(object):
    def __init__(self, *args):
        self.samples = {}
        self.idx = {}
        self.merge(*args)
    
    def merge(self, *samples):
        for group in samples:
            for sample, values in group.items():
                if sample in self.samples:
                    warn('duplicate sample across VCFs: {}'.format(sample))
                    continue
                
                self.idx[len(self.samples)] = sample
                self.samples[sample] = values
    
    def __len__(self):
        return len(self.samples)
    
    def __getitem__(self, key):
        if isinstance(key, int):
            assert 0 <= key < len(self)
            return self.samples[self.idx[key]]
        if key in self:
            return self.samples[key]
        
    def __iter__(self):
        for x in self.samples:
            yield x
    
    def __contains__(self, key):
        return key in self.samples or 0 <= key < len(self.samples)
    
    def keys(self):
        return self.samples.keys()
        
    def values(self):
        return self.samples.values()
        
    def pop(self, key):
        value = self[key]
        del self[key]
        return value
