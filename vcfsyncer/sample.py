
from itertools import chain
from warnings import warn

class MultiSamples(object):
    def __init__(self, *samples):
        
        self.samples = samples
        self.idx = {}
        self._index_groups()
    
    def _index_groups(self):
        ''' identify which VCF each sample is in
        
        Allows quick data access by sample ID, without changing data around.
        '''
        
        for i, group in enumerate(self.samples):
            for sample in group:
                if sample in self.idx:
                    warn('duplicate sample across VCFs: {}'.format(sample))
                
                self.idx[sample] = i
    
    def __len__(self):
        return sum(len(x for x in self.samples))
    
    def __getitem__(self, key):
        ''' access sample data either by sample ID (str), or index offset (int)
        '''
        if isinstance(key, int):
            assert 0 <= key < len(self)
            for group in self.samples:
                if key < len(group):
                    return group[key]
                key -= len(group)
        
        if key in self:
            group = self.idx[key]
            return self.samples[group][key]
    
    def __iter__(self):
        for sample_id in self.idx:
            yield sample_id
    
    def __contains__(self, key):
        return key in self.idx or 0 <= key < len(self)
    
    def keys(self):
        return self.idx
        
    def values(self):
        return chain(x.values() for x in self.samples)
        
    def pop(self, key):
        value = self[key]
        del self[key]
        return value
