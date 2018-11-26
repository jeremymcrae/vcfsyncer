## vcfsyncer: A package to synchronize multiple VCFs
This is a way to transparently use single sample VCFs as if they were a
multi-sample VCF.

### Installation
The simplest way to install vcfsyncer is through pip:
```sh
pip install vcfsyncer
```

### Usage
Use vcfsyncer within a python environment
```python

from vcfsyncer import VCFSyncer

synced = VCFSyncer(PATH_A, PATH_B, PATH_C)

# you can iterate through the synchronized VCFs
for var in synced:
    print(var.chrom, var.pos)

# alternatively, fetch regions from the synchronized VCFs (if tabix-indexed)
for var in synced.fetch(chrom='1', start=1000000, end=2000000):
    print(var.chrom, var.pos)
```
