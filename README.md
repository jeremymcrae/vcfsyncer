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
for group in synced:
    for match in group:
        print(match)

# alternatively, fetch regions from the synchronized VCFs (if tabix-indexed)
for group in synced.fetch(chrom='1', start=1000000, end=2000000):
    for match in group:
        print(match)
```
