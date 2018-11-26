from pkg_resources import get_distribution

name = 'vcfsyncer'
__version__ = get_distribution(name).version

from vcfsyncer.vcf_syncer import VCFSyncer
