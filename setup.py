
import io
from setuptools import setup

setup(name='vcfsyncer',
    description='Package to synchronize multiple VCFs',
    long_description=io.open('README.md', encoding='utf-8').read(),
    long_description_content_type='text/markdown',
    version='1.0.1',
    author='Jeremy McRae',
    author_email='jmcrae@illumina.com',
    license='MIT',
    url='https://github.com/jeremy_mcrae/vcfsyncer',
    packages=['vcfsyncer'],
    install_requires=['pysam'],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ])
