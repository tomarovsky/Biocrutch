__author__ = 'tomarovsky'

import os
from pathlib import Path
from setuptools import setup, find_packages

dependencies = ['scipy', 'numpy', 'pandas', 'matplotlib', 'matplotlib-venn',
                'venn', 'biopython', 'xmltodict', 'statsmodels', "ete3",
                "bcbio-gff", 'lxml', 'bs4', 'six']

setup(name='Biocrutch',
      version='0.1',
      packages=find_packages(),
      author='Andrey Tomarovsky',
      author_email='andrey.tomarovsky@gmail.com',
      install_requires=dependencies,
      long_description=open(os.path.join(os.path.dirname(__file__), 'README.md')).read(),
      scripts=list(map(str, sorted(Path('scripts/').rglob("*.py")))))