"""
Python-packaging for synbiomts

Copyright 2017 Alexander C. Reis, Howard M. Salis, all rights reserved.

"""

from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='synbiomts',
      version='1.0',
      description='Test suite for DNA sequence-function models',
      url='http://github.com/reisalex/SynBioMTS',
      author='Alexander C. Reis',
      author_email='alex.reis@psu.edu',
      license='MIT',
      packages=['synbiomts'],
      # install_requires=['numpy','scipy','pandas','biopython'],
      zip_safe=False)