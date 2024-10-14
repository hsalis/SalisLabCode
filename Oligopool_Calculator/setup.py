from setuptools import setup
from os import path
from io import open

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

# Arguments marked as "Required" below must be included for upload to PyPI.
# Fields marked as "Optional" may be commented out.

setup(
    name='oligopool',
    
    # Link: https://www.python.org/dev/peps/pep-0440/#version-scheme
    version='0.0.0',
    
    description='Oligopool Calculator - Automated design and analysis of oligopools for massively parallel assays',
    
    long_description=long_description,
    
    long_description_content_type='text/markdown',
    
    url='https://github.com/ayaanhossain/oligopool',
    
    author='Ayaan Hossain and Howard Salis',
    
    author_email='auh57@psu.edu, salis@psu.edu',  # Optional
    
    classifiers=[  # Optional
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',

        # Pick your license as you wish
        'License :: OSI Approved :: GPL Version 3',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        # These classifiers are *not* checked by 'pip install'. See instead
        # 'python_requires' below.
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
    keywords=' '.join([
        'synthetic',
        'computational',
        'biology',
        'genetic',
        'parts',
        'calculator',
        'non-repetitive',
        'design',
        'discovery',
        'algorithm',
        'stable',
        'systems',
        'nrp',
        'repeats',
        'vertex',
        'cover',
        'path',
        'finding']),

    packages=['oligopool'],

    package_dir={
        'oligopool': './oligopool'
    },

    python_requires=', '.join([
        '!=2.7',
        '!=3.0.*',
        '!=3.1.*',
        '!=3.2.*',
        '!=3.3.*',
        '!=3.4.*',
        '!=3.5.*',
        '>=3.6.*',
        '<4.0.*']),

    install_requires=[
        'numpy>=1.19.0',
        'biopython>=1.77',
        'leveldb>=0.201',
        'scipy>=1.5.1',
        'networkx>=2.4',
        'nrpcalc',
        'jupyter>=1.0.0',
        'scikit-learn>=0.23.1',
        'seaborn>=0.10.1',
        'statsmodels>=0.11.0'
        ],

    project_urls={  # Optional
        'Bug Reports': 'https://github.com/ayaanhossain/oligopool/issues',
        'Source'     : 'https://github.com/ayaanhossain/oligopool/tree/master/oligopool',
    },
)
