from setuptools import setup, find_packages

setup(
    name = 'bioiso',
    version = '0.0.1',
    package_dir = {'':'bioiso'},
    packages = find_packages('bioiso'),
    install_requires = ["cobra",
                        "numpy",
                        "pandas"],

    author = 'Fernando Cruz',
    author_email = 'fernando.cruz@ceb.uminho.pt',
    description = 'BioISO - Biomass constraint-based In Silico Optimization',
    license = 'GNU General Public License v3.0',
    keywords = 'metabolic model reconstruction and analysis',
    url = 'https://github.com/BioSystemsUM/bioiso',
    long_description = open('README.rst').read(),
    classifiers = [
        'Development Status :: 1 - Beta',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3.7',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Topic :: Software Development :: Libraries :: Python Modules'
    ],
)