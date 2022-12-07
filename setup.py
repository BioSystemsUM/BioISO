from setuptools import setup, find_packages

setup(
    name='bioiso',
    version='0.0.1',
    package_dir={'': 'src'},
    install_requires=["cobra",
                      "numpy",
                      "pandas"],

    author='Fernando Cruz',
    author_email='fernando.cruz@ceb.uminho.pt',
    description='BioISO - Biomass constraint-based In Silico Optimization',
    license='GNU General Public License v3.0',
    keywords='metabolic model reconstruction and analysis',
    url='https://github.com/BioSystemsUM/BioISO',
    long_description=open('README.md').read(),
    classifiers=[
        'Development Status :: 1 - Beta',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3.7',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Topic :: Software Development :: Libraries :: Python Modules'
    ],
)
