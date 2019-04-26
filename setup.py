# -*- coding: utf-8 -*-
from setuptools import setup

__author__ = "Jason Eu"
__license__ = "MIT"

about = {}
with open('sagar/version.py') as f:
    exec(f.read(), about)

setup(
    name = 'sagar',
    # Remenber change version in ababe.cmdline.runabalib.exec_from_cmdline()
    version = about['__version__'],
    description='Structures of Alloys Generation And Recognition',
    long_description=open('README.rst').read(),
    author='Jason Eu',
    author_email='morty.yu@yahoo.com',
    license=__license__,
    include_package_data=True,
    package_data={},
    keywords = 'crystal material strucutre DFT',
    install_requires=[
        "numpy",
        "spglib",
        "click>=6",
        "pniggli",
    ],
    extras_require={
        'dev': [
            'pip',
            'pytest',
            'pytest-cov',
            'twine',
            'ipython',
        ],
        'ase': [
            "ase==3.16",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 2.7",
        "Operating System :: POSIX :: Linux",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
        "Natural Language :: Chinese (Simplified)",
        "Natural Language :: English"],
    entry_points={
        'console_scripts': [
            'rexpand = sagar.rexpand:cli',
        ],
    },
    packages=[
        'sagar',
        'sagar.crystal',
        'sagar.element',
        'sagar.io',
        'sagar.molecule',
        'sagar.toolkit',
    ],
    test_suite='test',
)
