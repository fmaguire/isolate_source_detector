#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages
from isolate_source_detector import __version__

with open('README.md') as readme_file:
    readme = readme_file.read()

requirements = ['Click>=7.0', "Biopython", "tqdm", "ete3", "pandas"]

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest>=3', ]

setup(
    author="Finlay Maguire",
    author_email='finlaymaguire@gmail.com',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
    description="Script to extract likely isolate sources from augur trees and nextstrain metadata",
    entry_points={
        'console_scripts': [
            'isd=isolate_source_detector.cli:main',
        ],
    },
    install_requires=requirements,
    license="GNU General Public License v3",
    include_package_data=True,
    keywords='isolate_source_detector',
    name='isolate_source_detector',
    packages=find_packages(include=['isolate_source_detector', 'isolate_source_detector.*']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/fmaguire/isolate_source_detector',
    version=__version__,
    zip_safe=False,
)
