import os
import glob
import setuptools
from distutils.core import setup

with open("README.md", 'r') as readme:
    long_description = readme.read()

setup(
    name='vivarium-bioscrape',
    version='0.0.0.6',
    packages=[
        'vivarium_bioscrape',
        'vivarium_bioscrape.processes',
        'vivarium_bioscrape.composites',
        'vivarium_bioscrape.library'
    ],
    author='William Poole, Eran Agmon, Ryan Spangler',
    author_email='',
    url='https://github.com/vivarium-collective/vivarium-bioscrape',
    license='MIT',
    entry_points={
        'console_scripts': []},
    short_description='Integrates Vivarium with the Bioscrape CRN simulator, allowing Chemical Reaction Network Models to be Vivarium Processes.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    package_data={},
    include_package_data=True,
    install_requires=[
        'vivarium-core>=0.2.0',
        'bioscrape'
    ],
)
