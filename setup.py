from setuptools import setup, find_packages
from os import path as path

with open("Readme.md") as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

exclusions=""
if (path.isdir('CMakeFiles')):
    exclusions="CMakeFiles"

setup(
    name = 'mzMLTrace',
    version = '0.1.0',
    description = 'python package to extract Mass Traces from mzML files',
    long_description = readme,
    author = "Joseph O'Brien",
    author_email="joe@xialab.ca",
    url = "https://github.com/obrien951/mzMLTrace.git",
    license = license,
    packages = find_packages(where="src",exclude=(exclusions)),
    package_dir={"":"src"},
    install_requires=[
    'pybind11',
    ],
)
