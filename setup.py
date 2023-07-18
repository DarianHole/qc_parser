'''
Setup.py file for the qc_parser package.
'''
import setuptools
from qc_parser import __version__, _package_name

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name=_package_name,
    version=__version__,
    author="Darian Hole, Original: Richard J. de Borja",
    author_email="darian.hole@phac-aspc.gc.ca",
    description="A package for parsing viral qc analysis files",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/DarianHole/qc_parser",
    packages=setuptools.find_packages(),
    license="MIT",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.9',
    scripts=['bin/qc_parser']
)
