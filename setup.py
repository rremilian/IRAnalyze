from setuptools import setup, find_packages
from IRAnalyze import iranalyze

setup(
    name = "IRAnalyze",
    version = iranalyze.iranalyze.version,
    url = "https://github.com/rremilian/IRAnalyze",
    author = "rremilian",
    description = "A Python module to analyze experimental and theoretical spectra",
    packages=find_packages(),
    install_requires= ["numpy >= 1.21.5", "matplotlib >= 3.5.1", "scipy >= 1.7.3"],
)