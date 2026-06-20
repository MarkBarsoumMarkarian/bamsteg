from setuptools import setup, find_packages

setup(
    name="bamsteg",
    version="0.1.0",
    packages=find_packages(),
    install_requires=["pysam", "pycryptodome", "reedsolo"],
    entry_points={"console_scripts": ["bamsteg=bamsteg.cli:main"]},
)
