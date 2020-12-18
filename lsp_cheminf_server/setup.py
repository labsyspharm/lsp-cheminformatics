import sys

from setuptools import setup

setup(
    name="lspcheminf",
    version="0.4.1",
    packages=["lspcheminf"],
    include_package_data=True,
    install_requires=[
        "click",
        "flask",
        "molvs",
        "rdkit",
        "pandas",
        "gunicorn",
        "marshmallow",
        "apispec",
        "chembl_structure_pipeline@ssh://git@github.com/chembl/ChEMBL_Structure_Pipeline.git",
        "pebble",
    ],
    classifiers=[
        "Programming Language :: Python",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux  ",
    ],
    entry_points="""
        [console_scripts]
        tautomers=lspcheminf.tautomers_cli:cli
    """,
)
