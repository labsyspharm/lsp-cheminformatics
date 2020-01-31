import sys

from setuptools import setup

setup(
    name="lspcheminf",
    version="0.3",
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
