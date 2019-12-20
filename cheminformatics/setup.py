import sys

from setuptools import setup

setup(
    name="tas_cheminformatics",
    version="0.2",
    packages=["tas_cheminformatics"],
    include_package_data=True,
    install_requires=["click", "flask", "molvs", "rdkit", "pandas", "gunicorn", "marshmallow", "apispec"],
    classifiers=[
        "Programming Language :: Python",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux  ",
    ],
    entry_points="""
        [console_scripts]
        tautomers=tas_cheminformatics.tautomers_cli:cli
    """,
)
