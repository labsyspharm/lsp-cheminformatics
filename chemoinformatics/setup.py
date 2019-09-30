import sys

from setuptools import setup


py_version_requires = {2: ["chemfp"], 3: ["marshmallow", "apispec"]}

setup(
    name="tas_chemoinformatics",
    version="0.2",
    packages=["tas_chemoinformatics"],
    include_package_data=True,
    install_requires=["click", "flask", "molvs", "rdkit", "pandas", "gunicorn"]
    + py_version_requires[sys.version_info[0]],
    classifiers=[
        "Programming Language :: Python",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux  ",
    ],
    entry_points="""
        [console_scripts]
        tautomers=tas_chemoinformatics.tautomers_cli:cli
    """,
)
