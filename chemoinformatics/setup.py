from setuptools import setup

setup(
    name="tas_chemoinformatics",
    version="0.1",
    py_modules=["tautomers", "fingerprints"],
    include_package_data=True,
    install_requires=["click", "rdkit", "flask", "molvs"],
    entry_points="""
        [console_scripts]
        tautomers=tautomers:cli,
        fingerprints=fingerprints:cli,
    """,
)
