from setuptools import setup

setup(
    name="tas_chemoinformatics",
    version="0.1",
    packages=["tas_chemoinformatics"],
    include_package_data=True,
    install_requires=["click", "flask", "molvs"],
    entry_points="""
        [console_scripts]
        tautomers=tas_chemoinformatics.tautomers_cli:cli
    """,
)
