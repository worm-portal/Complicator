import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="Complicator",
    version="0.2.1",
    author="Grayson Boyer",
    author_email="gmboyer@asu.edu",
    description="Estimate thermodynamic properties of aqueous metal complexes with monovalent ligands.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    entry_points={},
    packages=['Complicator'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.10',
    install_requires=['pandas', 'WORMutils', 'chemparse'],
    include_package_data=True,
    package_data={'': ['*.txt', '*.csv']},
    zip_safe=False
)

