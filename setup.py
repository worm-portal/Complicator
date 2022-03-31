import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="Complicator",
    version="0.0.1",
    author="Grayson Boyer",
    author_email="gmboyer@asu.edu",
    description="Estimate thermodynamic properties of aqueous inorganic complexes with monovalent ligands.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    entry_points={},
    packages=['Complicator'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8',
    install_requires=['pandas'],
    include_package_data=True,
    package_data={'': ['*.txt']},
    zip_safe=False
)

