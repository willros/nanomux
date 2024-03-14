import setuptools
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()


setuptools.setup(
    name="nanomux",
    version="0.1.0",
    description="Demux reads from Nanopore sequencing",
    url="https://github.com/willros/nanomux",
    author="William Rosenbaum",
    author_email="william.rosenbaum@gmail.com",
    license="MIT",
    packages=setuptools.find_packages(),
    long_description=long_description,
    long_description_content_type="text/markdown",
    python_requires=">=3.9",
    install_requires=[
        "pandas==2.2.1",
        "regex==2023.12.25",
        "pyfastx==2.0.2",
    ],
    entry_points={"console_scripts": ["nanomux=nanomux.main:cli"]},
)
