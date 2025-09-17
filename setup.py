from setuptools import setup, find_packages
import os

# Read requirements from requirements.txt
with open("requirements.txt") as f:
    required = f.read().splitlines()


# Read long description from README.md
with open("README.md", encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="v4c",
    version="0.1.0",
    packages=find_packages(),
    install_requires=required,
    entry_points={
        "console_scripts": [
            "v4c-extract = v4c.extract:main",
            "v4c-plot = v4c.plot:main",
            "v4c-compare = v4c.compare:main",
            "v4c-convert = v4c.convert:main",
        ],
    },
    python_requires=">=3.6",
    author="Zoe Chen",
    author_email="zoe.chen2@ucsf.edu",
    description="Virtual 4C Analysis for Hi-C Data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/zoechen0717/V4C",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    project_urls={
        "Bug Reports": "https://github.com/zoechen0717/V4C/issues",
        "Source": "https://github.com/zoechen0717/V4C",
        "Documentation": "https://v4c.readthedocs.io/",
    },
    include_package_data=True,
    package_data={
        "v4c": ["genome/*.bed"],
    },
)
