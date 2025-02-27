from setuptools import setup, find_packages

setup(
    name="v4c",
    version="0.1.0",
    author="Zoe Chen",
    author_email="zoe.chen2@ucsf.edu",
    description="Virtual 4C analysis package for Hi-C contact frequency extraction and visualization.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/zoechen0717/V4C",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "pandas",
        "matplotlib",
        "cooler",
        "argparse"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    entry_points={
        "console_scripts": [
            "v4c-extract = v4c.extract:main",
            "v4c-plot = v4c.plot:main",
            "v4c-compare = v4c.compare:main"
        ]
    },
)
