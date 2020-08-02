from setuptools import setup, find_packages

with open("README.md", "r") as f:
    long_desc = f.read()

setup(
    name="autochem",
    version="0.1.0",
    description="Automates creation and post-processing of quantum chemical calculations",
    packages=find_packages(),
    py_modules=['autochem'],
    python_requires=">=3.6",
    scripts=["bin/chem_assist"],
    classifiers=[
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Natural Language :: English",
    ],
    long_description=long_desc,
    long_description_content_type="text/markdown",
    url="https://github.com/tommason14/autochem",
    author="Tom Mason",
    author_email="tom.mason14+pypi@gmail.com",
    install_requires=["pandas >= 1.0.1", "numpy >= 1.18.2", "dfply >= 0.3.3"],
)
