from setuptools import setup, find_packages
from distutils.command.install import install
import os
import subprocess
from glob import glob

with open("README.md", "r") as f:
    long_desc = f.read()


class CustomInstall(install):
    def run(self):
        install.run(self)
        compile_thermo()


def compile_thermo():
    print("Compiling fortran")
    cwd = os.getcwd()
    os.chdir("autochem/core")
    ret = subprocess.call("gfortran thermo.f -o thermo.exe", shell=True)
    if ret != 0:
        sys.exit("Fortran compilation failed")
    os.chdir(cwd)


setup(
    name="autochem",
    version="0.1.6",
    description="Automates creation and post-processing of quantum chemical calculations",
    packages=find_packages(),
    cmdclass={"install": CustomInstall},
    data_files=[
        ("autochem/core", "autochem/core/thermo.f"),
        ("autochem/templates", glob("autochem/templates/*")),
    ],
    py_modules=["autochem"],
    python_requires=">=3.6",
    scripts=["bin/autochem"],
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

compile_thermo()
