from os import path
from setuptools import setup

here = path.abspath(path.dirname(__file__))

with open(path.join(here, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

REQUIRES = ["biopython==1.76"]

setup(
    name="phykit",
    description="",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Jacob L. Steenwyk",
    author_email="jlsteenwyk@gmail.com",
    url="https://github.com/jlsteenwyk/phykit",
    packages=["phykit"],
    entry_points={"console_scripts": ["phykit = phykit.phykit:main"]},
    version="0.0.1",
    include_package_data=True,
    install_requires=REQUIRES,
)
