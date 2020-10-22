from os import path
from setuptools import setup, find_packages

here = path.abspath(path.dirname(__file__))

with open(path.join(here, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

about = {}
with open(path.join(here, 'phykit', 'version.py'), encoding="utf-8") as f:
    exec(f.read(), about)

CLASSIFIERS = [
    'Operating System :: OS Independent',
    'Intended Audience :: Science/Research',
    'Programming Language :: Python',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
    'Topic :: Scientific/Engineering',
]

REQUIRES = ["biopython==1.76", "scipy==1.4.1"]

setup(
    name="phykit",
    description="",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Jacob L. Steenwyk",
    author_email="jlsteenwyk@gmail.com",
    url="https://github.com/jlsteenwyk/phykit",
    packages=find_packages(),
    classifiers=CLASSIFIERS,
    entry_points={"console_scripts": ["phykit = phykit.phykit:main"]},
    version=about["__version__"],
    include_package_data=True,
    install_requires=REQUIRES,
)
