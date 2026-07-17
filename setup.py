from os import path
from setuptools import setup, find_packages

from phykit.version import __version__
from phykit.cli_registry import PUBLIC_COMMAND_TO_HANDLER

here = path.abspath(path.dirname(__file__))

with open(path.join(here, "README.md")) as f:
    long_description = f.read()

CLASSIFIERS = [
    'Operating System :: OS Independent',
    'Intended Audience :: Science/Research',
    'Programming Language :: Python',
    'Programming Language :: Python :: 3.10',
    'Programming Language :: Python :: 3.11',
    'Programming Language :: Python :: 3.12',
    'Programming Language :: Python :: 3.13',
    'Topic :: Scientific/Engineering',
]

REQUIRES = [
    "biopython>=1.82",
    "matplotlib>=3.7.0",
    "numpy>=1.24.0",
    "scipy>=1.11.3",
    "scikit-learn>=1.4.2",
    "umap-learn>=0.5.0",
    "tqdm>=4.65.0"
]

# Auto-generate console_scripts for every canonical command and alias.
# Each public name "foo" maps to handler "bar", producing
# "pk_foo = phykit.phykit:bar".
_console_scripts = ["phykit = phykit.phykit:main"]
for command, handler in sorted(PUBLIC_COMMAND_TO_HANDLER.items()):
    _console_scripts.append(f"pk_{command} = phykit.phykit:{handler}")

setup(
    name="phykit",
    description="",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Jacob L. Steenwyk",
    author_email="jlsteenwyk@gmail.com",
    url="https://github.com/jlsteenwyk/phykit",
    packages=find_packages(),
    python_requires=">=3.10",
    classifiers=CLASSIFIERS,
    entry_points={
        "console_scripts": _console_scripts,
    },
    version=__version__,
    include_package_data=True,
    install_requires=REQUIRES,
)

## push new version to pypi
# rm -rf dist
# python3 -m build
# python3 -m twine check dist/*
# python3 -m twine upload dist/* -r pypi
# then push to anaconda
#
