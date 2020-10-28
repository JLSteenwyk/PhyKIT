<p align="center">
  <a href="https://github.com/jlsteenwyk/phykit">
    <img src="https://raw.githubusercontent.com/JLSteenwyk/PhyKIT/master/docs/_static/img/phykit_logo.png" alt="Logo" width="400">
  </a>
  <p align="center">
    <a href="https://jlsteenwyk.com/PhyKIT/">Docs</a>
    ·
    <a href="https://github.com/jlsteenwyk/phykit/issues">Report Bug</a>
    ·
    <a href="https://github.com/jlsteenwyk/phykit/issues">Request Feature</a>
  </p>
    <p align="center">
        <a href="https://lbesson.mit-license.org/" alt="License">
            <img src="https://img.shields.io/badge/License-MIT-blue.svg">
        </a>
        <a href="https://pypi.org/project/phykit/" alt="PyPI - Python Version">
            <img src="https://img.shields.io/pypi/pyversions/phykit">
        </a>
        <a href="https://github.com/JLSteenwyk/PhyKIT/actions" alt="Build">
            <img src="https://img.shields.io/github/workflow/status/jlsteenwyk/phykit/CI/master">
        </a>
        <a href="https://codecov.io/gh/jlsteenwyk/phykit" alt="Coverage">
          <img src="https://codecov.io/gh/jlsteenwyk/phykit/branch/master/graph/badge.svg?token=0J49I6441V">
        </a>
        <a href="https://github.com/jlsteenwyk/phykit/graphs/contributors" alt="Contributors">
            <img src="https://img.shields.io/github/contributors/jlsteenwyk/phykit">
        </a>
        <a href="https://twitter.com/intent/follow?screen_name=jlsteenwyk" alt="Author Twitter">
            <img src="https://img.shields.io/twitter/follow/jlsteenwyk?style=social&logo=twitter"
                alt="follow on Twitter">
        </a>
    </p>
</p>

PhyKIT is a UNIX shell toolkit for processing and analyzing phylogenomic data.<br /><br />
If you found PhyKIT useful, please cite *PhyKIT: a UNIX shell toolkit for processing and analyzing phylogenomic data*. bioRxiv. doi: [10.1101/2020.10.27.358143](https://www.biorxiv.org/content/10.1101/2020.10.27.358143v1).
<br /><br />

---

This documentation covers downloading and installing PhyKIT. Details about each function as well as tutorials for using PhyKIT are available in the <a href="https://jlsteenwyk.com/PhyKIT/">online documentation</a>.

<br />

**Installation** <br />
To install, use the following commands:
```shell
pip install phykit
```

<br />

To install from source, use the following commands:
```shell
git clone https://github.com/JLSteenwyk/PhyKIT.git
cd PhyKIT/
make install
```

If you run into permission errors when executing make install, create a virtual environment for your installation:
```shell
git clone https://github.com/JLSteenwyk/PhyKIT.git
cd PhyKIT/
python -m venv .venv
source .venv/bin/activate
make install
```
Note, the virtual environment must be activated to use phykit.

<br />
To test phykit installation, launch the help message

```shell
phykit -h
```


