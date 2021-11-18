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
        <br />
        <a href="https://lbesson.mit-license.org/" alt="License">
            <img src="https://img.shields.io/badge/License-MIT-blue.svg">
        </a>
        <a href="https://pypi.org/project/phykit/" alt="PyPI - Python Version">
            <img src="https://img.shields.io/pypi/pyversions/phykit">
        </a>
        <a href="https://academic.oup.com/bioinformatics/article-abstract/37/16/2325/6131675?redirectedFrom=fulltext">
          <img src="https://zenodo.org/badge/DOI/10.1093/bioinformatics/btab096.svg">
        </a>
    </p>
</p>

PhyKIT is a UNIX shell toolkit for processing and analyzing phylogenomic data.<br /><br />
If you found PhyKIT useful, please cite *PhyKIT: a broadly applicable UNIX shell toolkit for processing and analyzing phylogenomic data*. Bioinformatics. doi: [10.1093/bioinformatics/btab096](https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/btab096/6131675).
<br /><br />

---

This documentation covers downloading and installing PhyKIT. Details about each function as well as tutorials for using PhyKIT are available in the <a href="https://jlsteenwyk.com/PhyKIT/">online documentation</a>.

<br />

**Installation** <br />

**If you are having trouble installing PhyKIT, please contact the lead developer, Jacob L. Steenwyk, via [email](https://jlsteenwyk.com/contact.html) or [twitter](https://twitter.com/jlsteenwyk) to get help.**

To install using *pip*, we strongly recommend building a virtual environment to avoid software dependency issues. To do so, execute the following commands:
```shell
# create virtual environment
python -m venv .venv
# activate virtual environment
source .venv/bin/activate
# install phykit
pip install phykit
```

**Note, the virtual environment must be activated to use phykit.**

After using PhyKIT, you may wish to deactivate your virtual environment and can do so using the following command:
```shell
# deactivate virtual environment
deactivate
```

<br />

Similarly, to install from source, we strongly recommend using a virtual environment. To do so, use the following commands:
```shell
# download
git clone https://github.com/JLSteenwyk/PhyKIT.git
cd PhyKIT/
# create virtual environment
python -m venv .venv
# activate virtual environment
source .venv/bin/activate
# install
make install
```
To deactivate your virtual environment, use the following command:
```shell
# deactivate virtual environment
deactivate
```
**Note, the virtual environment must be activated to use phykit.**

<br />

To install via anaconda, execute the following command:
```shell
conda install -c jlsteenwyk phykit
```
Visit here for more information:
https://anaconda.org/JLSteenwyk/phykit

<br />

To test phykit installation, launch the help message

```shell
phykit -h
```
