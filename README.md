![ONT_logo](/ONT_logo.png)
-----------------------------

Tools and software library developed by the ONT Applications group
==================================================================

[![CircleCI](https://circleci.com/gh/nanoporetech/wub.svg?style=svg)](https://circleci.com/gh/nanoporetech/wub) [![Documentation Status](https://readthedocs.org/projects/wub/badge/?version=latest)](http://wub.readthedocs.io/en/latest/?badge=latest) [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/wub/README.html)

1\. Introduction:
-----------------

Various tools for basespace analysis of long reads.

2\. Getting Started:
-------------------

## Installation:

Set up a virtual environment 

```
virtualenv  --system-site-packages wub_env
source wub_env/bin/activate
pip install --upgrade pip
pip install requests[security]
```

Then install the package via pip:

```
pip install git+https://github.com/nanoporetech/wub.git
```

If you installed the package in a virtual environment then do not forget to
load it before using the package:

```
source wub_env/bin/activate
```

Run the following to leave the virtual environment:

```
deactivate
```

You can also clone the repository and install using `setup.py`:

```
git clone https://github.com/nanoporetech/wub.git
cd wub
python setup.py install
```

Install the package in developer mode:

```
python setup.py develop
```

Run the tests:

```
make test
```

Build the documentation:

```
make docs
```

Issue `make help` to get a list of `make` targets.

3\. Documentation
-----------------

Online documentation is avalaible at [wub.readthedocs.io](http://wub.readthedocs.io/en/latest/).

4\. Contributing
----------------

- Please fork the repository and create a merge request to contribute.
- Please respect the structure outlined in `scripts/_template_script.py` from command line tools so documentation can be generated automatically.
- All non-trivial functions should have at least one test (with the exception of plotting functions).
- Use your best judgement when deciding whether to put a piece of code in a script or make it more reusable by incorporating it into the module.
- Use [bumpversion](http://bit.ly/2cSUryt) to manage package versioning.
- The code should be [PEP8](https://www.python.org/dev/peps/pep-0008) compliant, which can be tested by `make lint`.

5\. Help:
---------

## Licence and Copyright:

(c) 2016 Oxford Nanopore Technologies Ltd.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

