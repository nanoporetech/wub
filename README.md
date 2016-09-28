Tools and software library developed by the ONT Applications group
==================================================================


Installation
------------

Install the package:

```
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

Documentation
-------------

Documentation can be found at: https://apps.git.oxfordnanolabs.local/wub

Contributing
------------

- Please fork the repository and create a merge request to contribute.
- Avoid overlap of functionality with other ONT packages. First check if [tang](https://git.oxfordnanolabs.local/research/tang) or one of the [metrichor-bio](https://git.oxfordnanolabs.local/groups/metrichor-bio) modules has already implemented what you need.
- Please respect the structure outlined in `scripts/_template_script.py` from command line tools so documentation can be generated automatically.
- All functions should have at least one test (with the exception of plotting functions).
- Use [bumpversion](http://bit.ly/2cSUryt) to manage package versioning.
- The code should be [PEP8](https://www.python.org/dev/peps/pep-0008) compliant, which can be tested by `make lint`.
- For more guidance regarding coding style please refer to the [Tao of Tang](https://git.oxfordnanolabs.local/research/tang/blob/master/tao.md).

TODO
----


### What is a wub anyway?

It's a short, memorable name not on pypi yet. Also either [this](https://en.wikipedia.org/wiki/Beyond_Lies_the_Wub) or [this](http://www.urbandictionary.com/define.php?term=Wub).


