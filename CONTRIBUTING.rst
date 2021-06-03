Contributing
============

Welcome, friend!
Open-source software isn't open open-source without the community.
We appreciate your interest and welcome all contributions.
To help keep everything moving smoothly,
we have a few guidelines.

Bugs
----

If you think you've found a bug,
let us know `here`_.
We'll do our best to deal with it ASAP,
but please be patient as we also work many other projects!

.. _here: https://github.com/IMS-Bio2Core-Facility/BIC086-Sophie-Austin/issues

Developing
----------

If you think you can fix one of the bugs,
or would like to submit a new feature,
then let's get coding!

Once you've cloned the repository,
fork it,
and get your development environment set up.
We use `conda`_, `nox`_, and `pre-commit`_ to handle testing and linting.
Between them,
they make sure that all checks run in isolated environments.
Please make sure you activate them before making any changes,
using the commands below:

.. code-block:: shell

   conda env create -f environments/devel.yml
   conda activate bic086
   pre-commit install --install-hooks
   pre-commit install --hook-type pre-push

Now,
you don't even have to think about linting or testing.
When you commit a change,
pre-commit will automatically run `black`_,
`isort`_,
`mypy`_,
and a suite of `flake8`_-based linters.
When you push a change,
nox will run a more robust,
but slightly slower,
series of lints,
including security checks,
testing with `pytest`_,
and automatic doc building with `sphinx`_.

Speaking of documentation and testing -
if you add new code,
please add documentation and tests for it as well.
We use `napoleon numpy`_ doc strings.
Include the docstring with the code,
and sphinx will handle everything else!
Once your happy with your code,
open a pull-request,
and we will reveiw ASAP.

.. _conda: https://docs.conda.io/en/latest/
.. _nox: https://nox.thea.codes/en/stable/
.. _pre-commit: https://pre-commit.com/
.. _black: https://github.com/psf/black
.. _isort: https://pycqa.github.io/isort/
.. _mypy: https://github.com/cli/cli
.. _flake8: https://flake8.pycqa.org/en/latest/
.. _pytest: https://docs.pytest.org/en/6.2.x/
.. _sphinx: https://www.sphinx-doc.org/en/master/
.. _napoleon numpy: https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html

From the Command-Line
---------------------

If you haven't heard of it already,
give a peek to to `gh`_,
GitHub's official CLI.
It allows to manage all of the above steps from the command-line,
from forking,
to raising issues,
and checking on the status of your pull request.
Not a necessity,
but for you terminal warriors out there,
it just might help!

.. _gh: https://github.com/cli/cli
