Install
=======

The VAtools suite is written for Linux and Mac OS X.
If you are using Windows you will need to set up a
Linux environment, for example by setting up a virtual machine.

VAtools requires Python 3.5. Before running any
installation steps, check the Python version installed on your system:

.. code-block:: none

   python -V

If you don't have Python 3.5 installed, we recommend using `Conda
<http://conda.pydata.org/docs/py2or3.html>`_ to emulate a Python 3.5.
environment. We've encountered problems with users that already have Python
2.x installed when they also try to install Python 3.5. The defaults will
not be set correctly in that case. If you already have Python 2.x installed
we **strongly** recommmend using Conda instead of installing Python 3.5
locally.

Once you have set up your Python 3.5 environment correctly you can use
``pip`` to install VAtools. Make sure you have ``pip``
installed. ``pip`` is generally included in python distributions, but may
need to be upgraded before use. See the `instructions
<https://packaging.python.org/en/latest/installing/#install-pip-setuptools-and-wheel>`_
for installing or upgrading ``pip``.

After you have pip installed, type the following command on your Terminal:

.. code-block:: none

   pip install vatools

You can check that the ``vatools`` package has been installed
under the default environment by running this command:

.. code-block:: none

   pip show vatools

``pip`` will fetch and install VAtools and its dependencies for you.
After installing, each tool of the VAtools package is available in
its own command line tree directly from the Terminal.

If you have an old version of the vatools package installed you might
want to consider upgrading to the latest version:

.. code-block:: none

   pip install vatools --upgrade

Docker
------

A Docker container for VAtools is available on DockerHub using the
`griffithlab/vatools <https://hub.docker.com/r/griffithlab/vatools/>`_ repo.
