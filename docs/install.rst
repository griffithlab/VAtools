Install
=======

The VAtools suite is written for Linux and Mac OS X.
If you are using Windows you will need to set up a
Linux environment, for example by using WSL or setting up a virtual machine.

VAtools requires Python 3.10 or above. Before running any
installation steps, check the Python version installed on your system:

.. code-block:: none

   python --version

pip
---

Install VAtools using ``pip``:

.. code-block:: none

   pip install vatools

You can verify the installation with:

.. code-block:: none

   pip show vatools

To upgrade an existing installation:

.. code-block:: none

   pip install vatools --upgrade

Docker
------

A Docker container for VAtools is available on DockerHub using the
`griffithlab/vatools <https://hub.docker.com/r/griffithlab/vatools/>`_ repo.

.. code-block:: none

   docker pull griffithlab/vatools

Run any tool inside the container by passing it as the command. For example:

.. code-block:: none

   docker run griffithlab/vatools vcf-readcount-annotator --help

To annotate a VCF with files on your local filesystem, mount the directory
containing your data as a volume:

.. code-block:: none

   docker run -v /path/to/data:/data griffithlab/vatools \
     vcf-readcount-annotator /data/input.vcf /data/readcounts.tsv DNA \
     -o /data/output.vcf
