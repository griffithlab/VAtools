# Contributing

We welcome improvements, extensions, bug reports and bug fixes to VAtools!  

## Bug Reports
To make bug reports most useful, please include:

- a brief summary of the problem - what did you expect, and what actually happened?
- precise steps to reproduce the problem, and minimal VCF and annotation files on which the steps can be run, if possible.

## Fixes or Feature Additions
Direct pull requests are just fine - fork the repository, create an new branch, then make a PR against this repo.

Before investing time and effort, it is often best to create an issue in 
the github repository reporting the bug and/or proposing the change, and we'll be happy to discuss the 
details and feasibility with you. 

## Running Tests
The test suite uses Python's built-in `unittest` framework and requires Python 3.10 or newer.

### 1. Get the code
Clone the repository (or your fork of it) and move into the top-level directory, the one containing `setup.py`:

```
git clone https://github.com/griffithlab/VAtools.git
cd VAtools
```

### 2. Create an isolated Python environment
This keeps VAtools' dependencies separate from the rest of your system. Use either conda:

```
conda create -n vatools python=3.12
conda activate vatools
```

or a virtual environment:

```
python3 -m venv .venv
source .venv/bin/activate
```

### 3. Install VAtools in editable mode
With the environment active, and from the top-level `VAtools` directory:

```
pip install -e .
```

The `-e` ("editable") flag installs VAtools as a link to your checkout instead of a copy, so changes you make to the code take effect immediately without reinstalling. This step also installs all dependencies, including `testfixtures`, which the tests need.

### 4. Run the tests
Still in the top-level directory, run the full suite:

```
python -m unittest discover -v
```

To run the tests for a single tool, or a single test, pass its module or test name:

```
python -m unittest tests.test_vep_annotation_reporter
python -m unittest tests.test_vep_annotation_reporter.VcfExpressionEncoderTests.test_all_vep_fields
```

Test inputs and expected outputs are in `tests/test_data/`. Most tests run a tool and compare its output to the expected results in that directory, so if you change a tool's output, you'll have to update the matching expected files too.

The same tests run automatically on GitHub Actions (`.github/workflows/tests.yml`) for every push and pull request to `master`, on Python 3.10 through 3.14. Please make sure they pass locally before opening a PR.
