# Vivarium-Bioscrape

Visit [the Vivarium Core
documentation](https://vivarium-core.readthedocs.io/) to learn how to
use the core Vivarium engine to create computational biology models.

## Installation

To run this package from within a local vivarium-bioscrape directory, you need 
[Bioscrape](https://github.com/biocircuits/bioscrape) installed locally, and 
set the python path to that directory. This is required because Bioscrape is not 
available as a pip library.

We recommend setting up a [pyenv](https://github.com/pyenv/pyenv) in the bioscrape 
directory (with Python>=3.7), and then following the [Bioscrape installation instructions](
https://github.com/biocircuits/bioscrape/wiki/Installation).

Once bioscrape is set up, go to your vivarium-bioscrape directory and set your python path:

```
$ export PYTHONPATH="/path/to/bioscrape"
```

You can then run the vivarium-bioscrape process with the module format:

```
$ python -m vivarium_environment.processes.bioscrape
```