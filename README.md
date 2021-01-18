# Vivarium-Bioscrape

Vivarium-bioscrape is a python package connecting the multiscale simulation framework Vivarium to the chemical reaction network (CRN) simulator Bioscrape. In particular, this package wraps bioscrape as a vivarium Process and includes a number of Derivers in the form of OneWayMaps to help connect CRNs to eachother and other processes.

Visit [the Vivarium Core
documentation](https://vivarium-core.readthedocs.io/) to learn how to
use the core Vivarium engine to create computational biology models.

Visit [the Bioscrape Wiki](https://github.com/biocircuits/bioscrape/wiki) to learn about the bioscrape simulator and find bioscrape usage examples.

## Installation

The easiest way to install this package is via: 

    pip install vivarium-bioscrape 

(with Python>=3.7)

This command will automatically install the package along with vivarium-core and bioscrape and all othe dependencies. 

Please note that Bioscrape is a cython extension module and requires a C++ compiler to be set up on your computer for installation. Please visit the [Bioscrape Wiki](https://github.com/biocircuits/bioscrape/wiki/Installation) for more information.


You can then run the vivarium-bioscrape process with the module format:

```
$ python -m vivarium_bioscrape.processes.bioscrape
```