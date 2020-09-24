"""
====================
Bioscrape Composites
====================
"""

# TODO: Delete this file before publishing your project.

# vivarium core imports
from vivarium.core.experiment import Experiment
from vivarium.core.process import Generator

# import processes
from vivarium_bioscrape.processes.bioscrape import Bioscrape


class BioscrapeComposite(Generator):

    defaults = {
        'glucose_phosphorylation': {},
        'injector': {},
    }

    def __init__(self, config):
        super(BioscrapeComposite, self).__init__(config)

    def generate_processes(self, config):
        bioscrape_1 = Bioscrape()

        return {
            'bioscrape_1': bioscrape_1,
            'bioscrape_2': glucose_phosphorylation,
        }

    def generate_topology(self, config):
        return {
            'injector': {
                'internal': ('cell', ),
            },
            'glucose_phosphorylation': {
                'cytoplasm': ('cell', ),
                'nucleoside_phosphates': ('cell', ),
                'global': ('global', ),
            },
        }
