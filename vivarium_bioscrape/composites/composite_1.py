"""
======================
Bioscrape Composite #1
======================
"""
import os

# vivarium core imports
from vivarium.core.process import Generator
from vivarium.core.composition import (
    simulate_compartment_in_experiment,
    COMPARTMENT_OUT_DIR,
)
from vivarium.library.dict_utils import deep_merge
from vivarium.plots.simulation_output import plot_simulation_output

# import processes
from vivarium_bioscrape.processes.bioscrape import Bioscrape


NAME = 'composite_1'

class BioscrapeComposite(Generator):

    name = NAME
    defaults = {
        'model_1': {
            'sbml_file': 'Notebooks/model1.xml'},
        'model_3': {
            'sbml_file': 'Notebooks/model3.xml'
        },
    }

    def __init__(self, config):
        super(BioscrapeComposite, self).__init__(config)

    def initial_state(self, config={}):
        # TODO -- find and resolve conflicts

        # get the processes
        network = self.generate()
        processes = network['processes']

        # make initial state by merging the initial states of individual processes
        initial_state = {}
        for name, process in processes.items():
            initial_state = deep_merge(initial_state, process.initial_state())
        return initial_state

    def generate_processes(self, config):
        return {
            'model_1': Bioscrape(config['model_1']),
            'model_3': Bioscrape(config['model_3']),
        }

    def generate_topology(self, config):
        return {
            'model_1': {
                'species': ('species',),
                'delta_species': ('deltas',),
                'rates': ('rates',),
            },
            'model_3': {
                'species': ('species',),
                'delta_species': ('deltas',),
                'rates': ('rates',),
            },
        }

def main():
    '''Simulate the composite and plot results.'''

    # make an output directory to save plots
    out_dir = os.path.join(COMPARTMENT_OUT_DIR, NAME)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # make the composite
    composite = BioscrapeComposite({})

    initial_state = composite.initial_state()
    initial_state['species']['dna_G'] = 10.0

    # run a simulation
    sim_settings = {
        'total_time': 1000,
        'initial_state': initial_state}

    output = simulate_compartment_in_experiment(composite, sim_settings)

    # plot
    plot_settings = {}
    plot_simulation_output(output, plot_settings, out_dir)


if __name__ == '__main__':
    main()
