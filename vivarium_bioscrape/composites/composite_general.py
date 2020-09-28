"""
======================
Bioscrape Composite #1
======================
"""
import os

# vivarium core imports
from vivarium.core.experiment import Experiment
from vivarium.core.process import Generator
from vivarium.core.composition import (
    simulate_compartment_in_experiment,
    plot_simulation_output,
    COMPARTMENT_OUT_DIR,
)
from vivarium.library.dict_utils import deep_merge

# import processes
from vivarium_bioscrape.processes.bioscrape import Bioscrape


NAME = 'bioscrape_composer'

class BioscrapeComposer(Generator):

    name = NAME
    defaults = {
        'model_paths': []}

    def __init__(self, config=None, model_paths=None):
        if config is None:
            config = {}

        if 'model_paths' not in config:
            config['model_paths'] = model_paths

        super(BioscrapeComposer, self).__init__(config)

        self.topology = self.initial_topology()

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
            # TODO: pull out name from path for key
            path: Bioscrape({'sbml_file': path})
            for path in self.config['model_paths']}

    def generate_topology(self, config):
        return self.topology

    def initial_topology(self):
        return {
            path: {
                'species': ('species',),
                'rates': ('rates',)}
            for path in self.config['model_paths']}

    def add_mappings(self, config=None, model=None, species=None, rates=None):
        if config is None:
            config = {}

        if model not in config:
            if species:
                config[model] = {
                    'species': species}
            if rates:
                config[model] = {
                    'rates': rates}
        elif model:
            raise ValueError('config and model keywords are conflicting for {}'.format(model))
            
        for model_path, port in config.items():
            for port_key, mappings in port.items():
                for mapping in mappings:
                    self.insert_topology(model_path, port_key, mapping)

    def insert_topology(self, model_path, port_key, mapping):
        self.topology[model_path][port_key] = insert_topology(
            self.topology[model_path][port_key], mapping)
        

def insert_topology(path, mapping):
    if isinstance(path, tuple):
        path = {
            '_path': path}

    from_name, to_name = mapping
    path[from_name] = (to_name,)

    return path


def main():
    '''Simulate the composite and plot results.'''

    # make an output directory to save plots
    out_dir = os.path.join(COMPARTMENT_OUT_DIR, NAME)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # make the composite
    composite = BioscrapeComposer(
        model_paths=[
            'Notebooks/model1b.xml',
            'Notebooks/model3.xml'])

    composite.add_mappings({
        'Notebooks/model1b.xml': {
            'species': [
                ('rna_R', 'rna_T')]}})

    composite.add_mappings(
        model='Notebooks/model1b.xml',
        species=[
            ('rna_R', 'rna_T')])

    initial_state = composite.initial_state()
    initial_state['species']['dna_G'] = 10.0

    # run a simulation
    sim_settings = {
        'total_time': 100,
        'initial_state': initial_state}

    output = simulate_compartment_in_experiment(composite, sim_settings)

    # plot
    plot_settings = {}
    plot_simulation_output(output, plot_settings, out_dir)


if __name__ == '__main__':
    main()
