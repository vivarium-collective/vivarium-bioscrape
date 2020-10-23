"""
===================
Bioscrape Connector
===================
"""
import os

import numpy as np

# vivarium core imports
from vivarium.core.experiment import Experiment
from vivarium.core.process import Generator
from vivarium.core.composition import (
    simulate_compartment_in_experiment,
    COMPARTMENT_OUT_DIR,
)
from vivarium.plots.simulation_output import plot_simulation_output
from vivarium.library.dict_utils import deep_merge

# import processes
from vivarium_bioscrape.composites.composite_general import insert_topology
from vivarium_bioscrape.processes.bioscrape import Bioscrape, get_model_species_ids
from vivarium_bioscrape.processes.ports_map import LinearMap

NAME = 'bioscrape_connector'


class BioscrapeConnector(Generator):
    name = NAME
    defaults = {
        'models': {},
        'connections': [],
        }

    def __init__(self, config=None, models=None, connections=None):
        if config is None:
            config = {}

        if 'models' not in config:
            config['models'] = models
        if 'connections' not in config:
            config['connections'] = connections

        super(BioscrapeConnector, self).__init__(config)
        self.topology = self.initial_topology(self.config)

        # TODO -- assert that connectors

    def initial_state(self, config={}):
        # TODO -- find and resolve conflicts
        network = self.generate()
        processes = network['processes']

        # make initial state by merging the initial states of individual processes
        initial_state = {}
        for name, process in processes.items():
            initial = process.initial_state()
            initial_with_store = {
                f'{name}_{store}': values
                for store, values in initial.items()}
            initial_state = deep_merge(initial_state, initial_with_store)

        return initial_state

    def generate_processes(self, config):

        models = {
            name: Bioscrape({'sbml_file': path})
            for name, path in config['models'].items()}

        connections = {}
        for (m1, m2, projection) in config['connections']:
            m1_species = get_model_species_ids(config['models'][m1])
            m2_species = get_model_species_ids(config['models'][m2])

            connector_config = {
                'row_keys': m1_species,
                'column_keys': m2_species,
                'projection': projection}

            connections[f'{m1}_{m2}_connector'] = LinearMap(connector_config)

        return {**models, **connections}

    def generate_topology(self, config):
        return self.topology

    def initial_topology(self, config):

        # connect the models with their stores
        models = {
            name: {
                'species': (f'{name}_species',),
                'delta_species': (f'{name}_deltas',),
                'rates': (f'{name}_rates',)}
            for name, path in config['models'].items()}

        # make connections between model stores
        connections = {}
        for (m1, m2, projection) in config['connections']:
            connections[f'{m1}_{m2}_connector'] = {
                'row': (f'{m1}_species',),
                'row_deltas': (f'{m1}_deltas',),
                'column': (f'{m2}_species',),
                'column_deltas': (f'{m2}_deltas',),
            }

        return {**models, **connections}

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



def main():
    '''Simulate the composite and plot results.'''

    # make an output directory to save plots
    out_dir = os.path.join(COMPARTMENT_OUT_DIR, NAME)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    model1_keys = get_model_species_ids('Notebooks/model1.xml')
    model3_keys = get_model_species_ids('Notebooks/model3.xml')

    # define the projection
    model1_vector = np.array(['rna' in key for key in model1_keys])
    model3_vector = np.array(['rna' in key for key in model3_keys])
    projection = np.outer(model1_vector/np.sum(model1_vector), model3_vector/np.sum(model3_vector))

    # configuration
    models = {
        '1': 'Notebooks/model1.xml',
        '3': 'Notebooks/model3.xml',
    }
    connections = [
        ('1', '3', projection),
    ]

    # make the composite
    composite = BioscrapeConnector(
        models=models,
        connections=connections,
    )

    # initial state
    initial_state = composite.initial_state()

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
