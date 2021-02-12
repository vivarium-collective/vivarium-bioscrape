"""
===================
Bioscrape Connector
===================

run with:
> python -m vivarium_bioscrape.composites.bioscrape_connector
"""
import os

import numpy as np

# vivarium core imports
from vivarium.library.schema import array_from, array_to
from vivarium.core.experiment import pp
from vivarium.core.process import Composer
from vivarium.core.composition import (
    simulate_composer,
    COMPARTMENT_OUT_DIR,
)
from vivarium.plots.simulation_output import plot_simulation_output
from vivarium.library.dict_utils import deep_merge

# import processes
from vivarium_bioscrape.processes.bioscrape import Bioscrape
from vivarium_bioscrape.processes.one_way_map import OneWayMap


NAME = 'bioscrape_connector'


class BioscrapeConnector(Composer):
    name = NAME
    defaults = {
        'models': {},
        'connections': [],
        }

    def __init__(self, config=None, models=None, connections=None):
        if config is None:
            config = {}
        if models is None:
            models = {}
        if connections is None:
            connections = {}

        if 'models' not in config:
            config['models'] = models
        if 'connections' not in config:
            config['connections'] = connections

        self.models = {}
        self.connections = {}
        self.connection_counts = {} #Stores the number of connections between pairs of processes, for naming purposes
        super(BioscrapeConnector, self).__init__(config)
        self.topology = self.initial_topology(self.config)

    def initial_state(self, config=None):
        """
        TODO -- automated resolution of conflicts
        """
        if config is None:
            config = {}
        self.generate()

        # set model values to config
        for name, species_dict in config.items():
            self.models[name].initial_state(species_dict)

        initial_state = {}
        for name, node in self.connections.items():
            source = node.parameters['source']
            target = node.parameters['target']
            map_function = node.map_function
            source_state = self.models[source].initial_state()
            target_state = self.models[target].initial_state()

            # Assume that source is going from 0 to the source process' initial state
            modified_source_state = {
                'source_deltas': source_state['species'],
                'target_state': target_state['species']}
            target_deltas = map_function(modified_source_state)

            for species, value in target_deltas.items():
                target_value = target_state['species'][species]
                if value != 0.0 and target_value != value:
                    raise ValueError(
                        f"initial condition ['{target}_species']['{species}']={target_value} "
                        f"does not match mapping from model {source} to model {target}")

            # set initial state according source and target species
            initial_state[f'{source}_species'] = source_state['species']
            initial_state[f'{target}_species'] = target_state['species']

        #Reset internal variables because generate_processes will be called again
        return initial_state

    def generate_processes(self, config):
        # make bioscrape processes
        
        self.connection_counts = {} #models and connections may be overwritten, so reset the counts
        for name, parameters in config['models'].items():
            self.add_model(name = name, bioscrape_parameters = parameters)

        for connection in config['connections']:
            source = connection['source']
            target = connection['target']
            map_function = connection['map_function']
            self.add_connection(source, target, map_function)

        # combine model and connection processes

        return {**self.models, **self.connections}

    def generate_topology(self, config):
        return self.topology

    def initial_topology(self, config):
        # connect the models with their stores
        models = {
            name: {
                'species': (f'{name}_species',),
                'delta_species': (f'{name}_deltas',),
                'rates': (f'{name}_rates',),
                'globals': (f'{name}_globals',)}
            for name, path in config['models'].items()}

        # make connections between model stores
        connections = {}
        for connection in config['connections']:
            source = connection['source']
            target = connection['target']

            #Count is used to allow multiple connections between models if desired
            if (source, target) in self.connection_counts:
                self.connection_counts[(source, target)] += 1
                count = self.connection_counts[(source, target)]
            else:
                count = 1
                self.connection_counts[(source, target)] = count

            connections[f'{source}_{target}_connector_{count}'] = {
                'source_deltas': (f'{source}_deltas',),
                'target_state': (f'{target}_species',),
                'globals' : {
                    f'source_volume':(f'{source}', 'globals','volume',),
                    f'target_volume':(f'{target}', 'globals','volume',),

                }
            }

        return {**models, **connections}

    def add_model(self, name = None, bioscrape_process = None, bioscrape_parameters = None):
        if name is None:
            name = str(len(self.models))

        #This test doesn't allow generate to be called multiple times
        #if name in self.models:
        #    raise ValueError(f"A model named {name} already exists!")

        if bioscrape_process is not None and bioscrape_parameters is not None:
            raise ValueError("Recieved both a bioscrape_process and bioscrape_parameters! Please use one or the other.")
        elif bioscrape_process is not None:
            if isinstance(bioscrape_process, Bioscrape):
                self.models[name] =  bioscrape_process
            else:
                raise TypeError("bioscrape_process must be a Process of type Bioscrape.")
        elif bioscrape_parameters is not None:
            if isinstance(bioscrape_parameters, dict):
                self.models[name] = Bioscrape(bioscrape_parameters)
            else:
                raise TypeError("bioscrape_parameters must be a dictionary.")
        else:
            raise ValueError("Recieved neither bioscrape_process nor bioscrape_parameters keywords. Please use one or the other (not both).")

        #Add process to the topology
        self.topology[name] = {
                'species': (f'{name}_species',),
                'delta_species': (f'{name}_deltas',),
                'rates': (f'{name}_rates',),
                'globals':(f'volume',),}


    def add_connection(self, source, target, map_function, connector_config = None):
        if source not in self.models:
            raise KeyError(f"source {source} is not the name of a model.")
        if target not in self.models:
            raise KeyError(f"source {target} is not the name of a model.")
        
        source_species = self.models[source].get_species_names()
        target_species = self.models[target].get_species_names()

        if connector_config is None:
            connector_config = {
                'source': source,
                'target': target,
                'source_keys': source_species,
                'target_keys': target_species,
                'map_function': map_function}

        if (source, target) in self.connection_counts:
            count = self.connection_counts[(source, target)]
            count += 1
            self.connection_counts[(source, target)] = count
        else:
            count = 1
            self.connection_counts[(source, target)] = count

        self.connections[f'{source}_{target}_connector_{count}'] = OneWayMap(connector_config)

        #Add process to the topology
        self.topology[f'{source}_{target}_connector_{count}'] = {
            'source_deltas': (f'{source}_deltas',),
            'target_state': (f'{target}_species',),
            'globals': {
                'source_volume': (f'{source}_globals', 'volume'),
                'target_volume': (f'{target}_globals', 'volume'),
            }
        }

def test_connector():
    
    bioscrape_process_1 = Bioscrape(parameters = {
        'sbml_file': 'Notebooks/model1.xml'
        })
    bioscrape_process_3 = Bioscrape(parameters = {
        'sbml_file': 'Notebooks/model3.xml'
        })

    model1_keys = bioscrape_process_1.get_model_species_ids()
    model3_keys = bioscrape_process_3.get_model_species_ids()

    # define the projection
    model1_vector = np.array(['rna_T' == key for key in model1_keys])
    model3_vector = np.array(['rna_T' == key for key in model3_keys])
    projection_1_3 = np.outer(
        model1_vector/np.sum(model1_vector),
        model3_vector/np.sum(model3_vector))
    projection_3_1 = np.outer(
        model3_vector/np.sum(model3_vector),
        model1_vector/np.sum(model1_vector))

        # define map function
    def map_1_3(states):
        input_array = array_from(states['source_deltas'])
        output_array = np.dot(input_array, projection_1_3)
        return array_to(model3_keys, output_array)

    def map_3_1(states):
        input_array = array_from(states['source_deltas'])
        output_array = np.dot(input_array, projection_3_1)
        return array_to(model1_keys, output_array)


    # configuration
    time_step = 1
    models = {
        '1': {
            'time_step': time_step,
            'sbml_file': 'Notebooks/model1.xml'},
        '3': {
            'time_step': time_step,
            'sbml_file': 'Notebooks/model3.xml'},
    }
    connections = [
        {'source': '1', 'target': '3', 'map_function': map_1_3},
        {'source': '3', 'target': '1', 'map_function': map_3_1}
    ]
    # make the composite
    composite = BioscrapeConnector(
        models=models,
        connections=connections,
    )
    

    ## Run a simulation
    # initial state
    config = {
        '1': {
            'rna_T': 10.0
        }
    }
    initial_state = composite.initial_state(config)

    #This occurs after composite.generate() is called
    assert len(composite.models) == 2
    assert len(composite.connections) == 2

    # run a simulation
    sim_settings = {
        'total_time': 10,
        'initial_state': initial_state}
    output = simulate_composer(composite, sim_settings)

    #DNA should be constant
    assert all([output['1_species']['dna_G'][0] == g for g in output['1_species']['dna_G']])

    #RNA should be the same between the two simulations
    assert output['1_species']['rna_T'] == output['3_species']['rna_T']

    #total RNAase should be constant
    RNAase_tot = output['3_species']['protein_RNAase'][0]

    rnas = [output['3_species']['protein_RNAase'][i]+output['3_species']['complex_protein_RNAase_rna_T_'][i] for i in range(len(output['3_species']['complex_protein_RNAase_rna_T_']))]

    #Slight numerical errors are possible, but total RNAase should be nearly constant
    assert  all([np.abs(s-RNAase_tot) < 10**-6 for s in rnas])

def main():
    '''Simulate the composite and plot results.'''
    #make an output directory to save plots
    out_dir = os.path.join(COMPARTMENT_OUT_DIR, NAME)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    output = test_instantiation()

    # plot simulation output
    plot_settings = {}
    plot_simulation_output(output, plot_settings, out_dir)


if __name__ == '__main__':
    main()
