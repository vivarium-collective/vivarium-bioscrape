'''
Execute by running: ``python -m  bioscrape.processes.bioscrape.py``

TODO: Replace the bioscrape code to implement your own process.
'''

from __future__ import absolute_import, division, print_function

import os
import numpy as np

from vivarium.core.process import Process
from vivarium.core.composition import (
    simulate_process_in_experiment,
    plot_simulation_output,
    PROCESS_OUT_DIR,
)

from bioscrape.types import Model
from bioscrape.simulator import DeterministicSimulator, ModelCSimInterface


NAME = 'bioscrape'



def get_delta(before, after):
    # assuming before and after have the same keys
    return {
        key: after[key] - before_value
        for key, before_value in before.items()}
    

class Bioscrape(Process):
    '''
    This mock process provides a basic bioscrape that can be used for a new process
    '''

    # give the process a name, so that it can register in the process_repository
    name = NAME

    # declare default parameters as class variables
    defaults = {
        'sbml_file': 'model.xml',
        'internal_dt': 0.01}

    def __init__(self, parameters=None):
        if parameters is None:
            parameters = {}

        super(Bioscrape, self).__init__(parameters)

        # get the parameters out of initial_parameters if available, or use defaults
        self.sbml_file = self.parameters['sbml_file']
        self.internal_dt = self.parameters['internal_dt']

        # load the sbml file to create the model
        self.model = Model(sbml_filename = self.sbml_file, sbml_warnings = False)

        # create the interface
        self.interface = ModelCSimInterface(self.model)

        # create a Simulator
        self.simulator = DeterministicSimulator()

    def get_state(self, array):
        mapping = self.model.get_species2index()

        return {
            species: array[index]
            for species, index in mapping.items()}

    def initial_state(self):
        state = self.model.get_species_array()
        return self.get_state(state)

    def ports_schema(self):
        '''
        ports_schema returns a dictionary that declares how each state will behave.
        Each key can be assigned settings for the schema_keys declared in Store:

        * `_default`
        * `_updater`
        * `_divider`
        * `_value`
        * `_properties`
        * `_emit`
        * `_serializer`
        '''

        return {
            'species': {
                species: {
                    '_default': 0.0,
                    '_updater': 'accumulate',
                    '_emit': True}
                for species in self.model.get_species()},

            'rates': {}}

    def next_update(self, timestep, states):
        self.interface.py_prep_deterministic_simulation()
        self.model.set_species(states['species'])

        timepoints = np.arange(0, timestep, self.internal_dt)
        output = self.simulator.py_simulate(self.interface, timepoints)
        result = output.py_get_result()[-1]
        result_state = self.get_state(result)
        delta = get_delta(states['species'], result_state)

        update = {
            'species': delta}
        
        return update


def run_bioscrape_process():
    '''Run a simulation of the process.

    Returns:
        The simulation output.
    '''
    # initialize the process by passing initial_parameters
    initial_parameters = {
        'sbml_file': 'Notebooks/model1.xml'}
    bioscrape_process = Bioscrape(initial_parameters)
    state = bioscrape_process.initial_state()

    # run the simulation
    sim_settings = {
        'total_time': 10,
        'initial_state': {
            'species': state}}
    output = simulate_process_in_experiment(bioscrape_process, sim_settings)

    # Return the data from the simulation.
    return output


def test_bioscrape_process():
    '''Test that the process runs correctly.

    This will be executed by pytest.
    '''
    output = run_bioscrape_process()
    # TODO: Add assert statements to ensure correct performance.


def main():
    '''Simulate the process and plot results.'''
    # make an output directory to save plots
    out_dir = os.path.join(PROCESS_OUT_DIR, NAME)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    output = run_bioscrape_process()

    # plot the simulation output
    plot_settings = {}
    plot_simulation_output(output, plot_settings, out_dir)


if __name__ == '__main__':
    main()
