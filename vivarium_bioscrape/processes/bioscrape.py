'''
Execute by running: ``python -m  bioscrape.processes.bioscrape.py``

TODO: Replace the bioscrape code to implement your own process.
'''

import os
import numpy as np

from vivarium.core.process import Process
from vivarium.core.composition import (
    simulate_process,
    PROCESS_OUT_DIR,
)
from vivarium.plots.simulation_output import plot_simulation_output
from bioscrape.types import Model, Volume
from bioscrape.simulator import DeterministicSimulator, ModelCSimInterface, VolumeSSASimulator, SafeModelCSimInterface
from bioscrape.lineage import LineageCSimInterface, LineageModel, LineageSSASimulator, LineageVolumeCellState

NAME = 'bioscrape'

class Bioscrape(Process):
    '''
    This process provides a wrapper around a bioscrape model, interface, and simulator.
    It allows for stochastic or determinstic simulation, variable volume, and generates ports
    to access all bioscrape species and rate parameters.
    '''

    # give the process a name, so that it can register in the process_repository
    name = NAME

    # declare default parameters as class variables
    defaults = {
        'internal_dt': 0.01,
        'stochastic': False,
        'initial_volume': 1.0,
        'safe_mode':False,
        'lineage':False
    }

    def __init__(self, parameters=None):
        super(Bioscrape, self).__init__(parameters)

        # get the parameters out of self.parameters if available, or use defaults
        if 'sbml_file' not in self.parameters and 'bioscrape_model' not in self.parameters:
            raise ValueError("Bioscrape Process requires either an sbml_file or bioscrape_model parameter.")
        elif 'sbml_file' not in self.parameters and isinstance(self.parameters['bioscrape_model'], Model):
            # load the sbml file to create the model
            self.sbml_file = None
            self.model = self.parameters['bioscrape_model']
        elif isinstance(self.parameters['sbml_file'], str) and 'bioscrape_model' not in self.parameters:
            self.sbml_file = self.parameters['sbml_file']
            if self.parameters["lineage"]:
                self.model = LineageModel(sbml_filename = self.sbml_file,)
            else:
                self.model = Model(sbml_filename = self.sbml_file, sbml_warnings = False)
        elif isinstance(self.parameters['sbml_file'], str) and isinstance(self.parameters['bioscrape_model'], Model):
            raise ValueError("Bioscrape recieved an sbml_file and a bioscrape_model. Please use one or the other.")
        else:
            raise ValueError(f"Bioscrape did not recieve a valid bioscrape_model "
                             f"(recieved: {self.parameters['bioscrape_model']} or a "
                             f"valid sbml_file (recieved: {self.parameters['sbml_file']}).")
        
        self.internal_dt = self.parameters['internal_dt']
        self.stochastic = self.parameters['stochastic']
        
        
        #Toggle using Lineage Model
        if self.parameters["lineage"]:
            if not self.stochastic:
                raise ValueError("Bioscrape lineage only available with stochastic = True")
            self.simulator = LineageSSASimulator()
            self.interface = LineageCSimInterface(self.model)
            self.volume = LineageVolumeCellState() #Create an internal bioscrape Volume
        #Otherwise use normal bioscrape models
        else: 
            # create the interface
            if self.parameters["safe_mode"]:
                self.interface = SafeModelCSimInterface(self.model, max_species_count = 10**8)
            else:
                self.interface = ModelCSimInterface(self.model)

            #Stochastic
            if self.stochastic:
                self.simulator = VolumeSSASimulator()
            #Not Stochastic
            elif not self.stochastic:
                self.interface.py_prep_deterministic_simulation()
                # create a Simulator
                self.simulator = DeterministicSimulator()
            self.volume = Volume() #Create an internal bioscrape Volume

        #Set dt
        self.interface.py_set_dt(self.internal_dt)
        #Set the volume
        self.volume.py_set_volume(self.parameters['initial_volume']) 

    def get_species_names(self):
        #Gets the names of teh species in a bioscrape model
        model_species = self.model.get_species_dictionary()
        return list(model_species.keys())

    def get_state(self, array):
        #Gets the state of a bioscrape simulation
        mapping = self.model.get_species2index()

        return {
            species: array[index]
            for species, index in mapping.items()}

    def initial_state(self, config=None):
        #gets the current (or initial) state of a bioscrape simulation.
        if config is None:
            config = {}
        self.model.set_species(config)
        state = self.model.get_species_array()
        return {
            'species': self.get_state(state)}

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

        #Different divide settings between stochastic and determinsitic CRNs
        if self.stochastic:
            divider = "binomial"
        else:
            divider = "set" #division does not change concentrations

        return {
            'species': {
                species: {
                    '_default': 0.0,
                    '_updater': 'accumulate',
                    '_emit': True,
                    'divider':divider}
                for species in self.model.get_species()},
            'delta_species': {
                species: {
                    '_default': 0.0,
                    '_updater': 'set',
                    '_emit': False}
                for species in self.model.get_species()},
            'rates': {
                p: {
                    '_default': self.model.get_parameter_dictionary()[p],
                    '_updater':'set',
                }
                for p in self.model.get_param_list()},
            'globals': {
                'volume': {
                    '_default':self.parameters['initial_volume'],
                    '_updater': 'accumulate',
                    '_emit': True}
            }
        }

    def next_update(self, timestep, states):
        if 'species' in states:
            self.model.set_species(states['species'])
            self.interface.py_set_initial_state(self.model.get_species_array())
        if 'rates' in states:
            self.model.set_params(states['rates'])

        #Set Volume if needed
        if 'volume' in states:
            self.volume.py_set_volume(states['globals']['volume'])

        # create the interface

        timepoints = np.arange(0, timestep, self.internal_dt)
        if self.parameters["lineage"]:
            output = self.simulator.py_SimulateSingleCell(timepoints, Model = self.model, interface = self.interface, v = self.volume)
        elif self.stochastic:
            output = self.simulator.py_volume_simulate(self.interface, self.volume, timepoints)
        else:
            output = self.simulator.py_simulate(self.interface, timepoints)

        result = output.py_get_result()[-1]
        result_state = self.get_state(result)
        delta = get_delta(states['species'], result_state)
        rates = self.model.get_parameter_dictionary()

        #If the simulation is a volume simulation, return the change in volume
        if getattr(output, "py_get_volume", None) is not None:
            Vi = output.py_get_volume()[0]
            Vf = output.py_get_volume()[-1]
            deltaV = Vf-Vi
        else:
            deltaV = 0

        return {
            'species': delta,
            'delta_species': delta,
            'rates':rates,
            'globals': {'volume': deltaV}}

    def get_model(self):
        return self.model

    def get_model_species(self):
        return self.model.get_species_dictionary()

    def get_model_species_ids(self):
        return list(self.model.get_species_dictionary().keys())

    def get_volume(self):
        return self.volume.py_get_volume()


def get_delta(before, after):
    # assuming before and after have the same keys
    return {
        key: after[key] - before_value
        for key, before_value in before.items()}


def run_bioscrape_process():
    #Create a bioscrape process
    initial_parameters = {
        'sbml_file': 'Notebooks/model1.xml'}
    bioscrape_process = Bioscrape(initial_parameters)

    # run the simulation
    sim_settings = {
        'total_time': 10,
        'initial_state': bioscrape_process.initial_state()}
    output = simulate_process(bioscrape_process, sim_settings)
    

    # Return the data from the simulation.
    return output

def test_bioscrape_instantiation():

    #Deterministic Case
    bioscrape_process = Bioscrape({'sbml_file': 'Notebooks/model1.xml'})
    assert isinstance(bioscrape_process.model, Model)
    assert isinstance(bioscrape_process.interface, ModelCSimInterface)
    assert bioscrape_process.stochastic == False
    assert isinstance(bioscrape_process.simulator, DeterministicSimulator)

    #Stochastic Case
    bioscrape_process = Bioscrape({'sbml_file': 'Notebooks/model1.xml', "stochastic":True})
    assert isinstance(bioscrape_process.model, Model)
    assert isinstance(bioscrape_process.interface, ModelCSimInterface)
    assert bioscrape_process.stochastic == True
    assert isinstance(bioscrape_process.simulator, VolumeSSASimulator)

    #Custom Model Case
    M = Model(species = ["S"])
    bioscrape_process2 = Bioscrape({'bioscrape_model': M})
    assert M == bioscrape_process2.model
    assert bioscrape_process2.sbml_file is None

def test_next_update():
    initial_parameters = {
        'sbml_file': 'Notebooks/model1.xml'}
    bioscrape_process = Bioscrape(initial_parameters)

    initial_state = bioscrape_process.initial_state()
    output = bioscrape_process.next_update(1.0, initial_state)

    assert "species" in output
    assert "delta_species" in output
    assert all([output["species"][s] == output["delta_species"][s] for s in output["species"]])
    assert "rates" in output

    #set all the rates to 0
    state = {
        "rates": {p:0 for p in output["rates"]}, 
        "species":bioscrape_process.initial_state()["species"]
        }
    output2 = bioscrape_process.next_update(1.0, state)
    #nothing should change in the simulation
    assert all([output2["species"][s] == output2["delta_species"][s] ==0 for s in output2["species"]])

def test_bioscrape_process():
    '''Test that the process runs correctly.

    This will be executed by pytest.
    '''

    output = run_bioscrape_process()

    #DNA concentration should be constant
    assert all([v == output["species"]["dna_G"][0] for v in output["species"]["dna_G"]])
    #RNA concentration should be increasing
    assert all([output["species"]["rna_T"][i] < output["species"]["rna_T"][i+1] for i in range(len(output["species"]["rna_T"])-1)])
    #Protein concentration should be increasing
    assert all([output["species"]["protein_X"][i] < output["species"]["protein_X"][i+1] for i in range(len(output["species"]["protein_X"])-1)]) 

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

