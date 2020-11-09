"""
run with:
> python -m vivarium_bioscrape.processes.one_way_map
"""

import numpy as np

from vivarium.core.process import Deriver
from vivarium.core.experiment import pp
from vivarium.library.schema import array_from, array_to
from vivarium_bioscrape.processes.bioscrape import Bioscrape


class OneWayMap(Deriver):
    name = 'one_way_map'
    defaults = {
        'source': '',
        'target': '',
        'source_keys': [],
        'target_keys': [],
        'map_function': lambda states: {}}

    def __init__(self, parameters=None):
        super(OneWayMap, self).__init__(parameters)
        self.map_function = self.parameters['map_function']

    def initial_state(self):
        return {}

    def ports_schema(self):
        return {
            'source_deltas': {
                species_id: {
                    '_default': 0.0,
                } for species_id in self.parameters['source_keys']},
            'target_state': {
                species_id: {
                    '_default': 0.0,
                } for species_id in self.parameters['target_keys']}}

    def next_update(self, timestep, states):
        output = self.map_function(states)
        return {
            'target_state': output,
        }



def test_one_way_map():
    bioscrape_process_input = Bioscrape(parameters = {
        'sbml_file':'Notebooks/model3.xml'
        })
    bioscrape_process_output = Bioscrape(parameters = {
        'sbml_file':'Notebooks/model2.xml'
        })

    input_species = bioscrape_process_input.get_model_species() #get_model_species('Notebooks/model3.xml')
    output_species = bioscrape_process_output.get_model_species() #get_model_species('Notebooks/model2.xml')
    input_keys = list(input_species.keys())
    output_keys = list(output_species.keys())

    # define the projection
    input_vector = np.array(['rna' in key for key in input_keys])
    output_vector = np.array(['rna' in key for key in output_keys])
    projection = np.outer(input_vector/np.sum(input_vector), output_vector/np.sum(output_vector))

    # define map function
    def map_function(states):
        input_array = array_from(states['source_deltas'])
        output_array = np.dot(input_array, projection)
        return array_to(output_keys, output_array)

    one_way = OneWayMap({
        'map_function': map_function,
        'source_keys': input_keys,
        'target_keys': output_keys})

    # use input state
    state = {
        'source_deltas': {key: 100 for key in input_species.keys()},
        'target_state': {key: 0.0 for key in output_species.keys()}}

    transform = one_way.next_update(0, state)



def test_one_to_many():
    #One species is converted into many species, evenly.
    pass

def test_many_to_one():
    #Many species are converted into one species
    pass

def test_many_to_many():
    #Many species are converted into many species, evenly.
    pass

def test_one_to_many_stochiometrically_compatable():
    #One species is converted to many species in a way that is 
    #compatable with the stochiometric matrix of the manys species.
    pass

def test_output_proportional():
    #One species is converted into many species, proportionally to their amounts.
    pass

def test_positive_negative():
    #Negative deltas are treated mapped differently than positive deltas.
    pass

if __name__ == '__main__':
    test_one_way_map()
