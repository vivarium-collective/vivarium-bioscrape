"""
run with:
> python -m vivarium_bioscrape.processes.one_way_map
"""

import numpy as np

from vivarium.core.process import Deriver
from vivarium.core.experiment import pp
from vivarium.library.schema import array_from, array_to
from vivarium_bioscrape.processes.bioscrape import Bioscrape
from vivarium_bioscrape.library.mappings import *


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


def test_one_to_one():
    #Each species is converted into a single other species
    bsp1 = Bioscrape(parameters = {
        'sbml_file':'Notebooks/model1.xml'
        })
    bsp2 = Bioscrape(parameters = {
        'sbml_file':'Notebooks/model4.xml'
        })

    map_func = one_to_one_map(bsp1, bsp2, {'rna_T':'rna_RNA', 'protein_X':'protein_Protein'})

    f = lambda s: 0 if s != "rna_T" else 1
    r = map_func({'source_deltas':{str(s):f(s) for s in bsp1.get_model_species_ids()}})
    assert r['rna_RNA'] == 1.0
    assert r['protein_Protein'] == 0.0

    r = map_func({'source_deltas':{str(s):f(s)==0 for s in bsp1.get_model_species_ids()}})
    assert r['rna_RNA'] == 0.0
    assert r['protein_Protein'] == 1.0

    r = map_func({'source_deltas':{str(s):1 for s in bsp1.get_model_species_ids()}})
    assert r['rna_RNA'] == 1.0
    assert r['protein_Protein'] == 1.0




def test_one_to_many():
    #One species is converted into many species, evenly.
    #Create Bioscrape Processes
    #Model 1: is a simple transcription translation model
    bsp1 = Bioscrape(parameters = {
        'sbml_file':'Notebooks/model1b.xml'
        })

    #Model 4 adds degredation by dilution to proteins and RNAs
    bsp2 = Bioscrape(parameters = {
        'sbml_file':'Notebooks/model4.xml'
        })

    #rna_RNA --> rna_T1 + rna_T2
    #protein_Protein --> protein_X1 + protein_X2
    #These species are added evenly
    map_func_21_not_prop = one_to_many_map(bsp2, bsp1,  {"rna_RNA":["rna_T1", "rna_T2"], 
                               "protein_Protein":['protein_X']}, proportional = False)


    f = lambda s: int(s[-1]) if "rna_T" in s else .5
    states =  {
    'source_deltas':{s:10 for s in bsp2.get_model_species_ids()}, 
    'target_state':{s:f(s) for s in bsp1.get_model_species_ids()}}

    r_not_prop = map_func_21_not_prop(states)

    assert r_not_prop['rna_T1'] == 5.
    assert r_not_prop['rna_T2'] == 5.
    assert r_not_prop['protein_X'] == 10
    assert r_not_prop['dna_G1'] == 0

    #These species are added proprotionally to their target concentration
    map_func_21_prop = one_to_many_map(bsp2, bsp1,  {"rna_RNA":["rna_T1", "rna_T2"], 
                               "protein_Protein":['protein_X']}, proportional = True)

    #Degenerate Case: Target concentration is 0
    states0 =  {
    'source_deltas':{s:10 for s in bsp2.get_model_species_ids()}, 
    'target_state':{s:0 for s in bsp1.get_model_species_ids()}}
    r_0 = map_func_21_prop(states0)
    #Should result in the same as the evenly distributed case above
    assert r_0 == r_not_prop

    #Non-generate case
    r_prop = map_func_21_prop(states)
    assert np.allclose(r_prop["rna_T1"], 10./3.)
    assert np.allclose(r_prop["rna_T2"], 2*10./3.)
    assert np.allclose(r_prop["protein_X"], 10.)

def test_many_to_one():
    #Many species are converted into one species
    
    bsp1 = Bioscrape(parameters = {
        'sbml_file':'Notebooks/model1b.xml'
        })

    #Model 4 adds degredation by dilution to proteins and RNAs
    bsp2 = Bioscrape(parameters = {
        'sbml_file':'Notebooks/model4.xml'
        })

    #Create the map_function using one_to_one_map
    #rna_T1 + rna_T2 --> rna_RNA
    #protein_X1 + protein_X2 --> protein_Protein
    #the many_to_one_map takes the two processes and a dictionary {one_species:[list of many species]}
    map_func_12 = many_to_one_map(bsp1, bsp2, {"rna_RNA":["rna_T1", "rna_T2"], 
                                   "protein_Protein":['protein_X']})

    f = lambda s: int(s[-1]) if "rna_T" in s else .5
    states =  {
    'source_deltas':{s:f(s) for s in bsp1.get_model_species_ids()}, 
    'target_state':{s:0 for s in bsp2.get_model_species_ids()}}

    r = map_func_12(states)

    assert r['rna_RNA'] == 3
    assert r['protein_Protein'] == .5

def test_one_to_many_stochiometrically_compatable():
    #One species is converted to many species in a way that is 
    #compatable with the stochiometric matrix of the manys species.
    pass


if __name__ == '__main__':
    test_one_way_map()
    test_one_to_one()
    test_one_to_many()
    test_many_to_one()
