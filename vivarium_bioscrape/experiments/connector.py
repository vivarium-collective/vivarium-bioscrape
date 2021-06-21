
from vivarium.core.composition import simulate_composite

from vivarium_bioscrape.processes.bioscrape import Bioscrape
from vivarium_bioscrape.processes.one_way_map import one_to_one_map
from vivarium_bioscrape.composites.bioscrape_connector import BioscrapeConnector
import pytest


def test_one_to_one_map():
    # Create Bioscrape Processes
    # Model 1: is a simple transcription translation model
    bsp1 = Bioscrape(parameters={
        'sbml_file': 'Notebooks/model1.xml'
    })

    # Model 4 adds degredation by dilution to proteins and RNAs
    bsp2 = Bioscrape(parameters={
        'sbml_file': 'Notebooks/model4.xml'
    })

    # Create the map_function using one_to_one_map
    # rna_T --> rna_RNA
    # proteinX --> protein_Protein
    map_func_12 = one_to_one_map(bsp1, bsp2, {'rna_T': 'rna_RNA', 'protein_X': 'protein_Protein'})

    # Create the reverse map_function (which is also a one_to_one_map)
    map_func_21 = one_to_one_map(bsp2, bsp1, {'rna_RNA': 'rna_T', 'protein_Protein': 'protein_X'})

    # create a BioscrapeConnector called composite
    composer = BioscrapeConnector()

    # Add models to the composite
    composer.add_model(name="txtl", bioscrape_process=bsp1)
    composer.add_model(name="dilution", bioscrape_process=bsp2)

    # Add connections to the composite
    composer.add_connection(source="txtl", target="dilution", map_function=map_func_12)
    composer.add_connection(source="dilution", target="txtl", map_function=map_func_21)

    # Get the initial state

    with pytest.raises(ValueError):
        initial_state = composer.initial_state()

    #reset the initial state
    dilution_initial_state = {'dilution':{'rna_RNA':0, 'protein_Protein':0}}
    initial_state = composer.initial_state(dilution_initial_state)
    # make the composite
    composite = composer.generate()

    # Run a simulation
    sim_settings = {
        'total_time': 500,
        'initial_state': initial_state}
    output = simulate_composite(composite, sim_settings)




if __name__ == '__main__':
    test_one_to_one_map()
