
from vivarium.library.schema import array_from, array_to
import numpy as np


def one_to_one_map(bsp1, bsp2, species1, species2):
    #Takes two Bioscrape Processes, bsp1 and bsp2 and
    #two lists of species with species1 in bsp1 and species2 in bsp2.
    #Produces the mapping species1[i] --> species2[i]

    if len(species1) != len(species2):
        raise ValueError("species1 and species2 must be the same length.")

    all_species1 = bsp1.get_model_species_ids()
    all_species2 = bsp2.get_model_species_ids()

    projection = np.zeros((len(all_species1), len(all_species2)))

    #Create the projection matrix
    for i, s1 in enumerate(species1):
        if s1 not in all_species1:
            raise ValueError(f"{s1} not found in Bioscrape Process {bsp1.sbml_file}.")
        else:
            ind1 = all_species1.index(s1)

        s2 = species2[i]
        if s2 not in all_species2:
            raise ValueError(f"{s2} not found in Bioscrape Process {bsp2.sbml_file}.")
        else:
            ind2 = all_species2.index(s2)

        projection[ind1, ind2] = 1

    def map_function(states):
        source_delta_array = array_from(states['source_deltas'])
        output_delta_array = np.dot(source_delta_array, projection)
        return array_to(all_species2, output_delta_array)

    return map_function
