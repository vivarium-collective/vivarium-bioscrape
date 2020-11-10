
from vivarium.library.schema import array_from, array_to
import numpy as np

def one_to_one_map(bsp1, bsp2, species_map):
    #Takes two Bioscrape Processes, bsp1 and bsp2 and
    # species_map: dictionar(a species in bsp1: a species in bsp2)
    #Produces the mapping species1[i] --> species2[i]

    all_species1 = bsp1.get_model_species_ids()
    all_species2 = bsp2.get_model_species_ids()

    projection = np.zeros((len(all_species1), len(all_species2)))

    #Create the projection matrix
    for s1 in species_map:
        if s1 not in all_species1:
            raise ValueError(f"{s1} not found in Bioscrape Process {bsp1.sbml_file}.")
        else:
            ind1 = all_species1.index(s1)

        s2 = species_map[s1]
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


def one_to_many_map(bsp1, bsp2, species_map, proportional = True):
    #Takes two Bioscrape Processes, bsp1 and bsp2 and
    #species_map { species in bsp1: [list of species in bsp2] }
    #Produces the mapping species1[i] --> [species2[i]]
    #proprtional: True/False
    #    if True the amount of species1[i] mapped to species2[i][j] is proportianal to the amount of species2[i][j]
    #        one_to_many_map(species1[i]) ~ species2[i][j]/sum(species[2][i])
    #    if False: the amount of species1[i] mapped to species2[i][j] ~ 1/len(species[2][i])

    all_species1 = bsp1.get_model_species_ids()
    all_species2 = bsp2.get_model_species_ids()

    projection = np.zeros((len(all_species1), len(all_species2)))

    #Create the projection matrix
    for s1 in species_map:
        if s1 not in all_species1:
            raise ValueError(f"{s1} not found in Bioscrape Process {bsp1.sbml_file}.")
        else:
            ind1 = all_species1.index(s1)

        for s2 in species_map[s1]:
            if s2 not in all_species2:
                raise ValueError(f"{s2} not found in Bioscrape Process {bsp2.sbml_file}.")
            else:
                ind2 = all_species2.index(s2)

            projection[ind1, ind2] = 1


    def map_function(states):
        source_delta_array = array_from(states['source_deltas'])

        if not proportional:
            normalizer = 1/np.sum(projection, 1)
            normalized_proj = normalizer[:, np.newaxis]*projection
        else:
            
            output_array = array_from(states['target_state'])
            normalized_proj = projection*output_array
            normalizer = np.sum(normalized_proj, 1)[:, np.newaxis]*projection
            normalized_proj = normalized_proj/(normalizer+(normalizer == 0))

            #Cycle through projections, find degenerate cases
            #Degenerate cases occur when none of the many species are present
            #In this case, the normalizer is uniform.
            for i in range(projection.shape[0]):
                if np.sum(normalizer[i, :]) == 0 and np.sum(projection[i, :]) > 0:
                    normalizer[i, :] = sum(projection[i, :])
                    normalized_proj[i, :] = projection[i, :]/normalizer[i, :]

        print(normalized_proj)
        output_delta_array = np.dot(source_delta_array, normalized_proj)
        return array_to(all_species2, output_delta_array)

    return map_function


def many_to_one_map(bsp1, bsp2, species_map):
    #Takes two Bioscrape Processes, bsp1 and bsp2 and
    #species_map = { species in bsp2 : [list of many species in bsp1] }
    #Produces the mapping sum(species1[i]) --> species2[i]

    all_species1 = bsp1.get_model_species_ids()
    all_species2 = bsp2.get_model_species_ids()

    projection = np.zeros((len(all_species1), len(all_species2)))

    #Create the projection matrix
    for s2 in species_map:

        if s2 not in all_species2:
            raise ValueError(f"{s2} not found in Bioscrape Process {bsp2.sbml_file}.")
        else:
            ind2 = all_species2.index(s2)

        for s1 in species_map[s2]:
            if s1 not in all_species1:
                raise ValueError(f"{s1} not found in Bioscrape Process {bsp1.sbml_file}.")
            else:
                ind1 = all_species1.index(s1)

            projection[ind1, ind2] = 1

    def map_function(states):
        source_delta_array = array_from(states['source_deltas'])
        output_delta_array = np.dot(source_delta_array, projection)
        return array_to(all_species2, output_delta_array)

    return map_function

