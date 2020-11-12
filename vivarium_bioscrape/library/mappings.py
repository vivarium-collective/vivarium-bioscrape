
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
            #Normalize the projection matrix based upon the total concentrations of each target species
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


def stochiometric_map(parent_func, stochiometric_dictionary):
    #This converts any map function into a map that obeys stochiometric laws
    # parent_func is any parent mapping
    # stochiometric_dictionary is a dictionary: {species : {coupled species : stochiometric coefficient} }
    # Any species may have one or more coupled species (which always change along with it) 
    # using arbitrary stochiometric coefficients
    def map_function(states):
        deltas_dict = parent_func(states)
        new_deltas = dict(deltas_dict)
        for s in deltas_dict:
            if s in stochiometric_dictionary:
                for s2 in stochiometric_dictionary[s]:
                    coef = stochiometric_dictionary[s][s2]
                    if s2 in deltas_dict:
                        new_deltas[s2] += deltas_dict[s]*coef
                    else:
                        new_deltas[s2] = deltas_dict[s]*coef

        return new_deltas

    return map_function

def det_to_stoch_map(parent_func):
    #This converts any parent_func to produce integer changes in the species counts
    #using the volume of the underlying bioscrape processes. This allows for deterministic
    #and stochastic processes to be run together. Note: rounding errors may occur.

    def map_function(states):
        deltas_dict = parent_func(states)
        if 'global' in states and 'source_volume' in states['global']:
            V_source = states['global']['source_volume']
        else:
            V_source = 1.0

        for s in deltas_dict:
            sign = np.sign(deltas_dict[s])
            delta_int = np.floor(abs(deltas_dict[s]*V_source)) #Round down the change in number
            res = np.abs(deltas_dict[s]*V_source) - delta_int #Residual count

            #add/subtract an additional species with probability proportional to res.
            delta = sign*(delta_int + (np.random.random() < res))
            
            deltas_dict[s] = delta

        return deltas_dict

    return map_function

def stoch_to_det_map(parent_func):
    #This converts any parent_func that produces integer count changes into the species concentrations
    #using the volume of the underlying bioscrape processes. This allows for deterministic
    #and stochastic processes to be run together. 

    def map_function(states):
        deltas_dict = parent_func(states)
        if 'global' in states and 'target_volume' in states['global']:
            V_target = states['global']['target_volume']
        else:
            V_target = 1.0
        
        for s in deltas_dict:
            deltas_dict[s] = deltas_dict[s]/V_target #convert a change in count to a concentration

        return deltas_dict

    return map_function


