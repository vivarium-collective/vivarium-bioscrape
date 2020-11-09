
from vivarium.library.schema import array_from, array_to

def one_to_one_map(species1, species2, model1, model2):

    def map_function(states):
        source_delta_array = array_from(states['source_deltas'])
