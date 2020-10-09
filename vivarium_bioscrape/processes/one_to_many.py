import numpy as np

from vivarium_bioscrape.processes.bioscrape import Bioscrape
from vivarium.core.process import Deriver
from vivarium_bioscrape.library.schema import array_to
from vivarium.core.composition import (
    simulate_process_in_experiment,
    PROCESS_OUT_DIR,
)


class OneToMany(Deriver):
    name = 'one_to_many'

    # declare default parameters as class variables
    defaults = {
        'stoichiometry': np.array([[]]),
        'many_keys': [],
        'one_key': '',
    }

    def __init__(self, parameters=None):
        super(OneToMany, self).__init__(parameters)

    def ports_schema(self):
        return {
            'delta_one': {
                self.parameters['one_key']: {
                    '_default': 0.0,}
            },
            'one': {
                self.parameters['one_key']: {
                    '_default': 0.0,
                    '_updater': 'accumulate',
                }
            },
            'delta_many': {
                '*': {
                    '_default': 0.0,
                }
            },
            'many': {
                '*': {
                    '_default': 0.0,
                    '_updater': 'accumulate',
                }
            }
        }

    def next_update(self, timestep, states):
        delta_many = np.array(list(states['delta_many'].values()))
        delta_one = np.array(states['delta_one'][self.parameters['one_key']])
        many = np.array(list(states['many'].values()))

        # calculations
        new_one = float(np.sum(delta_many))
        new_many = delta_one / np.sum(many) * np.dot(delta_many, self.parameters['stoichiometry'])

        update = {
            'one': new_one,
            'many': array_to(states['many'].keys(), new_many)
        }

        import ipdb;
        ipdb.set_trace()

        return update



def test_adaptor():

    # parameters
    initial_parameters = {
        'sbml_file': 'Notebooks/model2.xml'}
    bioscrape_process = Bioscrape(initial_parameters)
    model = bioscrape_process.model
    stoichiometry = model.py_get_update_array()
    many_state = model.get_species_dictionary()
    rna_keys = [key for key in many_state.keys() if 'rna' in key]
    one_key = 'rna'

    # make the process
    adaptor = OneToMany({
        'stoichiometry': stoichiometry,
        'many_keys': rna_keys,
        'one_key': one_key,
        })

    # initial state
    state = {
        'one': {one_key: np.sum([value for key, value in many_state.items() if 'rna' in key])},
        'delta_one': {one_key: 0.0},
        'many': many_state,
        'delta_many': {
            mol_id: 0.0 for mol_id in many_state.keys()}
    }

    settings = {
        'total_time': 10,
        'initial_state': state}
    output = simulate_process_in_experiment(adaptor,settings)

    import ipdb; ipdb.set_trace()




if __name__ == '__main__':
    test_adaptor()

