import numpy as np

from vivarium.core.process import Deriver
from vivarium.core.experiment import pp
from vivarium_bioscrape.processes.bioscrape import Bioscrape
from vivarium_bioscrape.library.schema import array_to


class OneToMany(Deriver):

    name = 'one_to_many'
    defaults = {
        'stoichiometry': np.array([[]]),
        'one_key': '',
        'many_keys': [],
    }

    def __init__(self, parameters=None):
        super(OneToMany, self).__init__(parameters)

    def ports_schema(self):
        return {
            'one_delta': {
                self.parameters['one_key']: {
                    '_default': 0.0,
                }
            },
            'one': {
                self.parameters['one_key']: {
                    '_default': 0.0,
                    '_updater': 'accumulate',
                }
            },
            'many_deltas': {
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

        # pull the relevant values into arrays
        one_delta = np.array(states['one_delta'][self.parameters['one_key']])
        one = np.array(states['one'][self.parameters['one_key']])
        many_deltas = np.array([
            value if key in self.parameters['many_keys'] else 0.0
            for key, value in states['many_deltas'].items()])
        many = np.array([
            value if key in self.parameters['many_keys'] else 0.0
            for key, value in states['many'].items()])

        # transform
        one_update = np.sum(many_deltas)
        if one.any():
            many_updates = one_delta * many / one
        else:
            many_updates = np.zeros_like(many)

        # return update
        update = {
            'one': {self.parameters['one_key']: one_update},
            'many': array_to(states['many'].keys(), many_updates)
        }

        return update



def test_one_to_many():
    # get parameters from model2
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

    # test state 1
    state_1 = {
        'one': {
            one_key: np.sum([value
                             for key, value in many_state.items()
                             if 'rna' in key])},
        'one_delta': {
            one_key: 0.0},
        'many': many_state,
        'many_deltas': {
            mol_id: 0.5
            for mol_id in many_state.keys()}}

    transform_1 = adaptor.next_update(0, state_1)
    print('state 1 transform:')
    pp(transform_1)

    # test state 2
    state_2 = {
        'one': {
            one_key: 10},
        'one_delta': {
            one_key: 0.5},
        'many': many_state,
        'many_deltas': {
            mol_id: 0.0
            for mol_id in many_state.keys()}}

    transform_2 = adaptor.next_update(0, state_2)
    print('state 2 transform:')
    pp(transform_2)


    import ipdb; ipdb.set_trace()



if __name__ == '__main__':
    test_one_to_many()

