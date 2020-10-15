"""
run with:
> python -m vivarium_bioscrape.processes.linear_map
"""

import numpy as np

from vivarium.core.process import Deriver
from vivarium.core.experiment import pp
from vivarium_bioscrape.processes.bioscrape import Bioscrape
from vivarium_bioscrape.library.schema import array_to


class LinearMap(Deriver):
    name = 'linear_map'
    defaults = {
        'matrix': np.array([[]]),
        'row_keys': [],
        'column_keys': [],
        'projection_type': None,
    }

    def __init__(self, parameters=None):
        super(LinearMap, self).__init__(parameters)

    def ports_schema(self):
        return {
            'row_deltas': {
                '*': {
                    '_default': 0.0}},
            'row': {
                '*': {
                    '_default': 0.0,
                    '_updater': 'accumulate'}},
            'column_deltas': {
                '*': {
                    '_default': 0.0}},
            'column': {
                '*': {
                    '_default': 0.0,
                    '_updater': 'accumulate'}}}

    def projection(self, states):
        # return values / np.sum(values)
        return self.parameters['matrix']

    def next_update(self, timestep, states):

        # pull the relevant values into arrays
        row_deltas = np.array([
            value if key in self.parameters['row_keys'] else 0.0
            for key, value in states['row_deltas'].items()])
        row_values = np.array([
            value if key in self.parameters['row_keys'] else 0.0
            for key, value in states['row'].items()])
        column_deltas = np.array([
            value if key in self.parameters['column_keys'] else 0.0
            for key, value in states['column_deltas'].items()])
        column_values = np.array([
            value if key in self.parameters['column_keys'] else 0.0
            for key, value in states['column'].items()])

        # projections
        # NOTE: projects are missing normalization
        row_updates = np.dot(self.projection(states), column_deltas)
        column_updates = np.dot(row_deltas, self.projection(states))

        # return update
        update = {
            'column': array_to(states['column'].keys(), column_updates),
            'row': array_to(states['row'].keys(), row_updates),
        }

        return update


def CRNMap(LinearMap):
    def projection(self, states):
        return None


def test_linear_map():

    # get parameters from model2
    column_params = {
        'sbml_file': 'Notebooks/model2.xml'}
    column_process = Bioscrape(column_params)
    column_model = column_process.model
    columns = column_model.get_species_dictionary()
    column_keys = list(columns.keys())

    row_params = {
        'sbml_file': 'Notebooks/model3.xml'}
    row_process = Bioscrape(row_params)
    row_model = row_process.model
    rows = row_model.get_species_dictionary()
    row_keys = list(rows.keys())

    # make matrix
    column_vector = np.array(['rna' in key for key in columns.keys()])
    row_vector = np.array(['rna' in key for key in rows.keys()])
    matrix = np.outer(row_vector/np.sum(row_vector), column_vector/np.sum(column_vector))

    # make the process
    adaptor = LinearMap({
        'matrix': matrix,
        'column_keys': column_keys,
        'row_keys': row_keys,
    })

    # test state 1
    state_1 = {
        'row': rows,
        'row_deltas': {
            key: 0.5 for key in row_keys},
        'column': columns,
        'column_deltas': {
            key: 0.5 for key in column_keys},
    }

    transform_1 = adaptor.next_update(0, state_1)
    print('state 1 transform:')
    pp(transform_1)




    import ipdb;
    ipdb.set_trace()


if __name__ == '__main__':
    test_linear_map()

