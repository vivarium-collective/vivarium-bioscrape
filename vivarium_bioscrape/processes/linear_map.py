"""
run with:
> python -m vivarium_bioscrape.processes.ports_map
"""

import argparse

import numpy as np

from vivarium.core.process import Deriver
from vivarium.core.experiment import pp
from vivarium.library.schema import array_from, array_to
from vivarium_bioscrape.processes.bioscrape import get_model_species


class LinearMap(Deriver):
    name = 'linear_map'
    defaults = {
        'projection': np.array([[]]),
        'row_keys': [],
        'column_keys': [],
    }

    def __init__(self, parameters=None):
        super(LinearMap, self).__init__(parameters)

    def initial_state(self, config=None):
        return {}

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

    def map(self, states):
        return self.parameters['projection']

    def next_update(self, timestep, states):

        # pull the relevant values into arrays
        row_deltas = np.array([
            value if key in self.parameters['row_keys'] else 0.0
            for key, value in states['row_deltas'].items()])
        column_deltas = np.array([
            value if key in self.parameters['column_keys'] else 0.0
            for key, value in states['column_deltas'].items()])

        # projection
        row_updates = np.dot(self.map(states), column_deltas)
        column_updates = np.dot(row_deltas, self.map(states))

        # return update
        update = {
            'column': array_to(states['column'].keys(), column_updates),
            'row': array_to(states['row'].keys(), row_updates)}

        return update


def CRNMap(LinearMap):
    name = 'crn_map'
    defaults = {
        'stoichiometry': np.array([[]]),
        'row_keys': [],
        'column_keys': [],
    }

    def __init__(self, parameters=None):
        super(CRNMap, self).__init__(parameters)

    def map(self, states):
        return self.parameters['stoichiometry']


# testing
def test_linear_map():
    rows = get_model_species('Notebooks/model3.xml')
    columns = get_model_species('Notebooks/model2.xml')
    row_keys = list(rows.keys())
    column_keys = list(columns.keys())

    # define the projection
    column_vector = np.array(['rna' in key for key in column_keys])
    row_vector = np.array(['rna' in key for key in row_keys])
    projection = np.outer(row_vector/np.sum(row_vector), column_vector/np.sum(column_vector))

    # make the process
    adaptor = LinearMap({
        'projection': projection,
        'column_keys': column_keys,
        'row_keys': row_keys})

    # test state 1
    state_1 = {
        'row': rows,
        'row_deltas': {
            key: 0.5 for key in row_keys},
        'column': columns,
        'column_deltas': {
            key: 0.5 for key in column_keys}}

    transform_1 = adaptor.next_update(0, state_1)
    print('state 1 transform:')
    pp(transform_1)



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='tumor cells')
    parser.add_argument('--linear', '-l', action='store_true', default=False)
    args = parser.parse_args()

    if args.linear:
        test_linear_map()
