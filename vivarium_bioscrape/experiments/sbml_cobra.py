import os

import numpy as np
import pylab as plt


# vivarium imports
from vivarium.core.control import Control
from vivarium.core.composition import (
    simulate_compartment_in_experiment,
    compartment_in_experiment,
)
from vivarium.core.composition import simulate_process_in_experiment
from vivarium.core.process import Process, Deriver, Generator
from vivarium.library.units import units

# imported processes
from vivarium_cobra.processes.dynamic_fba import (
    DynamicFBA, get_iAF1260b_config, print_growth
)
from vivarium_bioscrape.processes.bioscrape import Bioscrape
from vivarium.processes.divide_condition import DivideCondition
from vivarium.processes.meta_division import MetaDivision

# plotting
from vivarium.plots.simulation_output import plot_simulation_output


def get_metabolism_config(
    volume=1e-5*units.L,
):
    config = get_iAF1260b_config()
    config.update({
        'bin_volume': volume})
    return config


# Flux deriver
class FluxDeriver(Deriver):
    defaults = {'time_step': 1}
    def __init__(self, parameters=None):
        super(FluxDeriver, self).__init__(parameters)
    def initial_state(self, config=None):
        return {}
    def ports_schema(self):
        return {
            'deltas': {
                '*': {'_default': 0.0}},
            'fluxes': {
                '*': {'_default': 0.0}}}
    def next_update(self, timestep, states):
        deltas = states['deltas']
        return {
            'fluxes': {
                reaction_id: delta/self.parameters['time_step']
                for reaction_id, delta in deltas.items()}}


# CRN-COBRA composite
class CRN_COBRA(Generator):
    name = 'crn_cobra'
    defaults = {
        'agent_id': np.random.randint(0, 100),
        'flux_names': {
            'Lactose_internal': 'EX_lcts_e',
            'Glucose_internal': 'EX_glc__D_e',
        },
        'name_map': {
            'Lactose_external': 'lcts_e',
            'Glucose_external': 'glc__D_e',
        },
        'crn': {
            'sbml_file': 'paper/LacOperon_simple.xml',
        },
        'cobra': get_metabolism_config(),
        'divide_condition': {
            'threshold': 2000 * units.fg
        },
        'agents_path': ('..', '..', 'agents',),
        'daughter_path': tuple(),
        '_schema': {
            'cobra': {
                'flux_bounds': {
                    '*': {
                        '_emit': True
                    }
                }
            }
        }
    }

    def __init__(self, config):
        super(CRN_COBRA, self).__init__(config)

        self.processes = {
            'crn': Bioscrape(self.config['crn']),
            'cobra': DynamicFBA(self.config['cobra'])}

    def generate_processes(self, config):
        # configure derivers
        flux_config = {
            'time_step': self.processes['crn'].local_timestep()}
        division_config = dict(
            daughter_path=config['daughter_path'],
            agent_id=config['agent_id'],
            compartment=self)

        derivers = {
            'flux_deriver': FluxDeriver(flux_config),
            #             'divide_condition': DivideCondition(config['divide_condition']),
            #             'globals_deriver': DeriveGlobals({}),
            #             'division': MetaDivision(division_config),
        }

        return {**self.processes, **derivers}

    def generate_topology(self, config):
        name_map = config['name_map']
        flux_names = config['flux_names']

        # map species names
        species_names = self.processes['crn'].get_species_names()
        species_mapping = {}
        for species in species_names:
            if 'external' in species:
                species_mapping[species] = ('external', name_map.get(species, species))
            else:
                species_mapping[species] = ('internal', name_map.get(species, species))

        # exchange deltas connect to flux_bounds in metabolism
        flux_mapping = {}
        for flux in species_names:
            if flux in flux_names.keys():
                flux_mapping[flux] = ('flux_bounds', flux_names.get(flux, flux))
            else:
                flux_mapping[flux] = ('delta_species', flux)

        return {
            'crn': {
                'species': species_mapping,
                'delta_species': flux_mapping,
                'rates': ('rates',),
                'globals': ('globals',),
            },
            'cobra': {
                'internal_counts': ('internal_counts',),
                'external': ('external',),
                'exchanges': ('exchanges',),
                'reactions': ('reactions',),
                'flux_bounds': ('flux_bounds',),
                'global': ('global',),
            },
            'flux_deriver': {
                'deltas': ('delta_species',),  # TODO -- need to add units!
                'fluxes': ('flux_bounds',),
            },
            # 'globals_deriver': {
            #     'global': ('boundary',)
            # },
            # 'divide_condition': {
            #     'variable': ('boundary', 'mass',),
            #     'divide': ('boundary', 'divide',),
            # },
            # 'division': {
            #     'global': ('boundary',),
            #     'agents': config['agents_path'],
            # }
        }



# experiments
def run_bioscrape(
        total_time=2500,
        time_step=1,
):
    # initialize Bioscrape process
    bioscrape_config = {
        'sbml_file': 'paper/LacOperon_simple.xml',
        'time_step': time_step}
    bioscrape_process = Bioscrape(bioscrape_config)

    # initial state
    initial_state = bioscrape_process.initial_state()

    # run simulation
    settings = {
        'total_time': total_time,
        'initial_state': initial_state,
        'display_info': False,
        'progress_bar': False}
    timeseries = simulate_process_in_experiment(bioscrape_process, settings)

    return timeseries


def run_metabolism(
        total_time=250,
        time_step=10,
        volume=1e-5 * units.L,
):
    # configure cobra process
    cobra_config = get_iAF1260b_config()
    cobra_config.update({
        'time_step': time_step,
        'bin_volume': volume})
    metabolism = DynamicFBA(cobra_config)

    # initial state
    initial_config = {}
    initial_state = metabolism.initial_state(
        config=initial_config)

    # run simulation
    sim_settings = {
        'initial_state': initial_state,
        'total_time': total_time}
    timeseries = simulate_process_in_experiment(metabolism, sim_settings)
    print_growth(timeseries['global'])

    return timeseries


def run_crn_cobra(
        total_time=200,
        time_step=10,
        agent_id='1',
        volume=1e-5 * units.L,
):

    ## initialize  he composite
    # COBRA config
    cobra_config = get_iAF1260b_config()
    cobra_config.update({
        'time_step': time_step,
        'bin_volume': volume})
    # CRN config
    crn_config = {'time_step': time_step}
    # flux map
    flux_map = {
        'Lactose_internal': 'EX_lcts_e',
        'Glucose_internal': 'EX_glc__D_e'
    }

    composite_config = {
        'agent_id': agent_id,
        'flux_map': flux_map,
        'crn': crn_config,
        'cobra': cobra_config}

    composite = CRN_COBRA(composite_config)

    # get initial state
    initial_state = {
        'agents': {
            '1': composite.initial_state()
        }
    }

    # make the experiment
    sim_settings = {
        'initial_state': initial_state,
        'outer_path': ('agents', agent_id,)
    }
    experiment = compartment_in_experiment(composite, sim_settings)

    # run and retrieve the data
    experiment.update(total_time)
    timeseries1 = experiment.emitter.get_timeseries()
    timeseries = timeseries1['agents']['1']
    timeseries['time'] = timeseries1['time']
    print_growth(timeseries['global'])
    return timeseries


# plotting

# TODO -- Control should handle figure saving.  run_plots should just return the fig.
# def save_figure(out_dir, filename=None):
#     column_width = 3
#
#     import ipdb; ipdb.set_trace()
#     if out_dir:
#         os.makedirs(out_dir, exist_ok=True)
#         if filename is None:
#             filename = 'simulation'
#         # save figure
#         fig_path = os.path.join(out_dir, filename)
#         plt.subplots_adjust(wspace=column_width/3, hspace=column_width/3)
#         plt.savefig(fig_path, bbox_inches='tight')


def plots_1(data, config, out_dir='out'):
    plot_simulation_output(
        data,
        settings={},
        out_dir=out_dir,
    )


# plotting function for metabolism output
def plot_metabolism(data, config, out_dir='out'):

    ncol = config.get('ncol', 2)
    original_fontsize = plt.rcParams['font.size']
    plt.rcParams.update({'font.size': 9})

    # initialize subplots
    n_rows = 2
    n_cols = 2
    fig = plt.figure(figsize=(n_cols * 7, n_rows * 3))
    grid = plt.GridSpec(n_rows, n_cols)

    time_vec = data['time']

    # mass
    ax = fig.add_subplot(grid[0, 0])
    ax.plot(time_vec, data['global'][('mass', 'femtogram')], label='mass')
    ax.set_title('total compartment mass (fg)')
    ax.set_xlabel('time (sec)')
    #     ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), ncol=ncol)

    # external
    ax = fig.add_subplot(grid[0, 1])
    for mol_id, series in data['external'].items():
        if sum(series) != 0.0:
            ax.plot(time_vec, series, label=mol_id)
    ax.set_title('external concentrations (log)')
    ax.set_yscale('log')
    ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), ncol=ncol)

    # internal
    ax = fig.add_subplot(grid[1, 1])
    for mol_id, series in data['internal_counts'].items():
        if sum(series) != 0.0:
            ax.plot(time_vec, series, label=mol_id)
    ax.set_title('internal molecule counts (log)')
    ax.set_xlabel('time (sec)')
    ax.set_yscale('log')
    fig.tight_layout()
    plt.rcParams.update({'font.size': original_fontsize})


    fig_path = os.path.join(out_dir, 'metabolism')
    plt.savefig(fig_path, bbox_inches='tight')
    # save_figure(out_dir, 'metabolism')


experiments_library = {
    '1': run_bioscrape,
    '2': {
        'experiment': run_metabolism,
        'config': {}
    },
    '3': run_crn_cobra,
}
plots_library = {
    '1': {
        'plot': plots_1,
        'config': {},
    },
    'plot_metabolism': {
        'plot': plot_metabolism,
        'config': {},
    },
}

workflow_library = {
    '1': {
        'name': 'crn_alone',
        'experiment': '1',
        'plots': ['1'],
    },
    '2': {
        'name': 'cobra_alone',
        'experiment': '2',
        'plots': [
            '1',
            'plot_metabolism',
        ],
    },
    '3': {
        'name': 'crn_cobra_composite',
        'experiment': '3',
        'plots': [
            '1',
            'plot_metabolism',
        ],
    },
}

if __name__ == '__main__':
    Control(
        experiments=experiments_library,
        plots=plots_library,
        workflows=workflow_library,
        out_dir='out/experiments',
        )
