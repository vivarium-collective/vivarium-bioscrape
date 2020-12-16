
import numpy as np

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
                '*': {
                    '_default': 0.0
                }
            },
            'fluxes': {
                '*': {
                    '_default': 0.0
                }
            }
        }
    def next_update(self, timestep, states):
        deltas = states['deltas']
        return {
            'fluxes': {
                reaction_id: delta/self.parameters['time_step']
                for reaction_id, delta in deltas.items()
            }
        }


# CRN-COBRA composite
class CRN_COBRA(Generator):
    name = 'crn_cobra'
    defaults = {
        'agent_id': np.random.randint(0, 100),
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

    def generate_processes(self, config):
        # TODO -- reaction names from CRN needs to be matched with COBRA

        processes = {
            'crn': Bioscrape(config['crn']),
            'cobra': DynamicFBA(config['cobra'])}

        # configure derivers
        flux_config = {
            'time_step': processes['crn'].local_timestep()}
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

        return {**processes, **derivers}

    def generate_topology(self, config):
        return {
            'crn': {
                'species': ('species',),
                'delta_species': ('delta_species',),
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
                'fluxes': ('flux_bounds_2',),
            },
            #             'globals_deriver': {
            #                 'global': ('boundary',)
            #             },
            #             'divide_condition': {
            #                 'variable': ('boundary', 'mass',),
            #                 'divide': ('boundary', 'divide',)
            #             },
            #             'division': {
            #                 'global': ('boundary',),
            #                 'agents': config['agents_path']
            #             }
        }



# experiments
def run_bioscrape(
        total_time=2500,
        time_step=1,
):
    # initialize Bioscrape process
    bioscrape_process = Bioscrape({
        'sbml_file': 'paper/LacOperon_simple.xml',
        'time_step': time_step,
    })

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
        total_time=2500,
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
        total_time=2500,
        time_step=10,
        agent_id='1',
        volume=1e-5 * units.L,
):

    # initialize composite
    cobra_config = get_iAF1260b_config()
    cobra_config.update({
        'time_step': time_step,
        'bin_volume': volume})
    crn_config = {'time_step': time_step}

    composite_config = {
        'agent_id': agent_id,
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
    timeseries = experiment.emitter.get_timeseries()
    print_growth(timeseries['agents']['1']['global'])
    return timeseries


# plotting
def plots_1(data, config, out_dir='out'):
    plot_simulation_output(
        data,
        settings={},
        out_dir=out_dir,
    )


experiments_library = {
    '1': run_bioscrape,
    '2': run_metabolism,
    '3': run_crn_cobra,
}
plots_library = {
    '1': {
        'plot': plots_1,
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
        'plots': ['1'],
    },
    '3': {
        'name': 'crn_cobra_composite',
        'experiment': '3',
        'plots': ['1'],
    },
}

if __name__ == '__main__':
    Control(
        experiments=experiments_library,
        plots=plots_library,
        workflows=workflow_library,
        out_dir='out/experiments',
        )
