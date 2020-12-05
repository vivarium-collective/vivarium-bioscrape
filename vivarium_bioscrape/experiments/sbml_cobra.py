from vivarium.core.control import Control
from vivarium.core.composition import simulate_compartment_in_experiment
from vivarium_cobra.processes.metabolism import Metabolism
from vivarium_bioscrape.processes.bioscrape import Bioscrape


def run_experiment():
    bioscrape = Bioscrape({'sbml_file': 'models/BIOMD0000000065_url.xml'})
    initial_state = bioscrape.initial_state()
    settings = {
        'total_time': 30,
        'initial_state': initial_state,
        }
    return simulate_compartment_in_experiment(bioscrape, settings)


experiments_library = {
    '1': run_experiment,
}
plots_library = {}
workflow_library = {}

if __name__ == '__main__':
    Control(
        experiments=experiments_library,
        plots=plots_library,
        workflows=workflow_library,
        )
