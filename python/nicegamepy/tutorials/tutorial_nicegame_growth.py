from nicegamepy.nicegame import GapFiller_nicegame
from pytfa.io import import_matlab_model
import pandas as pd
import os

# We load the model, database and thermo database (optionally)
model_path = "" # a model that does not grow under the desired media conditions
dbmodel_path = "" # your reaction database
#thermo_data_path = "./data/thermo_data.thermodb"

#thermo_data = load_thermoDB(thermo_data_path)
model = import_matlab_model(model_path)
output_path = f'output_{model.description}'
os.makedirs(output_path, exist_ok=True)

# we change the bounds in the dbmodel and later in the model to make them consistent with the bigM values
# we are using in the gap-filling constraints
dbmodel = import_matlab_model(dbmodel_path)
for rxn in dbmodel.reactions:
    if rxn.lower_bound < 0:
        rxn.lower_bound = -50
    if rxn.upper_bound > 0:
        rxn.upper_bound = 50

#################################################################################################################################
# set solver
CPLEX = 'optlang-cplex'
solver = CPLEX
model.solver = solver
model.solver.problem.parameters.reset()
model.solver.configuration.tolerances.feasibility = 1e-9
model.solver.configuration.tolerances.optimality = 1e-9
model.solver.configuration.tolerances.integrality = 1e-9
model.solver.problem.parameters.read.scale = -1
model.solver.problem.parameters.emphasis.numerical = 1
model.solver.problem.parameters.emphasis.mip.set(3)
model.solver.configuration.presolve = True

# add demand reaction for the biomass building blocks and check which ones cannot be produced
bbbs = model.reactions.get_by_id("Biomass_rxn").reactants
num_rxns = len(model.reactions)

data = []

for i in bbbs:
    model.add_boundary(model.metabolites.get_by_id(i.id), type="demand")

####################################################################################################################################
# identify the BBBs that cannot be produced in the draft model
for i, bbb_i in enumerate(bbbs):
    model.objective = model.reactions[-i-1].id
    sol = model.optimize()
    data.append([model.reactions[-i-1].id.replace('DM_', ''), sol.objective_value])

df_results = pd.DataFrame(data, columns=['bbb', 'objective_value'])
# problematic bbbs to gap fill for
df_results = df_results[df_results['objective_value'] < 1e-6]

output_file = os.path.join(output_path, f'problematic_bbbs_{model.description}.csv')
df_results['bbb'].to_csv(output_file, index=False)

# Extract reaction IDs as a list
problematic_bbbs = df_results['bbb'].tolist()

#################################################################################################################################
# merge model with database and check for which bbbs we can gap fill for
# Get the reaction IDs from the other model
# in case the model is build with the database used for gap-filling this will make mering must faster
reaction_ids_to_remove = [rxn.id for rxn in model.reactions if rxn.id in dbmodel.reactions]
# Remove matching reactions from dbmodel
dbmodel.remove_reactions(reaction_ids_to_remove)

for rxn in model.reactions:
    if rxn.lower_bound < 0:
        rxn.lower_bound = -50
    if rxn.upper_bound > 0:
        rxn.upper_bound = 50
model.reactions.get_by_id("ATPM").lower_bound = 8.39

merged_model = GapFiller_nicegame(
    model=model,
    dbmodel=dbmodel,
    imposeThermodynamics=False,
    run_gapfilling=False,
    check_reaction_names_only=True
)
######################################################################################################################
# check which of the BBBs that could not be produced in the draft model are now produced in the merged model
data = []
merged_model_i = merged_model.merged_model_tfa
for i, bbb_i in enumerate(problematic_bbbs):
    merged_model_i.objective = merged_model_i.reactions.get_by_id(f'DM_{bbb_i}').id
    sol = merged_model_i.optimize()
    data.append([bbb_i, sol.objective_value])

df_results_forgf = pd.DataFrame(data, columns=['bbb_rescued', 'objective_value'])
df_results_forgf = df_results_forgf[df_results_forgf['objective_value'] > 1e-3]

output_file = os.path.join(output_path, f'rescued_bbbs_{model.description}.csv')
df_results_forgf['bbb_rescued'].to_csv(output_file, index=False)

# Extract reaction IDs as a list
rescued_bbbs = df_results_forgf['bbb_rescued'].tolist()
##############################################################################################################################
# now let's gapfill
merged_model_tfa = merged_model.merged_model_tfa

gapfiller = GapFiller_nicegame.__new__(GapFiller_nicegame)
gapfiller.min_plus_k = 3
gapfiller.max_alternatives = 25
gapfiller.solver_timeout = 600

for i, bbb_i in enumerate(rescued_bbbs):
   # force lower bound for the production of the BBB you want to gap fill for
    merged_model_tfa.reactions.get_by_id(f'DM_{bbb_i}').lower_bound = 1e-3

    tmodel_with_constraints, bfuse_vars = gapfiller.add_constraints_gapfilling(merged_model_tfa)
    df_results, constrained_tmodel = gapfiller.alternatives_gapfilling(tmodel_with_constraints, bfuse_vars)
   # remove the lower bound for the next iteration
    merged_model_tfa.reactions.get_by_id(f'DM_{bbb_i}').lower_bound = 0

    print(bbb_i,df_results)
    # Save results to CSV
    output_file = os.path.join(output_path, f'gapfilling_results_{model.description}_{bbb_i}.csv')
    df_results.to_csv(output_file, index=False)
