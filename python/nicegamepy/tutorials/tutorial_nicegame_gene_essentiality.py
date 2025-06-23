from cobra.io import load_json_model
from nicegamepy.nicegame import GapFiller_nicegame
from pytfa.io import import_matlab_model, load_thermoDB,read_lexicon, annotate_from_lexicon
import pandas as pd

# We load the model, database and thermo database (optionally)
model_path = "./models/small_ecoli.mat"
dbmodel_path = "./models/iJO1366.json"
thermo_data_path = "./data/thermo_data.thermodb"


thermo_data = load_thermoDB(thermo_data_path)
model = import_matlab_model(model_path)
dbmodel = import_matlab_model(dbmodel_path)
#dbmodel = load_json_model(dbmodel_path)

# set the extracellular environment, etc...
model.reactions.get_by_id('EX_glc__D_e').lower_bound = -10
model.reactions.get_by_id('ATPM').lower_bound = 8

for rxn in dbmodel.reactions:
    if rxn.lower_bound < 0:
        rxn.lower_bound = -50
    if rxn.upper_bound > 0:
        rxn.upper_bound = 50


# set objective
model.objective = "Ec_biomass_iJO1366_WT_53p95M"
model.reactions.get_by_id("Ec_biomass_iJO1366_WT_53p95M").lower_bound = 1e-3


# set solver
CPLEX = 'optlang-cplex'
GUROBI = 'optlang-gurobi'
GLPK = 'optlang-glpk'
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

optimal_growth = model.optimize().objective_value

########################################################################################################################################################################################################
# here we identify gaps, in this example we consider as gaps essential reactions
results_reaction_del = []

# reaction deletion
for rxn in model.reactions:
    modeli = model.copy()
    modeli.reactions.get_by_id(rxn.id).upper_bound = 0
    modeli.reactions.get_by_id(rxn.id).lower_bound = 0
    sol = modeli.optimize()
    results_reaction_del.append({'reaction': rxn.id, 'objective_value': sol.objective_value})

df = pd.DataFrame(results_reaction_del)

# Filter for essential reactions (objective drops below 10% of wild-type optimal growth)
essential_df = df[df['objective_value'] < 0.1 * optimal_growth]
df.to_csv("essential_reactions.csv", index=False)            # Comma-separated
# Extract reaction IDs as a list
essential_reactions = essential_df['reaction'].tolist()

########################################################################################################################################################################################################
# we merge the model with the DB and check which gaps can be filled,
# in this case how many of the essential reactions are not essential any more
merged_model = GapFiller_nicegame(
    model=model,
    dbmodel=dbmodel,
    imposeThermodynamics=False,
    run_gapfilling=False
)

merged_model_for_essentiality = merged_model.merged_model
merged_model_for_essentiality.reactions.get_by_id("Ec_biomass_iJO1366_WT_53p95M").lower_bound = 1e-3
results_reaction_del_m = []
# reaction deletion
for rxn in essential_reactions:
    modeli = merged_model_for_essentiality.copy()
    modeli.reactions.get_by_id(rxn).upper_bound = 0
    modeli.reactions.get_by_id(rxn).lower_bound = 0
    sol = modeli.optimize()
    results_reaction_del_m.append({'reaction': rxn, 'objective_value': sol.objective_value})

df_m = pd.DataFrame(results_reaction_del_m)

# Filter for essential reactions (objective drops below 10% of wild-type optimal growth)
essential_df_m = df_m[df_m['objective_value'] < 0.1 * optimal_growth]
# Extract reaction IDs as a list
essential_reactions_m = essential_df_m['reaction'].tolist()
rescued_reactions = list(set(essential_reactions) - set(essential_reactions_m))
df_rescued = pd.DataFrame(rescued_reactions)
df_rescued.to_csv("rescued_reactions.csv", index=False)            # Comma-separated

#######################################################################################################################################################
# we now test the gapfilling algorithm by removing one the reactions
# that is essential in the model but not in the rescued and we gapfill
# #reaction = model.reactions.get_by_id('FMNAT')
reaction = model.reactions.get_by_id('IPPMIa')

model.remove_reactions([reaction])
model.optimize().objective_value
# we set a lower bound to the objective function of the original model
model.reactions.get_by_id("BIOMASS_Ec_iML1515_WT_75p37M").lower_bound = 1e-1

gapfilling_model = GapFiller_nicegame(
    model=model,
    dbmodel=dbmodel,
    imposeThermodynamics=False,
    k=0,  # The user-defined number threshold for the minimal size increase
    max_alternatives=10,
    run_gapfilling=True
)

# Retrieve and display the DataFrame of gap-filling alternatives
df_results = gapfilling_model.alternatives_results
print(df_results)

# Optional: Save results to CSV
df_results.to_csv("gapfilling_results.csv", index=False)
#######################################################################################################################################################
target_rxns = rescued_reactions[:5]
merged_model_tfa = merged_model.merged_model_tfa
merged_model_tfa.reactions.get_by_id("BIOMASS_Ec_iML1515_WT_75p37M").lower_bound = 1e-1

gapfiller = GapFiller_nicegame.__new__(GapFiller_nicegame)  # create instance without calling __init__
gapfiller.k = 2
gapfiller.max_alternatives = 10

for r in target_rxns:
    merged_model_c = merged_model_tfa.copy()
    # Manually copy the 'origin' attribute
    for rxn_orig, rxn_copy in zip(merged_model_tfa.reactions, merged_model_c.reactions):
        if hasattr(rxn_orig, "origin"):
            rxn_copy.origin = rxn_orig.origin
    # merged_model_c.reactions.get_by_id(r).lower_bound = 0
    # merged_model_c.reactions.get_by_id(r).upper_bound = 0
    merged_model_c.remove_reactions([r])

    # just reuse these methods on your existing TFA model
    tmodel_with_constraints, bfuse_vars = gapfiller.add_constraints_gapfilling(merged_model_c)
    df_results, constrained_tmodel = gapfiller.alternatives_gapfilling(tmodel_with_constraints, bfuse_vars)

    print(df_results)
    # Save results to CSV
    df_results.to_csv(f'gapfilling_results_{r}.csv', index=False)


