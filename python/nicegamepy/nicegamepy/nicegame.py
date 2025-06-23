import logging
from typing import Optional, TYPE_CHECKING
import pandas as pd
import pytfa
from pytfa import ThermoModel
from pytfa.io import load_thermoDB
from cobra.core import Model
from pytfa.optim.variables import ForwardBackwardUseVariable, ModelVariable, ReactionVariable
from pytfa.optim.constraints import ReactionConstraint, ModelConstraint, LinearizationConstraint
from pytfa.optim.utils import symbol_sum

if TYPE_CHECKING:
    from cobra import Reaction

# Set up logger for debugging and tracking
logger = logging.getLogger(__name__)

class ActiveReaction(ReactionConstraint):
    """
    Class to represent a constraint the implements
    """

    def __init__(self, reaction, expr, **kwargs):
        ReactionConstraint.__init__(self, reaction, expr, **kwargs)

    prefix = 'AR_'

class ActiveReaction2(ReactionConstraint):
    """
    Class to represent a constraint the implements
    """

    def __init__(self, reaction, expr, **kwargs):
        ReactionConstraint.__init__(self, reaction, expr, **kwargs)

    prefix = 'AR2_'

class GapFiller_nicegame:
    def __init__(
            self,
            model: Model,  # The model we are working on
            dbmodel: Model,  # The database model (DB model)
            thermo_data: Optional[dict] = None,  # Thermo data is a dict or None
            integer_threshold: float = 1e-9,
            imposeThermodynamics: Optional[bool] = None,
            min_plus_k: int = 0,  # The user-defined number threshold for the minimal size increase
            max_alternatives: int = 1,  # The maximum number of alternative constraints to add
            run_gapfilling: bool = True,
            solver_timeout: int = 600,
            check_reaction_names_only = True,
            **kwargs,  # Any other keyword arguments
    ) -> None:
        # Store the original model and create a copy
        self.original_model = model
        self.model = model.copy()
        self.dbmodel = dbmodel  # Store the database model
        self.integer_threshold = integer_threshold
        self.min_plus_k = min_plus_k
        self.max_alternatives = max_alternatives
        self.solver_timeout = solver_timeout
        self.check_reaction_names_only = check_reaction_names_only

        # Determine whether to impose thermodynamic constraints
        if thermo_data is None:
            print("No thermodynamic database was provided, continuing without imposing thermodynamic constraints.")
            self.imposeThermodynamics = False
            self.thermo_data = {
                'name': "empty_db",
                'units': "kcal/mol",
                'metabolites': [],
                'cues': []
            }
        else:
            self.imposeThermodynamics = True if imposeThermodynamics is None else imposeThermodynamics
            self.thermo_data = thermo_data

        # Merge the DB model into the original model
        self.merged_model = self.merge_db_model(self.model,self.dbmodel)

        # Initialize the TFA model
        self.merged_model_tfa = ThermoModel(self.thermo_data, self.merged_model) if self.imposeThermodynamics else ThermoModel({'name': "empty_db", 'units': "kcal/mol", 'metabolites': [], 'cues': []}, self.merged_model)

       # if self.imposeThermodynamics:
        self.merged_model_tfa.prepare()  # Prepare the model for TFA
        self.merged_model_tfa.convert()  # Add thermodynamic displacement constraints

        # set solver

        CPLEX = 'optlang-cplex'
        solver = CPLEX
        self.merged_model_tfa.solver = solver
        self.merged_model_tfa.solver.problem.parameters.reset()
        self.merged_model_tfa.solver.configuration.tolerances.feasibility = 1e-9
        self.merged_model_tfa.solver.configuration.tolerances.optimality = 1e-9
        self.merged_model_tfa.solver.configuration.tolerances.integrality = 1e-9
        # self.merged_model_tfa.solver.problem.parameters.read.scale = -1
        #self.merged_model_tfa.solver.problem.parameters.emphasis.numerical = 1
        # self.merged_model_tfa.solver.problem.parameters.emphasis.mip.set(3)
        # self.merged_model_tfa.solver.configuration.presolve = True

        # Run gap-filling if requested
        if run_gapfilling:
            self.tmodel_gf_constrained, self.bfuse_variables = self.add_constraints_gapfilling(self.merged_model_tfa)
            self.alternatives_results, self.tmodel_gf_alt_constrained = self.alternatives_gapfilling(self.tmodel_gf_constrained, self.bfuse_variables)
        else:
            self.alternatives_results = None
            self.tmodel_gf_constrained = None

    def merge_db_model(self, merged_model,dbmodel):
        """
        Merge the database model (DB model) into the working model.
        Each reaction is given an 'origin' attribute to track whether it comes from the DB model or the source model.
        """

        metabolite_id = ['C00080_c', 'cpd00067_c','h_c']

        # Now compare and add reactions, ignoring 'C00080_c' from reactants and products
        for rxn in dbmodel.reactions:
            new_rxn = rxn.copy()
            if self.check_reaction_names_only:
                # Check only by reaction ID
                reaction_exists = any(existing_rxn.id == new_rxn.id for existing_rxn in merged_model.reactions)
                if reaction_exists: # here there should not be any match if we remove them same names in the main script
                    print(f"Reaction {new_rxn.id} already exists by ID. Not adding.")
                else:
                    new_rxn.id = f"{new_rxn.id}_db"  # Modify the reaction ID to indicate it's from the DB model
                    merged_model.add_reactions([new_rxn])
                    new_rxn.origin = "db"  # Mark the new reaction as coming from the DB model
            else:
                # Get the reactants and products of the new reaction, excluding C00080_c
                new_reactants = sorted([met.id for met in new_rxn.reactants if met.id not in metabolite_id])
                new_products = sorted([met.id for met in new_rxn.products if met.id not in metabolite_id])

                # Check if a reaction with the same reactants and products already exists in the model
                reaction_exists = False
                for existing_rxn in merged_model.reactions:
                    # Get the reactants and products of the existing reaction, excluding C00080_c
                    existing_reactants = sorted([met.id for met in existing_rxn.reactants if met.id != metabolite_id])
                    existing_products = sorted([met.id for met in existing_rxn.products if met.id != metabolite_id])

                    # Compare both reactants and products in both directions (direct and reversed)
                    if (new_reactants == existing_reactants and new_products == existing_products) or \
                            (new_reactants == existing_products and new_products == existing_reactants):
                        reaction_exists = True
                        print(f"Reaction {new_rxn.id} is the same as existing reaction {existing_rxn.id}. Not adding it.")
                        break  # Stop the loop if we find a match

                # If no matching reaction exists, add the new reaction
                if not reaction_exists:
                    new_rxn.id = f"{new_rxn.id}_db"  # Modify the reaction ID to indicate it's from the DB model
                    merged_model.add_reactions([new_rxn])
                    new_rxn.origin = "db"  # Mark the new reaction as coming from the DB model
                        #print(f"Reaction {new_rxn.id} added to the model.")

            # Mark all reactions without the 'origin' attribute as 'source'
            for rxn in merged_model.reactions:
                if not hasattr(rxn, 'origin'):
                    rxn.origin = "source"
        return merged_model


    def add_constraints_gapfilling(self, tmodel):
        """
        Add gap-filling constraints to reactions that originate from the DB model.
        These constraints enforce that forward + reverse flux is conditional on a binary use variable.
        """
        bfuse_variables = []
        for rxn in tmodel.reactions:
            if rxn.origin == 'db':

                F_ = tmodel.reactions.get_by_id(rxn.id).forward_variable
                B_ = tmodel.reactions.get_by_id(rxn.id).reverse_variable
                BFUSE = tmodel.add_variable(ForwardBackwardUseVariable,rxn)
                bfuse_variables += [BFUSE]

                expression = F_+B_+100.0*BFUSE
                tmodel.add_constraint(
                    ActiveReaction,
                    rxn,
                    expression,
                    lb=0.0,
                    ub=100.0
                )

                expression = F_+B_+BFUSE
                tmodel.add_constraint(
                    ActiveReaction2,
                    rxn,
                    expression,
                    lb=1e-7,
                    ub=50.0
                )

        expression_objective = symbol_sum(bfuse_variables)
        tmodel.objective = expression_objective
        tmodel.objective_direction = 'max'
        return tmodel, bfuse_variables

    def alternatives_gapfilling(self, tmodel,bfuse_variables):
        """
        Add alternative constraints until no feasible solution is found,
        minimal solution size exceeds initial + k, or max_alternatives is reached.
        """
        from optlang.interface import INFEASIBLE

        CPLEX = 'optlang-cplex'
        solver = CPLEX
        tmodel.solver = solver
        tmodel.solver.problem.parameters.reset()
        tmodel.solver.configuration.timeout = self.solver_timeout
        tmodel.solver.configuration.tolerances.feasibility = 1e-9
        tmodel.solver.configuration.tolerances.optimality = 1e-9
        tmodel.solver.configuration.tolerances.integrality = 1e-9
        # tmodel.solver.problem.parameters.read.scale = -1
       # tmodel.solver.problem.parameters.emphasis.numerical = 1
       #  tmodel.solver.problem.parameters.emphasis.mip.set(3)
       #  tmodel.solver.configuration.presolve = True

        results = []
        sol = tmodel.optimize()

        BFUSE_var_names = [x.name for x in bfuse_variables]

        df_sol = sol.raw
        df_bfuse = df_sol.loc[BFUSE_var_names]
        idx_active_bfuse = df_bfuse[df_bfuse < 0.5].index
        active_vars = [tmodel._var_dict[k] for k in idx_active_bfuse]
        active_rxn_id = [i.reaction.id for i in active_vars]
        initial_min_size = len(active_rxn_id)
        print("Min Size:")
        print(initial_min_size)
        alternatives_added = 0
        while True:

            reaction_formulas_ids = []
            reaction_formulas_names = []
            reaction_directions = []

            for var in active_vars:
                rxn = var.reaction
                # Formula with metabolite IDs
                formula_ids = rxn.reaction
                # Formula with metabolite names
                formula_names = rxn.build_reaction_string(use_metabolite_names=True)

                # Determine directionality
                Flux_value = sol.get_primal_by_id(rxn.id)
                if Flux_value > 0:
                    direction = 'F'
                elif Flux_value < 0:
                    direction = 'R'

                reaction_formulas_ids.append(formula_ids)
                reaction_formulas_names.append(formula_names)
                reaction_directions.append(direction)

            results.append({
                'iteration': alternatives_added,
                'min_size': len(active_rxn_id),
                'active_rxn_ids': active_rxn_id,
                'reaction_formulas_ids': reaction_formulas_ids,
                'reaction_formulas_names': reaction_formulas_names,
                'directions': reaction_directions,
            })
            #
            # results.append({
            #     'iteration': alternatives_added,
            #     'min_size': len(active_rxn_id),
            #     'active_rxn_ids': active_rxn_id,
            # })

            if len(active_rxn_id) > initial_min_size + self.min_plus_k:
                print("Size limit has been reached")
                break
            if alternatives_added >= self.max_alternatives:
                print("Alternative number limit has been reached")
                break

            constraint_expr = symbol_sum([var for var in active_vars])
            constraint_id = f'BFUSE_sum_constraint_{alternatives_added}'
            tmodel.add_constraint(
                ModelConstraint,
                tmodel,
                constraint_expr,
                lb=0.5,
                ub=len(bfuse_variables),
                id_=constraint_id
            )
            try:
                sol = tmodel.optimize()
            except Exception:
                print("No more solutions found")
                break

            df_sol = sol.raw
            df_bfuse = df_sol.loc[BFUSE_var_names]
            idx_active_bfuse = df_bfuse[df_bfuse < 0.5].index
            active_vars = [tmodel._var_dict[k] for k in idx_active_bfuse]
            active_rxn_id = [i.reaction.id for i in active_vars]

            alternatives_added += 1
            print("Alternative: ",alternatives_added)

        return pd.DataFrame(results), tmodel

