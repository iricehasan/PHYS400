from sympy.physics.quantum import Commutator, Operator
import sympy as sp
from itertools import permutations, product
from fractions import Fraction
from bch import (
    calculate_phi_n,
    calculate_big_phi_m_for_2_operators,
    calculate_big_phi_m,
    bch_expansion_up_to_order,
    convert_coeffs_to_fractions
)

def group_terms_by_commutator(expr, commutators):
    terms = expr.as_ordered_terms()
    grouped_terms = {comm: 0 for comm in commutators}
    
    for term in terms:
        for comm in commutators:
            if comm in term.args:
                grouped_terms[comm] += term.coeff(comm)
                break
        else:
            if None in grouped_terms:
                grouped_terms[None] += term
            else:
                grouped_terms[None] = term

    return grouped_terms

T = Operator('T')
V = Operator('V')
t = sp.symbols('t', commutative = True)

order = 3  # Example order, adjust as necessary
Z = bch_expansion_up_to_order([t/2*T, t*V, t/2*T], order)
Z = convert_coeffs_to_fractions(Z)

# Generate alpha and beta symbols for the given order
alpha = sp.symbols(f'alpha1:{order+1}')
beta = sp.symbols(f'beta1:{order+1}')

# Create the list of operators for \Pi exp(alpha_i t T) exp(beta_i t V)
expansion_operators = []
for a, b in zip(alpha, beta):
    expansion_operators.append(a * t * T)
    expansion_operators.append(b * t * V)

# Compute Z' for \Pi exp(alpha_i t T) exp(beta_i t V)
Z_prime = bch_expansion_up_to_order(expansion_operators, order)
Z_prime = convert_coeffs_to_fractions(Z_prime)

# Expand and collect terms
Z_collected = sp.collect(Z.expand(commutator = True), t)

Z_prime_collected = sp.collect(Z_prime.expand( commutator = True), t)
sp.pprint(Z_prime_collected, use_unicode=True)

# Extract coefficients for each power of t
coeff_Z = {}
coeff_Z_prime = {}
for i in range(1, order+1):
    coeff_Z[f't{i}'] = Z_collected.coeff(t, i)
    coeff_Z_prime[f't{i}'] = Z_prime_collected.coeff(t, i)

# this could be obtained from a function as a new feature
commutator_terms = [Commutator(T, V), Commutator(Commutator(T, V), T), Commutator(Commutator(T, V), V), Commutator(Commutator(V, V), T), Commutator(Commutator(V, V), V) ]

grouped_terms_Z = {}
grouped_terms_Z_prime = {}
for i in range(2, order+1):
    grouped_terms_Z[f't{i}'] = group_terms_by_commutator(coeff_Z[f't{i}'], commutator_terms)
    grouped_terms_Z_prime[f't{i}'] = group_terms_by_commutator(coeff_Z_prime[f't{i}'], commutator_terms)

eq1 = sp.Eq(coeff_Z['t1'].coeff(T), coeff_Z_prime['t1'].coeff(T))
eq2 = sp.Eq(coeff_Z['t1'].coeff(V), coeff_Z_prime['t1'].coeff(V))

equations = []

for i in range(2, order+1):
    for comm in grouped_terms_Z_prime[f't{i}'].keys():
        if comm in grouped_terms_Z[f't{i}']:
            eq = sp.Eq(grouped_terms_Z[f't{i}'][comm], grouped_terms_Z_prime[f't{i}'][comm])
            equations.append(eq)

equations.append(eq1)
equations.append(eq2)

# filter equations by removing True values
equations = [eq for eq in equations if eq != True]

variables = alpha + beta

solutions = {}
remaining_equations = equations

while remaining_equations:
    # Solve for the next variable
    result = sp.solve(remaining_equations, variables, dict=True)
    
    if not result:
        break  # Exit if no solutions are found
    
    # If result is a list of dictionaries, iterate through each solution
    if isinstance(result, list):
        for res in result:
            if isinstance(res, dict):
                solutions.update(res)
            else:
                print("Unexpected result format:", res)
                break
    elif isinstance(result, dict):
        solutions.update(result)
    else:
        print("Unexpected result format:", result)
        break

    # Generate new equations from unsolved variables
    remaining_equations = []
    for res in result:
        for var in variables:
            if var in res:
                remaining_equations.append(sp.Eq(var, res[var]))
    
    # Substitute found solutions into remaining equations
    remaining_equations = [eq.subs(solutions) for eq in remaining_equations]
    remaining_equations = [eq for eq in remaining_equations if eq != True]  # Remove solved equations
    
    # Update variables for the next iteration
    variables = [var for var in variables if var not in solutions]

for var in solutions:
    solutions[var] = solutions[var].subs(solutions)

# Check if all variables are solved
for var in list(alpha) + list(beta):
    if var not in solutions:
        print(f"Warning: {var} was not solved.")

# Print the final solutions
print("coefficients: ", solutions)
