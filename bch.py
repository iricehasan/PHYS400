import sympy as sp
from sympy.physics.quantum import Commutator
from itertools import permutations, product
from math import factorial
from fractions import Fraction

# calculates the descent number in permutations
def descent_count(m, perm):
    """Count the number of descents in the permutation with a dynamic hierarchy."""
    hierarchy = {f"X{i+1}": i + 1 for i in range(len(perm))}  
    return sum(1 for i in range(len(perm) - 1) if hierarchy[str(perm[i])] > hierarchy[str(perm[i + 1])])

def right_nested_commutator(elements):
    # Compute right-nested commutator from a list of elements.
    if len(elements) == 1:
        return elements[0]
    else:
        return Commutator(elements[0], right_nested_commutator(elements[1:]))

def binomial_coeff(n, k):
    # Calculate binomial coefficient.
    return factorial(n) // (factorial(k) * factorial(n - k))

# Replaces the X_1, X_2, X_3 terms with operators given in the variable in usage
def replace_terms(terms, original, indexed):
    # Replace indexed terms with original terms in the expression.
    replacement_dict = {indexed[i]: original[i] for i in range(len(original))}
    #print(replacement_dict)
    return [term.subs(replacement_dict) for term in terms]

# equation 12 in the paper
def calculate_phi_n(elements):
    # Calculate the phi_n using the given elements.
    n = len(elements)
    sum_term = sp.S.Zero  # sympy expression for zero
    terms = []
    for perm in permutations(elements):
        d = descent_count(n, perm)
        commutator = right_nested_commutator(list(perm))
        coefficient = (-1)**d / binomial_coeff(n - 1, d)
        term = coefficient * commutator
        terms.append(term / n**2)
    return terms

# equation 10 in the paper
def calculate_big_phi_m_for_2_operators(X, Y, m):
    # Calculate the Phi_m(X, Y) for given m using the equation 10.
    sum_term = sp.S.Zero  # Initialize the sum as a sympy expression for zero
    # Sum over all i and j such that i + j = m and i, j >= 1
    all_terms = []
    for i in range(1, m):
        j = m - i
        terms = [X] * i + [Y] * j
        indexed_terms = [sp.Symbol(f'X{i+1}', commutative=False) for i in range(len(terms))] # for the calculation of descents
        phi_n_terms = calculate_phi_n(indexed_terms)
        replaced_terms = replace_terms(phi_n_terms, terms, indexed_terms)
        facc = 1 / (factorial(i) * factorial(j))
        for term in replaced_terms:
            all_terms.append(term / (factorial(i) * factorial(j)))
    return all_terms

# general formula given in equation 7
def calculate_big_phi_m(operators, m):
    all_terms = []
    num_operators = len(operators)

    for i in find_combinations(m, num_operators):
        terms = []
        indexed_terms = []
        for op_index, count in enumerate(i):
            terms.extend([operators[op_index]] * count)
            indexed_terms.extend([sp.Symbol(f'{str(operators[op_index])}{i+1}', commutative=False) for i in range(count)])

        indexed_terms = [sp.Symbol(f'X{i+1}', commutative=False) for i in range(len(terms))]
        phi_n_terms = calculate_phi_n(indexed_terms)
        replaced_terms = replace_terms(phi_n_terms, terms, indexed_terms)

        facc = 1 / sp.prod(factorial(x) for x in i)

        for term in replaced_terms:
            all_terms.append(term * facc)
            
    return all_terms

# function to get i_1, i_2,i_3 terms in equation 7 in the paper
def find_combinations(total, split_number, current=None):
    """
    Recursively find all ordered combinations of positive integers that sum to a given total,
    with exactly 'split_number' parts.
    """
    if current is None:
        current = []

    if len(current) == split_number:
        if sum(current) == total:
            yield current
        return
    
    # Start from 0 for each recursive call
    start = 0
    # The maximum value we consider is total - sum of current + 1
    # This is because we need at least `split_number - len(current) - 1` more parts to complete the combination
    max_val = total - sum(current) + 1
    for i in range(start, max_val):
        # Continue to the next part of the split
        yield from find_combinations(total, split_number, current + [i]) 

# function to convert coefficients to fractions to better visibility
def convert_coeffs_to_fractions(expr):
    if expr.func == sp.Mul and expr.args[0].is_Float:
        # Convert the first argument if it is a float (coefficient in Mul expressions)
        fraction = Fraction(str(expr.args[0])).limit_denominator()
        # Recreate the multiplication with the fraction and remaining arguments
        return sp.Mul(fraction, *expr.args[1:], evaluate=False)
    else:
        # Recursively apply this to all arguments in other expressions
        return expr.func(*(convert_coeffs_to_fractions(arg) for arg in expr.args), evaluate=False)

# Example usage
m = 3  # For example, calculate Φ3(X, Y)

operators = sp.symbols('X Y', commutative=False) 

for m in range(1,4):
    phi_3 = calculate_big_phi_m(operators, m)
    phi_3_summed = sum(phi_3)
    converted_expression = convert_coeffs_to_fractions(phi_3_summed)
    print(f"Φ{m}(X, Y) =", phi_3_summed, "\n")
    sp.pprint(converted_expression, use_unicode=True)

