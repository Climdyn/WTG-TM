
##############################################################################################
#
# This script generates the Weak Temperature Gradient (WTG) tropical model LaTeX equations
# using the LayerCake library
#
##############################################################################################

from model_definition import define_model
from layercake.utils.symbolic_tensor import get_coords_from_index
from sympy import Add


def convert_to_latex(tensor, field=None, number_of_terms_per_line=5):
    """Convert a model symbolic tendencies terms tensor to a list of the ODE equations in LaTeX format.

    Parameters
    ----------
    tensor: ~sympy.tensor.array.ImmutableSparseNDimArray
        Symbolic tendencies terms tensor to convert.
    field: str, optional
        Latex expression for the field of the equation.
        Default to `\\psi`.
    number_of_terms_per_line: int, optional
        Number of additive terms per line in the latex equations.
        Default to `4`.

    Returns
    -------
    list(str)
        List of ODE equations in LaTeX format.
    """
    if field is None:
        field = r'\psi'
    ndim = tensor.shape[0]
    shape_len = len(tensor.shape)
    equations_list = list()
    for i in range(ndim):
        equations_list.append(r'\dot ' + field + r'_{' + str(i) + r'} & = & ')
    for i in range(1, ndim):
        k = 0
        split_final = ""
        for n, val in tensor[i]._args[0].items():
            k += 1
            equations_list[i] += split_final
            coords = get_coords_from_index(n, ndim, shape_len-1)
            new_term = val._repr_latex_()[15:-1] + ' '
            if isinstance(val, Add):
                new_term = '+ (' + new_term + ')'
            elif new_term[0] != '-':
                new_term = '+' + new_term
            for c in coords:
                if c != 0:
                    new_term += r'\, ' + field + '_{' + str(c) + r'}'
            equations_list[i] += new_term
            if k == number_of_terms_per_line:
                k = 0
                split_final = r'\nonumber \\' + ' \n     & & '
            else:
                split_final = ""

        equations_list[i] += r' \\'

    return equations_list[1:]

# Guarding the main script to deal with multiprocessing import issue
# in case the start method is 'spawn' or 'forkserver'.
# See https://docs.python.org/3/library/multiprocessing.html#the-spawn-and-forkserver-start-methods
# 'Safe importing of main module' section for more details.
if __name__ == "__main__":

    print('Creating the model...')
    # constructing the model
    model_definition = define_model(2, 2)

    # computing the tensor (might take a long time depending on the resolution
    model_definition.compute_tensor(numerical=False, compute_inner_products=True, compute_inner_products_kwargs={'timeout': None})

    latex_equations_list = convert_to_latex(model_definition.tensor, number_of_terms_per_line=3)

    for eq in latex_equations_list:
        print(eq)
