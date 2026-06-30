##############################################################################################
#
# This script defines the Weak Temperature Gradient (WTG) tropical model using
# the LayerCake library
#
##############################################################################################


# Importing libraries
#####################

import numpy as np
from sympy import Symbol, sin, cos

# importing all that is needed to create a cake
from layercake import *

# importing specific modules to create the model basis of functions
from layercake.basis.centered_planar_fourier import contiguous_channel_basis
from layercake.inner_products.definition import StandardSymbolicInnerProductDefinition
from layercake.arithmetic.terms.gradient import vorticity_gradients_product
from layercake.arithmetic.terms.operations import ProductOfTerms
from layercake.arithmetic.terms.operators import OperatorTerm, ComposedOperatorsTerm


def define_model(nx, ny, chi=0.07):
    """Construct the WTG-TN model using LayerCake.

    Parameters
    ----------
    nx: int
        The truncation of the model basis in the zonal x direction.
    ny: int
        The truncation of the model basis in the meridional y direction.
    chi: float, optional
        The non-dimensional value of the WTG forcing chi_2 and chi_5 components.
    """

    # Setting some parameters
    ##########################

    # Characteristic length scale (L_y / pi)
    L_symbol = Symbol('L')
    L = Parameter(1273239.5447351628, symbol=L_symbol, units='[m]')

    # Domain aspect ratio
    n_symbol = Symbol('n')
    n = Parameter(0.2, symbol=n_symbol)

    # Meridional gradient of the Coriolis parameter at phi_0
    beta_symbol = Symbol(u'β')
    beta = Parameter(2.88e-11, symbol=beta_symbol, units='[m^-1][s^-1]')

    # Atmosphere bottom friction coefficient
    r_symbol = Symbol('r')
    r = Parameter(1.e-06, symbol=r_symbol, units='[s^-1]')

    # Mixed layer height
    H_symbol = Symbol('H')
    H = Parameter(200., symbol=H_symbol, units='[m]')

    # Gravity acceleration
    g_symbol = Symbol('g')
    g = Parameter(9.81, symbol=g_symbol, units='[m][s^-2]')

    # Forcing parameters (dimensionless for now)
    Ap = Parameter(3*np.pi*chi/4, symbol=Symbol("A'"))
    Bp = Parameter(3*np.pi*chi/4, symbol=Symbol("B'"))

    # Defining the domain
    ######################

    parameters = [n]
    atmospheric_basis = contiguous_channel_basis(nx, ny, parameters)

    # coordinates
    x = atmospheric_basis.coordinate_system.coordinates_symbol_as_list[0]
    y = atmospheric_basis.coordinate_system.coordinates_symbol_as_list[1]

    # creating an inner product definition with an optimizer for trigonometric functions
    inner_products_definition = StandardSymbolicInnerProductDefinition(coordinate_system=atmospheric_basis.coordinate_system,
                                                                       optimizer='trig', kwargs={'conds': 'none'})

    # Derived (non-dimensional) parameters
    #######################################

    # wave speed
    c = Parameter(np.sqrt(g * H), symbol=Symbol('c'), units='[m][s^-1]')

    # time units
    T = Parameter(1./np.sqrt(beta * c), symbol=Symbol('T'), units='[s]')

    # inverse of timescale for Laplacian to power 2 (nabla^4)
    T4 = Parameter(3 * 24 * 3600, symbol=Symbol('T_4'), units='[s]')
    iT4 = Parameter(T / T4, symbol=Symbol('T_4^{-1}'))

    # nondimensional bottom friction
    # rp_symbol = Symbol("r'")
    rp_symbol = Symbol("r")
    rp = Parameter(r * T, symbol=rp_symbol, units='')

    # Non-dimensional Meridional gradient of the Coriolis parameter at phi_0
    # betap_symbol = Symbol(u"β'")
    betap_symbol = Symbol(u"β")
    beta_nondim = Parameter(beta * L * T, symbol=betap_symbol, units='')

    # Defining the fields
    #######################
    p = u'ψ'
    psi = Field("psi", p, atmospheric_basis, inner_products_definition, units="[m^2][s^-2]", latex=r'\psi')
    ch = u'χ'
    ch_expression = Ap.symbol * cos(y)**2 * cos(n_symbol * x) + Bp.symbol * sin(y) * cos(n_symbol * x)
    chi = FunctionField("chi", ch, ch_expression,  atmospheric_basis, [Ap, Bp], inner_products_definition, units="[m^2][s^-2]", latex=r'\chi')

    # --------------------------------
    #
    #   Tropical field equation
    #
    # --------------------------------

    # defining the LHS as the time derivative of the vorticity
    vorticity = OperatorTerm(psi, Laplacian, atmospheric_basis.coordinate_system)
    tropical_equation = Equation(psi, lhs_terms=vorticity)

    # Defining the advection term
    advection_term = vorticity_advection(psi, psi, atmospheric_basis.coordinate_system, sign=-1)
    tropical_equation.add_rhs_terms(advection_term)

    # adding the beta term
    beta_term = OperatorTerm(psi, D, x, prefactor=beta_nondim, sign=-1)
    tropical_equation.add_rhs_term(beta_term)

    # adding the friction with the ground
    friction = OperatorTerm(psi, Laplacian, atmospheric_basis.coordinate_system, prefactor=rp, sign=-1)
    tropical_equation.add_rhs_term(friction)

    # Adding power 2 of Laplacian (nabla^4) dissipation
    operators = (Laplacian,) * 2
    operators_args = (atmospheric_basis.coordinate_system,) * 2
    lap2 = ComposedOperatorsTerm(psi, operators, operators_args, prefactor=iT4)
    tropical_equation.add_rhs_term(lap2)

    # forcing
    vorticity_gradients = vorticity_gradients_product(chi, psi, atmospheric_basis.coordinate_system, sign=-1)
    tropical_equation.add_rhs_terms(vorticity_gradients)

    beta_chi = OperatorTerm(chi, D, y, prefactor=beta_nondim, sign=-1)
    tropical_equation.add_rhs_term(beta_chi)

    chi_laplacian_term = OperatorTerm(chi, Laplacian, atmospheric_basis.coordinate_system)
    psi_laplacian_term = OperatorTerm(psi, Laplacian, atmospheric_basis.coordinate_system)
    chi_psi_laplacian_product = ProductOfTerms(psi_laplacian_term, chi_laplacian_term, sign=-1)
    tropical_equation.add_rhs_term(chi_psi_laplacian_product)

    y_expression = Expression(betap_symbol*y, latex=r"\beta y", expression_parameters=[beta_nondim])
    y_chi_laplacian_term = OperatorTerm(chi, Laplacian, atmospheric_basis.coordinate_system, prefactor=y_expression, sign=-1)
    tropical_equation.add_rhs_term(y_chi_laplacian_term)

    # --------------------------------
    #
    #   Constructing the layer
    #
    # --------------------------------

    layer = Layer()
    layer.add_equation(tropical_equation)

    # --------------------------------
    #
    #   Constructing the cake
    #
    # --------------------------------

    cake = Cake()
    cake.add_layer(layer)

    return cake

