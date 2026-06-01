from firedrake import *
from firedrake.cython import dmcommon
from firedrake.output import VTKFile
#from firedrake.petsc import PETSc
#from firedrake.constant import Constant
import ufl
#from firedrake.__future__ import interpolate
from firedrake import interpolate

import numpy# as np
from tabulate import tabulate
import argparse
## To deal with the division by zero issue using except RunTimeWarning: i.e. to catch a 
## warning as if it were an exception
import warnings
warnings.filterwarnings("error", category=RuntimeWarning)
import matplotlib.pyplot as plt
import time
import argparse

parser = argparse.ArgumentParser(add_help = False)
parser.add_argument("--problem", type = str, choices = ["traction", "necking", "inclusion"], default = "traction")
parser.add_argument("--clscale", type = float, default = 0.1) ## default should be something quite coarse, so can run quickly
parser.add_argument("--n_iterates", type = int, default = 10)
parser.add_argument("--incomp_strain_elt", type = str, choices = ["HCT", "Regge"], default = "HCT")
parser.add_argument("--incomp_strain_deg", type = int, default = 3)
parser.add_argument("--theta_deg", type = int, default = 0)
parser.add_argument("--theta_space", type = str, choices = ["DG", "CG"], default = "DG")
parser.add_argument("--method", type = str, choices = ["Newton", "bisection"], default = "Newton")
parser.add_argument("--nonlinear_nmax", type = int, default = 10)
parser.add_argument("--scale_displacement", type = float, default = 10) ## factor by which to scale the displacement, just for plotting
parser.add_argument("--print_all_args", type = str, choices = ["True", "False"], default = "True")
## TODO option of whether or not calculate the total work?
args, unknown = parser.parse_known_args()
if (len(unknown) > 2): ## There will always be the 2 unrecognised arguments "-options_left" "0"
    if not unknown[0] == "-options_left":
        ## Unrecognised arguments correspond, in our setting, to typos/misspecified arguments, so we 
        ## want to error.
        ## The alternative, args = parser.parse_args(), also does this, but then doesn't recognise 
        ## -options_left 0 which we want to pass to PETSc.
        raise Exception("You may have passed an unrecognised argument to argparse: possibly " + str(unknown[0]))
white  = "\033[0m"
yellow = "\033[33m"
green = "\033[92m"
blue = "\033[34m"
red = "\033[31m"
orange = "\033[33m"
purple = "\033[35m"
cyan = "\033[0;36m"
if eval(args.print_all_args):
    to_print = ""
    for arg in vars(args):
        #print(arg + " " + cyan + str(getattr(args, arg)) + white)
        to_print += arg + " " + cyan + str(getattr(args, arg)) + white + ", "
    to_print = to_print[:len(to_print) - 2] ## remove final comma
    print(to_print)

start_time = time.time()

#problem = "traction"
#problem = "necking"
problem = args.problem

if problem == "traction":
    meshname = "perforated_plate"
    #clscale = 0.01 ## SO SMALL that the iterations did not begin overnight (27/7/25)
    #clscale = 0.015 ## Killed (8/8/25)
    #clscale = 0.02 ## time: 2h48, or 6h55min, or 7h25min
    #clscale = 0.03 ## time: 
    clscale = 0.04 ## time: 2h8min ## -clscale 0.04 gives 4,644 cells, which is close to the traction example's 4,704
    #clscale = 0.05 ## time: I think this is the largest for which the asymmetry of the mesh does not make the solution asymmetric. But may have been solved by fixing the theta loop.
    #clscale = 0.1 ## time: 12min39, or 19min33
    #clscale = 0.2 ## time: 3min55, or 8min04
    #clscale = 0.3 ## time: 3min41
    #clscale = 0.4 ## time: 1min53

    sigma_Y = 1 ## the yield (according to Amstutz code)
    n_iterates = 20 ## 20 in the paper
    scale_displacement = 10 ## 10 in the paper
elif problem == "necking":
    meshname = "necking_plate"
    #clscale = 0.01 ## time: 
    #clscale = 0.0132 ## time: 3h35 ## -clscale 0.0132 gives 6,744 cells, which is close to the traction necking _refined_ example's 6,750 
    #clscale = 0.018 ## time: 1h58 ## -clscale 0.018 gives 3,802 cells, which is close to the traction necking _coarser_ example's 3,806
    #clscale = 0.02
    #clscale = 0.05
    clscale = 0.1 ## time: 2min08

    sigma_Y = numpy.sqrt(0.9)
    n_iterates = 10 ## 50 in the paper
    scale_displacement = 5 ## 5 in the paper
else: ## inclusion, corresponding to Amstutz's inhomogeneity.edp
    meshname = "inclusion_plate"
    #clscale = 0.018 ## time: 57min20 ## -cscale 0.018 gives 15,377 cells, which is close to the inclusion example's 15,046, IF fine_length is taken as 1.0 in inclusion_plate.geo
    #clscale = 0.04
    #clscale = 0.0454 ## -cscale 0.0454 gives 15,163 cells, which is close to the inclusion example's 15,046, IF fine_length is taken as 0.1 in inclusion_plate.geo
    clscale = 0.1 ## time: 2min31, or 20min24 on laptop
    
    sigma_Y = 1 ## ?
    n_iterates = 10 ## ?
    scale_displacement = 0.1 ## in the paper, 0.1/[Linfty norm of u[0]]

## Leave all of the above (despite the use of argparse), since it's useful to know the runtimes, and what values were used in the paper
clscale = args.clscale
n_iterates = args.n_iterates
scale_displacement = args.scale_displacement

#N_base = 10
#msh = UnitSquareMesh(N_base, N_base) 

print("The " + problem + " problem. Mesh clscale is " + str(clscale))

msh = Mesh("meshes/" + meshname + "_clscale" + str(clscale) + ".msh") ## Don't call it mesh, else will clash with firedrake.mesh
print("Loaded mesh file:", "meshes/" + meshname + "_clscale" + str(clscale) + ".msh")
print("Number of cells: " + str(FunctionSpace(msh, "DG", 0).dim())) ## TODO very weirdly, including this causes an error in the "vertices = " line below, at least with -clscale 0.04 for the traction example
print("Mesh topological dimension:", msh.topological_dimension())
dm = msh.topology_dm
sec = dm.getCoordinateSection()
coords = dm.getCoordinatesLocal()
faces = dm.getStratumIS("exterior_facets", 1).indices
print("All closure sizes:", [dm.vecGetClosure(sec, coords, f).size for f in faces])
x, y = SpatialCoordinate(msh)

## Obtain boundary labels
LEFT = 1
RIGHT = 2
TOP = 3
BOTTOM = 4
HOLE = 5
dm = msh.topology_dm
sec = dm.getCoordinateSection()
coords = dm.getCoordinatesLocal()
dm.removeLabel(dmcommon.FACE_SETS_LABEL)
dm.createLabel(dmcommon.FACE_SETS_LABEL)
faces = dm.getStratumIS("exterior_facets", 1).indices
if problem == "traction" or problem == "inclusion":
    for face in faces:
        #vertices = dm.vecGetClosure(sec, coords, face).reshape(2, 2)
        closure = dm.vecGetClosure(sec, coords, face)
        if len(closure) != 4:
            continue  # skip degenerate boundary entities (points, not edges)
        vertices = closure.reshape(2, 2)

        if numpy.allclose(vertices[:, 0], 0):
            dm.setLabelValue(dmcommon.FACE_SETS_LABEL, face, LEFT)
        elif numpy.allclose(vertices[:, 0], 1.0):
            dm.setLabelValue(dmcommon.FACE_SETS_LABEL, face, RIGHT)
        elif numpy.allclose(vertices[:, 1], 0):
            dm.setLabelValue(dmcommon.FACE_SETS_LABEL, face, BOTTOM)
        elif numpy.allclose(vertices[:, 1], 1.0):
            dm.setLabelValue(dmcommon.FACE_SETS_LABEL, face, TOP)
        else:
            dm.setLabelValue(dmcommon.FACE_SETS_LABEL, face, HOLE)
else: ## necking
    PERFORATIONS = 5
    for face in faces:
        vertices = dm.vecGetClosure(sec, coords, face).reshape(2, 2)
        if numpy.allclose(vertices[:, 0], 0):
            dm.setLabelValue(dmcommon.FACE_SETS_LABEL, face, LEFT)
        elif numpy.allclose(vertices[:, 0], 1.0):
            dm.setLabelValue(dmcommon.FACE_SETS_LABEL, face, RIGHT)
        elif numpy.allclose(vertices[:, 1], 0):
            dm.setLabelValue(dmcommon.FACE_SETS_LABEL, face, BOTTOM)
        elif numpy.allclose(vertices[:, 1], 0.5):
            dm.setLabelValue(dmcommon.FACE_SETS_LABEL, face, TOP)
        else:
            dm.setLabelValue(dmcommon.FACE_SETS_LABEL, face, PERFORATIONS)


## Refresh Firedrake's cached exterior-facet marker list after manual relabelling.
marker_ids = dm.getLabelIdIS(dmcommon.FACE_SETS_LABEL)

if marker_ids is None:
    raise RuntimeError("No exterior boundary markers were created.")

msh.exterior_facets.unique_markers = sorted(
    int(marker) for marker in marker_ids.indices
)

print(
    "Boundary markers recognised by Firedrake:",
    msh.exterior_facets.unique_markers
)


if problem == "inclusion":
    mu_A_inside = 1e-3
    mu_A_outside = 0.3
    mu_D_inside = (1e-5)/2
    mu_D_outside = (1e3)/2
    
    ## I presume this could instead be done with labels
    R = 0.15
    difference = sqrt( (x - 0.5)**2 + (y - 0.5)**2 ) - R
    indicator = conditional(ge(difference, 0), 0.0, 1.0)
    inside_indicator = project(indicator, FunctionSpace(msh, "DG", 0))
    outside_indicator = project(1.0 - indicator, FunctionSpace(msh, "DG", 0))

e1 = Constant([1, 0])
#e1 = as_vector([1, 0])
e2 = Constant([0, 1])
#e2 = as_vector([0, 1])

## Choice of elements, and fix all polynomial degrees
#incomp_strain_elt = "HCT"
#incomp_strain_elt = "Regge"
incomp_strain_elt = args.incomp_strain_elt
## The lowest order for which HCT is defined is 3
#incomp_strain_deg = 3
incomp_strain_deg = args.incomp_strain_deg

if incomp_strain_elt == "HCT":
    ## The paper says P2 for u and p (where HCT is cubic), but could also do disp degree of incomp_strain_deg + 1, due to strain = incomp strain Ei + epsilon(u). Disp deg = strain deg + 1 is also suggested by the Regge complex.
    ## Amstutz says he chose 1 higher degree for E_i because we 'work with' inc E_i and epsilon(u), but these are resp. [the inc of the incompatible strain], and [the compatible strain]
    disp_deg = incomp_strain_deg - 1
    pressure_deg = incomp_strain_deg - 1
else: ## Regge
    disp_deg = incomp_strain_deg + 1
    pressure_deg = incomp_strain_deg + 1 # ? TODO. It was indeed + 1 in my notes/slides

## This choice of degree is since 2 derivs fall on Ei, then 2 derivs fall on epsilon(u) -- ALTHOUGH in fact latter vanishes since inc(epsilon(u)) is always 0.
incE_deg = max(1, incomp_strain_deg - 2, disp_deg - 3)
#theta_deg = 0 
theta_deg = args.theta_deg
## Choice of degree here: |dev(Ei + epsilon(u))|^2 has degree 2*(incomp_strain_deg + disp_deg), and mu_A(theta) is rational (the reciprocal of something affine in theta), so just make large.
sigma_M_squared_deg = 2*(incomp_strain_deg + disp_deg) + 10
## E = Ei + 1 deriv of u
E_deg = max(1, incomp_strain_deg, disp_deg - 1)

## Vector-valued column-wise rot acting on matrix fields
## (should coincide with row-wise rot by symmetry of the matrix)
def Rot(e):
    #return as_vector([curl([e[0, 0], e[0, 1]]), curl([e[1, 0], e[1, 1]])])
    return as_vector([curl(dot(e, e1)), curl(dot(e, e2))])
    #return [curl(dot(e, e1)), curl(dot(e, e2))]

def inc(e):
    #return curl(Rot(e))
    return Dx(Dx(e[1, 1], 0), 0) + Dx(Dx(e[0, 0], 1), 1) - 2*Dx(Dx(e[0, 1], 1), 0)

def L2norm(fn):
    return sqrt(assemble(inner(fn, fn)*dx))

## Linearised strain tensor
def epsilon(u):
    return 0.5*(grad(u) + grad(u).T)

## Jump operator for the vector-valued quantity A 
## dotted with N, with N = n or t.
def Jump(A, N):
    #return 2*avg(outer(A, N))
    return 2*avg(inner(A, N))

## Total strain = incompatible strain + epsilon(u) (Beltrami decomp)
def fullstrain(Ei, u):
    return Ei + epsilon(u)

n = FacetNormal(msh)
## Unit tangent defined by ACW rotation of the normal
#t = as_vector([- dot(n, e2), dot(n, e1)])
## Adapted from Guosheng's code hdg2d.py, doing not tangent but tangential COMPONENT
#def tang_component(vec):
#    return vec - dot(vec, n)*n

h = CellDiameter(msh) 

## Note to self: dS is interior edges
d_all_edges = ds + dS

def discrete_airy_D_inc(Ei, Fi, theta): # D depends on theta
    ## Discrete weak form of airy(D inc E), which in 3D is inc(D inc E)
    if incomp_strain_elt == "HCT":
        return exact_airy_D_inc(Ei, Fi, theta)
    else: ## Regge
        ## TODO see how the D coefficient should be incorporated into the rest of the DG version of airy_D_inc
        airy_D_inc = (
#            inner(inc(Ei), inc(Fi))*dx
            inner((lam_D + 2*mu_D(theta))*inc(Ei), inc(Fi))*dx ## the coeff here is just 1/theta
            ## TODO Rest of discrete form for Regge to be settled by FA.
        )
        return airy_D_inc

def exact_airy_D_inc(Ei, Fi, theta): ## D depends on theta
    a = (
        #inner(inc(Ei), inc(Fi))*dx
        inner((lam_D + 2*mu_D(theta))*inc(Ei), inc(Fi))*dx ## the coeff here is just 1/theta
    )
    return a

I = Identity(2)

## Bulk modulus; always assumed constant (near top of p18)
kappa_A = Constant(83)
## Shear modulus
mu_0 = 38.46 ## This one cannot be a Constant() due to being used in numpy etc.
k = 1e4 ## Similarly cannot be Constant()
#print("k/mu_0 = " + str(k/mu_0)) ## since this is a term in h'(r)
## Eq. (39)
gamma = (sigma_Y**2)/(4*k) ## Similarly cannot be Constant()

theta_min = 1e-3

## top of p18
lam_D = 0.0

## Eq. (6)
def lam_A(theta):
    return kappa_A - (2/3)*mu_A(theta)

def tr(Ei):
    return Ei[0, 0] + Ei[1, 1]

## Eq. (38)
def tilde_mu(theta):
    return k/theta

## Eq. (37)
## This is decreasing in theta.
def mu_A(theta):
    if problem == "inclusion":
        return mu_A_inside*inside_indicator + mu_A_outside*outside_indicator
    else:
        mu_A_inv = (1/mu_0) + (1/tilde_mu(theta))
        return 1.0/mu_A_inv

## NEW for HCT formulation
def as_sym_matrix(Ei1, Ei2, Ei3):
    return as_tensor([[Ei1, Ei2], [Ei2, Ei3]])

## Eq. (5)
def A(Ei, theta):
    return 2*mu_A(theta)*Ei + lam_A(theta)*tr(Ei)*I

## Above Eq. (34)
def mu_D(theta):
    if problem == "inclusion":
        return mu_D_inside*inside_indicator + mu_D_outside*outside_indicator
    else:
        return 1.0/(2*theta)

## Unused
## Few lines above Eq. (34), inverting the above method
#def Theta(mu_D):
#    return 1.0/(2*mu_D)

## Eq. (5)
def D(Ei, theta):
    return 2*mu_D(theta)*Ei + lam_D*tr(Ei)*I

## Deviator
def dev(sig):
    #return sig - 0.5*tr(sig)*I
    ## 1/3 here since that's used for calculation of the von Mises stress in AS's code;
    ## I guess it also means this is all technically 3D.
    return sig - (1/3)*tr(sig)*I 

## Total work functional, top of p18
## This is a squared norm of the full strain E, which is equivalent
## to the H(inc) norm (with constants depending on theta).
def W(Ei, u, theta):
    E = fullstrain(Ei, u)
    integral = 0.5*(
        inner(kappa_A*tr(E), tr(E))*dx
        + inner(2*mu_A(theta)*dev(E), dev(E))*dx
        ## TODO come up with Regge/DG equivalent of this. Strictly, ill-defined since you can't square the Dirac delta (even as a distribution).
        ## But this is still well-defined numerically since Firedrake will just calculate using the piecewise inc of Regge.
        + inner(2*mu_D(theta)*inc(E), inc(E))*dx
    )
    return assemble(integral)

## Overall energy fnl.
## Eq. (36)
def energy(Ei, u, theta):
## TODO what's the meaning of "Of course, K has to be constructed in the form equation (54)." ?
    K_integral = (
    ## This is the same as -int(K:E), by the decomposition of E below (40), and comparing (26)
    ## with the 2nd identity on p5. This is also used just above (41).
        - inner(f, u)*dx
        - inner(g_from_left, u)*ds(LEFT)
        - inner(g_from_right, u)*ds(RIGHT)
    )
    if problem == "inclusion":
        K_integral += (
            - inner(g_from_bottom, u)*ds(BOTTOM)
            - inner(g_from_top, u)*ds(TOP)
        )

    theta_integral = (
        inner(gamma*theta, Constant(1))*dx ## phi(theta) = gamma*theta from Eq. (38)
    )
    return W(Ei, u, theta) + assemble(K_integral + theta_integral)

## Save a plot of all variables of interest
def save_plot(Ei, u, p, theta, time):
    ## Plot the solutions
    E = fullstrain(Ei, u)
    E = project(E, FullStrainTensorSpace)
    incE = inc(E) ## Could do with only Ei here, since inc(epsilon(u)) == 0.
    incE = project(incE, FunctionSpace(msh, FiniteElement("DG", msh.ufl_cell(), incE_deg)))

    ## above eq. (39). Note this is of the full strain E, not just the incomp strain Ei
    sigma_M_squared = 4*mu_A(theta)*mu_A(theta)*inner(dev(E), dev(E))
    sigma_M_squared = project(sigma_M_squared, FunctionSpace(msh, FiniteElement("DG", msh.ufl_cell(), sigma_M_squared_deg)))
    ## Losslessly project/interpolate the incomp_strain into matrix Lagrange space of same degree, so can name the variable. This is also needed in order to plot it. It is lossless if project into DG.

    ## Plot the sign of the von Mises bound (39). This quantity has the same polynomial degree as sigma_M_squared.
    yield_squared_minus_sigma_M_squared = project(Constant(sigma_Y**2) - sigma_M_squared, FunctionSpace(msh, FiniteElement("DG", msh.ufl_cell(), sigma_M_squared_deg)))

    u = project(scale_displacement*u, VectorFunctionSpace(msh, "CG", disp_deg))

    Ei = project(Ei, IncompStrainTensorSpace)
    mu_A_theta = assemble(interpolate(mu_A(theta), thetaspace)) ## lossy, although a rational function of DG0 would still be DG0.
    Ei.rename("incompatible strain Ei")
    u.rename("displacement")
    p.rename("pressure")
    theta.rename("theta")
    mu_A_theta.rename("mu_A(theta)")
    E.rename("full strain")
    incE.rename("inc of full strain")
    sigma_M_squared.rename("von Mises stress squared")
    yield_squared_minus_sigma_M_squared.rename("yield stress squared minus von Mises stress squared")

    outfile.write(Ei, u, p, theta, mu_A_theta, E, incE, sigma_M_squared, yield_squared_minus_sigma_M_squared, time = time) ## This command is independent of HCT vs Regge

## Define the methods to be used in bisection/Newton for the theta update, so they aren't 
## redefined for each DOF.
## TODO they had been defined WITHIN update_theta_dof since they depended on r_star, which 
## depends on the DOF, but could remove this dependence, at the cost of not being able to ALSO 
## print out r_star upon failure. The cost of passing an extra parameter is presumably less 
## than the cost of redefining these 2 or 3 methods with each DOF.
## I REALISE NOW they also depend on a, b, hence do depend on the dof.
## So TODO define these instead via lambdas, within the update_theta_dof method?
def h(r, a, b): # on p28
    ## No need for this strictly, but may as well test its value
    return 2*sqrt(b*r) - (sqrt(a*mu_0) - sqrt(k*(1 - r)/mu_0))**2

def h_prime(r, a, b, r_star): ## easily checked by hand/Wolfram Alpha
    ## This causes division by zero for r near 0 or 1.
    try:
        #hprime = sqrt(b/r) - sqrt(a*k/(1 - r)) + k/mu_0
        return sqrt(b/r) - sqrt(a*k/(1 - r)) + k/mu_0
        #hprime = sqrt(b/r) - (sqrt(a*mu_0) - sqrt(k*(1 - r)/mu_0))*sqrt(k/(mu_0*(1 - r))) # alternative form
        #print("h' worked for r* = " + str(r_star) + " and r = " + str(r))
    except:
        print("h' failed for r* = "+ str(r_star) + " and r = " + str(r))
        #hprime = 1e10
        return 1e10
    #return hprime

def h_prime_prime(r, a, b, r_star): ## easily checked by hand/Wolfram Alpha
    ## This causes division by zero for r near 0 or 1.
    try:
        #hprimeprime = - (1/2)*sqrt(b/(r**3)) - (1/2)*sqrt(k*a/((1 - r)**3))
        return - (1/2)*sqrt(b/(r**3)) - (1/2)*sqrt(k*a/((1 - r)**3))
        #print("h'' worked for r* = " + str(r_star) + " and r = " + str(r))
    except:
        print("h'' failed for r* = "+ str(r_star) + " and r = " + str(r))
        #hprimeprime = 1e10
        return 1e10
    #return hprimeprime

## The Newton loop (or bisection) at each DOF of theta
bisection_Nmax = 50 ## 50 was used in Amstutz' MATLAB code
bisection_iterate_tol = 1e-2 ## Was previously taking both these bisection tolerances as 1e-12, then bisection_EPS as 1e-13
bisection_residual_tol = 1e-9
bisection_EPS = 1e-13 ## 1e-8 in Amstutz's MATLAB bisection code, but I've found r can become within 1e-8 of 1. NOTE it does not make sense for this to be larger than bisection_iterate_tol!
bisection_left_EPS = 1e-8
bisection_right_EPS = 1e-8

newton_Nmax = 10 ## 10 was used in Amstutz' FreeFEM code (for the traction and necking examples), [TODO? with inclusion], 20 in the MATLAB code
newton_iterate_tol = 1e-7 ## Previously was taking both these Newton tolerances as 1e-12, then newton_EPS as 1e-13
newton_residual_tol = 1e-7
newton_EPS = 1e-8 ## As above, it does not make sense for this to be larger than newton_iterate_tol!
############ Once again, it's useful to keep all the above for the comments (e.g. what values were used in the paper etc.)
bisection_Nmax = args.nonlinear_nmax
newton_Nmax = args.nonlinear_nmax
############

## tol just for checking whether a, b are zero
tol = 1e-12
def update_theta_dof(a, b, method):
    ## Get rid of the trivial cases that either a or b is zero
    if (abs(a) < tol) and (abs(b) < tol):
        return 0
    elif (abs(a) < tol) and (b > 0): ## Writing it this way just to make it more readable.
        ## Of course, this technically can overlap with the previous if condition.
        return sqrt(b)
    elif (a > 0) and (abs(b) < tol):
        ## Note: if a is still essentially 0, then this returns -(k/mu_0) which is negative,
        ## so is ignored in the theta update due to the max(0, .)
        return sqrt(a*k) - (k/mu_0)
    ## If got to THIS stage, then h' should admit a root in (0, 1)
    ## although possibly not in [r*, 1).

    ## Regularise b
    #b += bisection_left_EPS

    ## Here a, b are as in the abstract form from the paper appendix.
    
    #r_star = -1 ## This is set within the choice of method (bisection or Newton).
    ## In bisection, to mimic Amstutz's MATLAB code flowrule.m, while in Newton to mimic Amstutz's FreeFEM code traction.edp
    ## -> NO changing to check whether h'(r*) <= 0, since that just means within [r*, 1), h' has no roots and r* is the maximiser of h. Use the generic value tol as "zero" to regularise.
    r_star = max(tol, 1 - (a*(mu_0**2)/k))

    if h_prime(r_star, a, b, r_star) <= 0:
        r = r_star
    elif method == "bisection":
        ## Bisection method to find the root of h_prime, to maximise h
        ## Slight regularisation in order to avoid division by zero when initially evaluating h'(0) or h'(1), i.e. EPS should be zero
        #r_star = max(bisection_EPS, 1 - (a*(mu_0**2)/k) - bisection_EPS) # This 2nd instance of EPS was not there in Amstutz's MATLAB bisection code; it's to deal with the case a = 0, giving r* = 1
        r_star = max(bisection_left_EPS, 1 - (a*(mu_0**2)/k)) ## from Amstutz's MATLAB bisection code
        
        ## Define endpoints of the interval
        left = r_star
        #print("r_star = " + str(r_star))
        right = 1 - bisection_right_EPS
        
        ## Check whether the sign condition for starting bisection is fulfilled. This should be true, 
        ## now that a, b > 0 and h'(r*) > 0.
        #sgn = numpy.sign(h_prime(left, a, b, r_star)*h_prime(right, a, b, r_star))
        ##print(sgn)
        #try:
        #    assert sgn < 0
        #except:
        #    print("Sign assumption for bisection failed, with left = " + str(left) + ", b = " + str(b) + ", h'(left) = " + str(h_prime(left, a, b, r_star)) + ", right = " + str(right) + ", h'(right) = " + str(h_prime(right, a, b, r_star)) + ", sqrt(ak/(1 - left)) = " + str(sqrt(a*k/(1 - left))) + ", a = " + str(a) + ", sqrt(b/left) = " + str(sqrt(b/left)))
        #    #raise Exception(
        #    ## Plot h' failing to go through 0
        #    x = numpy.linspace(left, right, 1000)
        #    plt.figure()
        #    plt.plot(x, [h_prime(t, a, b, r_star) for t in x])
        #    plt.show()

        j = 0
        while True:
        #for j in range(Nmax):
            m = (left + right)/2
            #print("h(r) at this point = " + str(h(m)))
            #print("h'(m) = " + str(h_prime(m, a, b, r_star)) + ", h'(left) " + str(h_prime(left, a, b, r_star)) + ", h'(right) = " + str(h_prime(right, a, b, r_star)))
            #print("m = " + str(m))
            j = j + 1
            H_PRIME = h_prime(m, a, b, r_star)
            #if (abs((m - left)/2) < bisection_iterate_tol) or (abs(h_prime(m, a, b, r_star)) < bisection_residual_tol):
            #if j >= bisection_Nmax: 
            if (j >= bisection_Nmax) or (abs(H_PRIME) < bisection_residual_tol): ## Had to add this 2nd criterion to what Amtutz's MATLAB code had, since h'(m) was literally becoming 0.0
                ## Could do a for (instead of while) loop, but then have to re-code this initial if statement
                ## Alternative stopping criterion: relative error abs((m - left)/m) < tol
                break
            elif numpy.sign(h_prime(right, a, b, r_star)) == numpy.sign(h_prime(m, a, b, r_star)):
            ## NOTE THIS SIGN CONDITION SHOULD BE CHECKED FIRST since h' is decreasing
            #elif H_PRIME < 0: # Similarly: used in Amstutz's bisection code
                right = m
            elif numpy.sign(h_prime(left, a, b, r_star)) == numpy.sign(h_prime(m, a, b, r_star)):
            #elif H_PRIME > 0: # This is used in Amstutz's bisection code
                left = m
            else:
                ## NOTE none of these if statements pass if h' takes complex values, which is true for a or b being < 0. So
                ## take their positive part before passing them to this method.
                raise Exception("bisection failed, with left = " + str(left) + ", m = " + str(m) + ", right = " + str(right) + ", h'(m) = " + str(H_PRIME) + ", h'(left) " + str(h_prime(left, a, b, r_star)) + ", h'(right) = " + str(h_prime(right, a, b, r_star)))
            #print(j)
        r = m
    else: # elif method == "Newton":
        #r_star = max(newton_EPS, 1 - (a*(mu_0**2)/k)) ## TODO could add newton_left_EPS and _right_EPS
        r_star = max(0.0, 1 - (a*(mu_0**2)/k)) # To try to best mimic the FreeFEM code (where the max was taken when defining r, there called q): there was no regularisation of r* (there called q*) or the initial r
        ## There WAS eps = 1e-8 in the first argument here in Amstutz's MATLAB code flowrule.m (though the point of that was to demonstrate bisection)

        r = (r_star + 1)/2 ## The starting value used in traction.edp

        j = 0
        while True:
        #for j in range(newton_Nmax):
            r_prev = r
            r = r_prev - h_prime(r_prev, a, b, r_star)/h_prime_prime(r_prev, a, b, r_star)
            ## Project onto [r*, 1] if go outside it:
            #r = max(r_star, min(1 - newton_EPS, r))
            r = max(max(r_star, newton_EPS), min(1 - newton_EPS, r)) ## To best mimic the FreeFEM code
            j = j + 1
            #if (abs(r - r_prev) < newton_iterate_tol) or (abs(h_prime(r, a, b, r_star)) < newton_residual_tol):
            if j >= newton_Nmax:
                break
        
    #print("r = " + str(r))
    try:
        ## Eq. (55)
        theta_dof = max(0, sqrt(a*k/(1 - r)) - (k/mu_0))
    except:
        print("Updating theta by (55) failed, with r = " + str(r) + "and ak/(1 - r) = " + str(a*k/(1 - r)) + ", a = " + str(a) + ", k = " + str(k))
        theta_dof = 1e10

    ## theta appeared to be getting large in a few random cells, although that should happen for some of them
    #if theta_dof > 100:
    #    print("large theta_dof = " + str(theta_dof) + ", r = " + str(r) + ", r* = " + str(r_star) + ", h'(r) = " + str(h_prime(r, a, b, r_star)) + ", h'(r*) = " + str(h_prime(r_star, a, b, r_star)) + ", a = " + str(a) + ", b = " + str(b))

    return theta_dof

outfile = VTKFile(incomp_strain_elt + "_" + problem + "/output" + str(clscale) + ".pvd")

lu = {
        "snes_type": "ksponly",
        "snes_linesearch_type": "basic",
        "snes_max_it": 100,
        "snes_monitor": None,
        "snes_converged_reason": None,
        "ksp_type": "preonly",
        "pc_type": "lu",
        "pc_factor_mat_solver_type": "mumps",
        "mat_mumps_icntl_14": "1000",
        "mat_type": "aij",
        "snes_stol": 0.0,
}
direct = {
        "ksp_type": "preonly",
        "pc_type": "lu",
        #"snes_monitor": None, 
}
from_scratch = {
        "ksp_type": "preonly", ## NOTE THIS FORCES IT TO DO A LINEAR SOLVE since the Ei-u-p system is linear. I think to be precise, this forces a direct solve when you have solve(bilinear == linear, ...), in which case it's a linear solve in any case.
        "pc_type": "lu",
        "ksp_monitor": None,
        "snes_monitor": None, 
        #"snes_max_it": 100, ## Some had been going beyond 10
        #"snes_rtol": 1e-8, ## 1e-8 is default
}

## RHS forcing, and boundary terms
f = Constant([0.0, 0.0])
g_from_left = Constant([-1.0, 0.0])
#g_from_left = Constant([-0.5, 0.0])
g_from_right = Constant([1.0, 0.0])
#g_from_right = Constant([0.5, 0.0])
if problem == "inclusion":
    g_from_bottom = Constant([0.0, -1.0])
    g_from_top = Constant([0.0, 1.0])

## Stabilisation parameters
eps_u = Constant(1e-6)
eps_p = Constant(1e-6)

## Strain-displacement-pressure mixed space
if incomp_strain_elt == "HCT":
    IncompStrainComponentElement = FiniteElement("HCT", triangle, incomp_strain_deg)
    ## Note here the first of the 3 spaces is itself 3 functions, corresponding to the 3 components of the symmetric strain in 2D
    triplet_element = MixedElement([IncompStrainComponentElement, IncompStrainComponentElement, IncompStrainComponentElement, VectorElement("CG", msh.ufl_cell(), disp_deg), VectorElement("CG", msh.ufl_cell(), pressure_deg)])
else: ## Regge elt
    IncompStrainElement = FiniteElement("Regge", triangle, incomp_strain_deg)
    triplet_element = MixedElement([IncompStrainElement, VectorElement("CG", msh.ufl_cell(), disp_deg), VectorElement("CG", msh.ufl_cell(), pressure_deg)])

triplet_space = FunctionSpace(msh, triplet_element)
Eiup = Function(triplet_space)
## Compared with Section 9:
## E_i there corresponds to Ei here, \hat{u} to v, \hat{E}_i to Fi, \hat{p} to q
if incomp_strain_elt == "HCT":
    (Fi1, Fi2, Fi3, v, q) = TestFunctions(triplet_space)
    Fi = as_sym_matrix(Fi1, Fi2, Fi3)
else: ## Regge
    (Fi, v, q) = TestFunctions(triplet_space)

## Want a FE space for which each DOF has an associated physical point (see below).
## Brubeck confirmed on fd slack this IS true by default for DG.
## I think best to stick with DG0 for now, since that was working in Amstutz's code.
#if theta_deg == 0:
#    # in Amstutz's code this was DG0.
#    # I've found it has to be DG0 rather than e.g. CG2, then the linear solves fail not at all
#    # (or perhaps on a later iteration).
#    thetaspace = FunctionSpace(msh, "DG", theta_deg)
#else:
#    #thetaspace = FunctionSpace(msh, "CG", theta_deg)
#    thetaspace = FunctionSpace(msh, "DG", theta_deg)
thetaspace = FunctionSpace(msh, args.theta_space, args.theta_deg)

## The starting value used in the FreeFEM code was 1e-3
theta_init = Constant(1e-3)
#theta_init = Constant(1e1)
#theta_init = 1000*sin(pi*x)*sin(pi*y) # this is symmetric about the centre of the unit square [0, 1]^2

theta_init = assemble(interpolate(theta_init, thetaspace))
theta = Function(thetaspace) ## Note: theta has/needs no associated test functions.
theta.assign(theta_init)

## NOTE strictly this space, and projection onto it, is only needed in the HCT case, but that's ok
#IncompStrainTensorSpace = TensorFunctionSpace(msh, "CG", incomp_strain_deg, shape = (2, 2)) ## This space should contain the sym-matrix-valued HCT space -> NOPE CG is not H2-conforming.
IncompStrainTensorSpace = TensorFunctionSpace(msh, "DG", incomp_strain_deg, shape = (2, 2)) ## Changed this to DG so that projection onto it is lossless.
#IncompStrainTensorSpace = TensorFunctionSpace(msh, "HCT", incomp_strain_deg, shape = (2, 2)) ## Could have used this space and removed as_sym_matrix - but this errors when try to outfile.write

FullStrainTensorSpace = TensorFunctionSpace(msh, "DG", E_deg, shape = (2, 2)) ## DG so that projection onto it is lossless.

## TODO could put this into the plotting method
## Plot the initial conditions (but only theta has an initial condition)
if incomp_strain_elt == "HCT":
    (Ei1, Ei2, Ei3, u, p) = Eiup.subfunctions
    Ei = as_sym_matrix(Ei1, Ei2, Ei3)
else: ## Regge
    (Ei, u, p) = Eiup.subfunctions

## All of these except theta will be identically zero, but I have to plot them all since
## they are to be saved in the same .vtu files
save_plot(Ei, u, p, theta, 0.0)

## TODO try the formulation without p, Eq. (44) -- but then have to discretise the space 
## Y in Eq. (19)
## TODO INSTEAD - try (42), but passing the nullspace of RM into the solver for both u and p.
## Stopping criterion could instead be the Cauchy property (i.e. update between iterations becomes small).

## Energy (and total work) plotting: want initial energy, then 0.5 later after Ei-u-p update, 
## then after 1 overall iteration (having updated theta), etc.
energies = [0.0 for j in range(1 + 2*n_iterates)]
total_works = [0.0 for j in range(1 + 2*n_iterates)]
energy_xticks = [j/2 for j in range(2*n_iterates)] + [n_iterates]
## TODO could put this into a method which calculates and prints the energy and work; its output
## can be saved to the arrays
initial_energy = energy(Ei, u, theta)
initial_work = W(Ei, u, theta)
print(green + "Initial energy is " + str(initial_energy) + white)
print(blue + "Initial total work is " + str(initial_work) + white)
energies[0] = initial_energy
total_works[0] = initial_work

## This is here outside the loop, so then have the outputs of .subfunctions as new_Ei1, new_Ei2 etc.
#(Ei1, Ei2, Ei3, u, p) = split(Eiup) 
#Ei = as_sym_matrix(Ei1, Ei2, Ei3)
## use of TrialFunction rather than Function is necessary for when doing solve(bilinear == linear, ...)
#Eiup = TrialFunction(triplet_space)
## I think the below is equivalent to split(TrialFunction(triplet_space))
if incomp_strain_elt == "HCT":
    (Ei1, Ei2, Ei3, u, p) = TrialFunctions(triplet_space)
    Ei = as_sym_matrix(Ei1, Ei2, Ei3)
else: ## Regge
    (Ei, u, p) = TrialFunctions(triplet_space)

coordinate_space = VectorFunctionSpace(msh, "DG", theta_deg)
X = SpatialCoordinate(msh)
dof_coordinates_function = Function(coordinate_space).interpolate(X)
dof_coordinates = dof_coordinates_function.dat.data_ro

#method = "bisection"
#method = "Newton"
method = args.method

for increment in range(n_iterates):
    print("\nIteration number " + str(increment + 1))
    ## Update Ei, u, p. Amstutz: we can make a rough initial guess for theta but not for E, so 
    ## it makes sense to minimise first wrt E.

    ## TODO if make this linear solve a numerically nonlinear solve using solve(): interpolate 
    ## new_Ei etc. into the trial functions here

    ## Eq. (42)
    #Form = (
    bilin = (
    ## NOTE inc need not commute with the application of D! So cannot take inc of D(Ei) here.
        #airy_inc(D(Ei, theta), Fi)
        #inner(D(inc(Ei), theta), inc(Fi))*dx
    ## NOPE it takes a specific form in 2D; see p23
        #inner((lam_D + 2*mu_D(theta))*inc(Ei), inc(Fi))*dx ## the coeff here is just 1/theta
        discrete_airy_D_inc(Ei, Fi, theta)
        + inner(A(Ei, theta), Fi)*dx
        + inner(A(epsilon(u), theta), Fi)*dx
        + inner(epsilon(p), Fi)*dx
        + inner(A(Ei, theta), epsilon(v))*dx
        + inner(A(epsilon(u), theta), epsilon(v))*dx
        + inner(eps_u*u, v)*dx ## stabilisation
        + inner(Ei, epsilon(q))*dx
        - inner(eps_p*p, q)*dx ## stabilisation
    )
    rhs = ( ## These should have minus sign out front if do solve(F == 0, ...)
        + inner(f, v)*dx
        + inner(g_from_left, v)*ds(LEFT)
        + inner(g_from_right, v)*ds(RIGHT)
    )
    if problem == "inclusion":
        rhs += (
            + inner(g_from_bottom, v)*ds(BOTTOM)
            + inner(g_from_top, v)*ds(TOP)
        )

    print("Solving incomp_strain-disp-pressure problem:")
    #solve(Form == 0, Eiup, solver_parameters = from_scratch)
    #solve(Form == 0, Eiup, solver_parameters = direct)
    #solve(bilin == rhs, Eiup, solver_parameters = from_scratch)
    solve(bilin == rhs, Eiup, solver_parameters = direct)

    ## The ones renamed 'new' are from the outputs of solve()
    if incomp_strain_elt == "HCT":
        (new_Ei1, new_Ei2, new_Ei3, new_u, new_p) = Eiup.subfunctions
        new_Ei = as_sym_matrix(new_Ei1, new_Ei2, new_Ei3)
    else: ## Regge
        (new_Ei, new_u, new_p) = Eiup.subfunctions

    ## Plot intermediate energy and total work
    intermediate_energy = energy(new_Ei, new_u, theta)
    intermediate_total_work = W(new_Ei, new_u, theta)
    print(green + "New overall energy (after updating Ei-u-p) is " + str(intermediate_energy) + white)
    print(blue + "New total work (after updating Ei-u-p) is " + str(intermediate_total_work) + white)
    energies[1 + 2*increment] = intermediate_energy
    total_works[1 + 2*increment] = intermediate_total_work

    ## Plot all variables at this intermediate stage
    save_plot(new_Ei, new_u, new_p, theta, time = increment + 0.5)
    
    ## Update theta
    ## Need this quantity in the form of a Function. Lossless projection onto scalar DG
    ## Note this is the dev squared and inc squared of the full strain E, not just the incomp strain Ei
    E = fullstrain(new_Ei, new_u)
    dev_E_squared = project(inner(dev(E), dev(E)), FunctionSpace(msh, "DG", 2*E_deg))
    ## Similarly for |inc(E)|^2, but could do with only Ei here, since inc(epsilon(u)) == 0.
    inc_E_squared = project(inner(inc(E), inc(E)), FunctionSpace(msh, "DG", 2*incE_deg))

    ## NOTE since I'm changing the DOFs of the theta Function directly, there's no need for a new_theta (as was done for the other variables Ei, u, p).

    print("Solving the theta problem by the " + method + " method at each of its DOFs.")
    ## Loop over the global DOF indices. Note that since I use the values
    ## of |dev(E)|^2 etc. "at" each DOF of theta, theta needs to be in a Lagrange space
    ## (or more generally any in which each DOF has an associated physical point, rather
    ## than being e.g. an integral over an edge, since I need somewhere to evaluate
    ## |dev(E)|^2 etc.). One could do a theta space not satisfying this, but then would need
    ## to do Newton for the entire function theta, not just at each DOF.
    for j in range(len(dof_coordinates)):
        pt = dof_coordinates[j]
        #print("DOF number " + str(j) + " out of " + str(len(dof_coordinates)))
        #print(pt)
        ## Comparing Appendix I and Eq. (36) to infer what a, b are.
        ## Taking the positive part of it just to be absolutely safe, since found it was slightly 
        ## negative when taking square root of it while updating theta
        ## Amstutz advised to not add an EPS to a.
        a = max(0.0, dev_E_squared.at(tuple(pt))/gamma) # Note this incorporates 2*(1/2).
        #print("a = " + str(a))
        # Taking positive part of it, as above, to avoid sqrt of a negative number.
        b = max(0.0, inc_E_squared.at(tuple(pt))/(2*gamma))
        #print("b = " + str(b))

        ## Careful here: a near 0 gives r* (and hence r) near 1

        ## This step takes a while, so comment out BOTH steps below if want to see only 1 
        ## linear solve without updating theta.
        new_theta_dof = update_theta_dof(a, b, method)
        theta.dat.data[j] = max(theta_min, new_theta_dof)
        #print("new_theta_dof = " + str(new_theta_dof) + ", theta = " + str(theta.at(pt))) ## the latter should evaluate to theta_min, except when new_theta_dof is larger, then they should coincide.

    new_energy = energy(new_Ei, new_u, theta)
    new_total_work = W(new_Ei, new_u, theta)
    ## Printing this out so I can read it, in lieu of having a massive try statement around the for loop 
    ## in case of divergence of one of these overall iterations. Same for the total work.
    print(green + "New overall energy (after updating theta) is " + str(new_energy) + white)
    print(blue + "New total work (after updating theta) is " + str(new_total_work) + white)
    energies[2 + 2*increment] = new_energy
    total_works[2 + 2*increment] = new_total_work

    ## For the inclusion example, have to do this rescaling after the plotting of the ICs, because
    ## the ICs have zero displacement, which would cause division by zero here.
    if problem == "inclusion":
        ## TODO is it right that this should be an linfty norm just in space?
        u_x = project(new_u[0], FunctionSpace(msh, "DG", disp_deg))
        #linfty_norm = (u_x.dat.vec_ro).max()[1] ## errors
        linfty_norm = max(u_x.vector().array())
        #print(linfty_norm)
        scale_displacement = 0.1/linfty_norm
        scale_displacement = float(scale_displacement) ## This is so that scale_displacement*u can be projected before plotting

    save_plot(new_Ei, new_u, new_p, theta, increment + 1) ## increment + 1 since it starts at 0

def convert(seconds): ## from https://www.geeksforgeeks.org/python/python-program-to-convert-seconds-into-hours-minutes-and-seconds/
    mn, s = divmod(seconds, 60)
    hr, mn = divmod(mn, 60)
    return "%d:%02d:%02d" % (hr, mn, s)
print("Elapsed time: " + str(convert(time.time() - start_time)))

## Plot the overall energy eq. (36), which should monotonically decrease with increments.
## The total work should be close to the energy (compare energy (36), with total work at 
## top of p18)
plt.figure()
#plt.plot(range(n_iterates), energies, marker = "o")
plt.plot(energy_xticks, energies, marker = "o", label = "energy (eq. (36))")

## Confusing to have this extra printout and graph, and quite similar to the total energy 
## anyway, so commenting out for now (as for the printouts above)
plt.plot(energy_xticks, total_works, marker = "*", label = "total work (top of p18)")

plt.xlabel("increment")
#plt.ylabel("Overall energy (eq. (36))")
## Force x-axis labels to be integers
axes = plt.gca()
axes.set_xticks(energy_xticks)

## To avoid xticks overlapping, although in practice you can stop that by expanding the window (this doesn't help the saved PDF, though)
plt.tight_layout() 
plt.tick_params(axis = "x", labelsize = 8)
plt.xticks(rotation = 45, ha = "right") ## Rotate x-axis labels by 45 degrees

plt.legend()
plt.savefig(incomp_strain_elt + "_" + problem + "/energy.pdf")
plt.show()

