# -*- coding: utf-8 -*-
"""

Governing equations: Double porosity/permeability
Formulation: Stabilized mixed CG formulation 
Problem: 2D convergence - Square domain 
Run as:
  python <>.py <rtol> <show_monitor> <solver_type> <amg_type> <nx>

  where:
  <rtol> = relative solver tolerance
  <show_monitor> = Turn off KSP residual monitor (0), turn on (1)
  <solver_type> = naive (1), split by fields ("fields"), split by scales ("scales")
  <amg_type> = algebraic multigrid type, either gamg, hypre, or ml
  <nx> = number of cells in each spatial direction
  <element> = T3/Q4 

  I recommend starting off with:
  <rtol> = 1e-5
  <show_monitor> = 1
  <solver_type> = fields or scales
  <amg_type> = hypre
  <element> = T3

  Play around with these as necessary. Few notes:
  - Turn off monitor to improve timing metrics
  - Penaties eta_u and eta_p affect performances
  - Solver count will fluctuate with number of MPI processes
  - If solving 2D problems, comment out the hypre boomeramg parameters
  - Preliminary runs suggest gamg is the worst amg_type
  - Preliminary runs suggest that splitting by scales is the best
"""
from firedrake import *
import numpy as np
from pyop2.profiling import timed_region
from mpi4py import MPI
import math,sys,time
rank = MPI.COMM_WORLD.Get_rank()
size = MPI.COMM_WORLD.Get_size()
rtol = float(sys.argv[1])
show_monitor = int(sys.argv[2])
solvestrategy = sys.argv[3]
amg_type = sys.argv[4]
nx = int(sys.argv[5])
element = sys.argv[6]

#=====================================;
#  Create mesh and identify boundary  ;
#=====================================;
if element == 'T3':
	mesh = UnitSquareMesh(nx,nx)
else:
	mesh = RectangleMesh(nx,nx,1,1,quadrilateral = True)
#====================================================;
#  Define function spaces and mixed (product) space  ;
#====================================================;
vSpace = VectorFunctionSpace(mesh,"CG",1)
pSpace = FunctionSpace(mesh,"CG",1)
wSpace = MixedFunctionSpace([vSpace,pSpace,vSpace,pSpace])

#===================================;
#  Define trial and test functions  ;
#===================================;
(v1,p1,v2,p2) = TrialFunctions(wSpace)
(w1,q1,w2,q2) = TestFunctions(wSpace)

#======================;
#  Define body forces  ;
#======================;
rhob1, rhob2 = Constant((0.0,0.0)), Constant((0.0,0.0))

#============================;
#  Define medium properties  ;
#============================;
mu = Constant(1.0)
beta = Constant(1.0)                     
fact = 1.0  #  fact = beta / mu 
k1, k2 = Constant(1.0), Constant(0.1)
alpha1, alpha2 = Constant(mu/k1), Constant(mu/k2)

eta = Constant(3.316625)

#----------------;
#  Pressure BCs  ;
#----------------;
x, y = SpatialCoordinate(mesh)
p1_left = Function(pSpace).interpolate((1/pi)*sin(pi*y) - exp(3.316625*y))
p1_right = Function(pSpace).interpolate((1/pi)*exp(pi)*sin(pi*y) - exp(3.316625*y))
p1_bottom = Constant(-1.0)
p1_top = Constant(-27.5671484)
p2_left = Function(pSpace).interpolate((1/pi)*sin(pi*y) + 10.0 *exp(3.316625*y))
p2_right = Function(pSpace).interpolate((1/pi)*exp(pi)*sin(pi*y) + 10.0*exp(3.316625*y))
p2_bottom = Constant(10.0)
p2_top = Constant(275.671484)


bcs = []
#=====================================;
# Define normal vector and mesh size  ;
#=====================================;
n = FacetNormal(mesh)
h = CellSize(mesh)
h_avg = (h('+') + h('-'))/2

#penalty parameter
eta_u, eta_p = Constant(100.), Constant(1.)

#eta_Nitsche = Constant(100.0)

#===========================;
#  Define variational form  ;
#===========================;
invalpha1 = 1.0 / alpha1
invalpha2 = 1.0 / alpha2


a = dot(w1, alpha1*v1)*dx + dot(w2, alpha2*v2)*dx \
    - div(w1) * p1 *dx - div(w2) * p2 * dx \
    + q1 * div(v1) * dx + q2 * div(v2) * dx +\
    q1 * fact * (p1 - p2) * dx -\
    q2 * fact * (p1 - p2) * dx -\
    0.5 * dot( alpha1 * w1 - grad(q1), \
               invalpha1 * (alpha1 * v1 + grad(p1)) ) * dx -\
    0.5 * dot( alpha2 * w2 - grad(q2), \
               invalpha2 * (alpha2 * v2 + grad(p2)) ) * dx


L = dot(w1,rhob1)*dx +\
    dot(w2,rhob2)*dx -\
    0.5 * dot( alpha1 * w1 - grad(q1), \
               invalpha1 * rhob1 ) * dx -\
    0.5 * dot( alpha2 * w2 - grad(q2), \
               invalpha2 * rhob2 ) * dx -\
    dot(w1,n) * p1_left * ds(1) -\
    dot(w2,n) * p2_left * ds(1) -\
    dot(w1,n) * p1_right * ds(2) -\
    dot(w2,n) * p2_right * ds(2) -\
    dot(w1,n) * p1_bottom * ds(3) -\
    dot(w2,n) * p2_bottom * ds(3) -\
    dot(w1,n) * p1_top * ds(4) -\
        dot(w2,n) * p2_top * ds(4)
#==================;
#  Solver options  ;
#==================;
if show_monitor:
  monitor_view = True
else:
  monitor_view = False
parameters_standard = {
  "ksp_type": "gmres",
  "pc_type": "bjacobi",
  "ksp_monitor_true_residual": show_monitor,
  "ksp_converged_reason": True,
  "ksp_rtol": rtol
}
parameters_twofields = {
  "ksp_type": "gmres",
  #"mat_type": "aij",
  "pc_type": "fieldsplit",
  # first split contains velocities
  # second split contains pressures
  "pc_fieldsplit_0_fields": "0,2",
  "pc_fieldsplit_1_fields": "1,3",
  # Use schur complement
  "pc_fieldsplit_type": "schur",
  "pc_fieldsplit_schur_fact_type": "full",
  "pc_fieldsplit_schur_precondition": "selfp",
  # Velocities
  "fieldsplit_0_ksp_type": "preonly",
  "fieldsplit_0_pc_type": "bjacobi",
  # Pressures - fieldsplit into scales
  "fieldsplit_1_ksp_type": "preonly",
  "fieldsplit_1_pc_type": "fieldsplit",
  "fieldsplit_1_pc_fieldsplit_type": "additive",
  # Macro scale pressure
  "fieldsplit_1_fieldsplit_0_ksp_type": "preonly",
  "fieldsplit_1_fieldsplit_0_pc_type": amg_type,
  # Only needed if using hypre for 3D problems
  # Comment out for 2D problems
  #"fieldsplit_1_fieldsplit_0_pc_hypre_boomeramg_strong_threshold": 0.75,
  #"fieldsplit_1_fieldsplit_0_pc_hypre_boomeramg_agg_nl": 2,
  # Micro scale pressure
  "fieldsplit_1_fieldsplit_1_ksp_type": "preonly",
  "fieldsplit_1_fieldsplit_1_pc_type": amg_type,
  # Only needed if using hypre for 3D problems
  # Comment out for 2D problems
  #"fieldsplit_1_fieldsplit_1_pc_hypre_boomeramg_strong_threshold": 0.75,
  #"fieldsplit_1_fieldsplit_1_pc_hypre_boomeramg_agg_nl": 2,
  "ksp_monitor_true_residual": monitor_view,
  "ksp_converged_reason": True,
  "ksp_rtol": rtol
}
parameters_twoscales = {
  "ksp_type": "gmres",
  #"mat_type": "aij",
  "pc_type": "fieldsplit",
  # First split contains macro scale
  # Second split contains micro scale
  "pc_fieldsplit_0_fields": "0,1",
  "pc_fieldsplit_1_fields": "2,3",
  "pc_fieldsplit_type": "additive", # Can change this to multiplicative
  # Macro scale Darcy
  "fieldsplit_0_ksp_type": "preonly",
  "fieldsplit_0_pc_type": "fieldsplit",
  "fieldsplit_0_pc_fieldsplit_type": "schur",
  "fieldsplit_0_pc_fieldsplit_schur_fact_type": "full",
  "fieldsplit_0_pc_fieldsplit_schur_precondition": "selfp",
  # Fieldsplitting of macro scale
  "fieldsplit_0_fieldsplit_0_ksp_type": "preonly",
  "fieldsplit_0_fieldsplit_0_pc_type": "bjacobi",
  "fieldsplit_0_fieldsplit_1_ksp_type": "preonly",
  "fieldsplit_0_fieldsplit_1_pc_type": amg_type,
  # Only needed if using hypre for 3D problems
  # Comment out for 2D problems
  #"fieldsplit_0_fieldsplit_1_pc_hypre_boomeramg_strong_threshold": 0.75,
  #"fieldsplit_0_fieldsplit_1_pc_hypre_boomeramg_agg_nl": 2,
  # Micro scale Darcy
  "fieldsplit_1_ksp_type": "preonly",
  "fieldsplit_1_pc_type": "fieldsplit",
  "fieldsplit_1_pc_fieldsplit_type": "schur",
  "fieldsplit_1_pc_fieldsplit_schur_fact_type": "full",
  "fieldsplit_1_pc_fieldsplit_schur_precondition": "selfp",
  # Fieldsplitting of micro scale
  "fieldsplit_1_fieldsplit_0_ksp_type": "preonly",
  "fieldsplit_1_fieldsplit_0_pc_type": "bjacobi",
  "fieldsplit_1_fieldsplit_1_ksp_type": "preonly",
  "fieldsplit_1_fieldsplit_1_pc_type": amg_type,
  # Only needed if using hypre for 3D problems
  # Comment out for 2D problems
  #"fieldsplit_1_fieldsplit_1_pc_hypre_boomeramg_strong_threshold": 0.75,
  #"fieldsplit_1_fieldsplit_1_pc_hypre_boomeramg_agg_nl": 2,
  "ksp_monitor_true_residual": monitor_view,
  "ksp_converged_reason": True,
  "ksp_rtol": rtol
}

#=================;
#  Solve problem  ;
#=================;
solution = Function(wSpace)
with timed_region('solve_firedrake'):
  initialtime = time.time()
  A = assemble(a, bcs=bcs, mat_type='aij')
  b = assemble(L)
  assembletime = time.time()
  # Option 1: default solver
  if solvestrategy == "naive":
    solver = LinearSolver(A,P=None,options_prefix="standard_",solver_parameters=parameters_standard)
  # Option 2: split by fields
  elif solvestrategy == "fields":
    solver = LinearSolver(A,P=None,options_prefix="twofields_",solver_parameters=parameters_twofields)
  # Option 3: split by scales
  elif solvestrategy == "scales":
    solver = LinearSolver(A,P=None,options_prefix="twoscales_",solver_parameters=parameters_twoscales)
  solver.solve(solution,b)
  solvetime = time.time()


#=======================;
#  Performance metrics  ;
#=======================;
with solution.dat.vec as solution_vec:
  DoF = solution_vec.getSize()
totaltime = solvetime - initialtime
dofsec = DoF/totaltime
dofsecLOG = np.log10(dofsec)
totaltimeLOG = np.log10(totaltime)

if rank == 0:
  print("================================")
  print("%d MPI processes" % size)
  print("AMG and solver type: %s %s" % (amg_type,solvestrategy))
  print("Relative residual: %1.1e" % (rtol))
  print("Nx = %d" % nx)
  print("Dof = %d\n" % DoF)
  print("Assembly time = %1.3e seconds" % (assembletime - initialtime))
  print("Solve time = %1.3e seconds" % (solvetime - assembletime))
  print("Totaltime = %1.3e seconds" % totaltime)
  print("DofSec = %1.3e" % dofsec)
  print("DofSecLOG = %1.3e" % dofsecLOG) # used for static scaling
  print("TotaltimeLOG = %1.3e log(seconds)" % totaltimeLOG) # used for static scaling , True static scaling, and DoE plots
  #print("L2 error = %1.3e" % L2error) # True static-scaling
  #print("TrueDof = %1.3e" % truedof) # True static-scaling
  #print("DoE = %1.3e" % DoE) # Digits of efficacy
  print("================================")
#=======================================;
#  Dump solution to file in VTK format  ;
#=======================================;
v1sol,p1sol,v2sol,p2sol = solution.split()

#file = File('Output/v1_2D_CG.pvd')
#file.write(v1sol)

#file = File('Output/v2_2D_CG.pvd')
#file.write(v2sol)

#file = File('Output/p1_2D_CG.pvd')
#file.write(p1sol)

#file = File('Output/p2_2D_CG.pvd')
#file.write(p2sol)

#==========================;
#  Define exact solutions  ;
#==========================;
p1_ex = Function(pSpace)
p2_ex = Function(pSpace)
v1_ex = Function(vSpace)
v2_ex = Function(vSpace)
p1_exact = Expression("(1/pi)*exp(pi*x[0])*sin(pi*x[1]) - (1/(1.0*1.0))*exp(3.316625*x[1])", degree = 5)
p2_exact = Expression("(1/pi)*exp(pi*x[0])*sin(pi*x[1]) + (1/(1.0*0.1))*exp(3.316625*x[1])", degree = 5)	
v1_exact = Expression(("-1*exp(pi*x[0])*sin(pi*x[1])","-1*exp(pi*x[0])*cos(pi*x[1]) + (3.316625)*exp(3.316625*x[1])"), degree = 5)
v2_exact = Expression(("-0.1*exp(pi*x[0])*sin(pi*x[1])","-0.1*exp(pi*x[0])*cos(pi*x[1]) - (3.316625)*exp(3.316625*x[1])"), degree = 5)	
p1_ex = interpolate(p1_exact, pSpace)
v1_ex = interpolate(v1_exact, vSpace)
p2_ex = interpolate(p2_exact, pSpace)
v2_ex = interpolate(v2_exact, vSpace)

#file = File('Output/v1ex_2D_Q4.pvd')
#file.write(v1_ex)

#file = File('Output/v2ex_2D_Q4.pvd')
#file.write(v2_ex)

#file = File('Output/p1ex_2D_Q4.pvd')
#file.write(p1_ex)

#file = File('Output/p2ex_2D_Q4.pvd')
#file.write(p2_ex)

##plot(p2_ex,title = "Exact solution_microP")
L2_p1 = errornorm(p1_ex,p1sol,norm_type='L2',degree_rise= 3)
#H1_p1 = errornorm(p1_ex,p1sol,norm_type='H1',degree_rise=3)
L2_v1 = errornorm(v1_ex,v1sol,norm_type='L2',degree_rise= 3)
#print ("H1 error in p1 = ", H1_p1)
L2_p2 = errornorm(p2_ex,p2sol,norm_type='L2',degree_rise= 3)
#H1_p2 = errornorm(p2_ex,p2sol,norm_type='H1',degree_rise=3)
L2_v2 = errornorm(v2_ex,v2sol,norm_type='L2',degree_rise= 3)
#print ("H1 error in p2 = ", H1_p2)	
#L2error = ...       
#-------------------;              
# CONVERGENCE PLOTS ;
#-------------------;
#DoA = -np.log10(L2error)
DoA_p1 = -np.log10(L2_p1)
DoA_p2 = -np.log10(L2_p2)
DoA_v1 = -np.log10(L2_v1)
DoA_v2 = -np.log10(L2_v2)
DoS = np.log10(DoF)
#--------------;
# TRUE SCALING ;
#--------------;	
#truedof = DoA/DoS*dofsec
truedof_p1 = DoA_p1/DoS*dofsec
truedof_p2 = DoA_p2/DoS*dofsec	
truedof_v1 = DoA_v1/DoS*dofsec
truedof_v2 = DoA_v2/DoS*dofsec
#-----;
# DoE ;
#-----;
#DoE = -np.log10(L2error*totaltime)
DoE_p1 = -np.log10(L2_p1*totaltime)
DoE_p2 = -np.log10(L2_p2*totaltime)
DoE_v1 = -np.log10(L2_v1*totaltime)
DoE_v2 = -np.log10(L2_v2*totaltime)

if rank == 0:
	print ("L2 error in p1 =", L2_p1)
	print ("L2 error in v1 =", L2_v1)
	print ("L2 error in p2 =", L2_p2)
	print ("L2 error in v2 =", L2_v2)
	print ("DoA_p1 =" , DoA_p1)
	print ("DoA_p2 =" , DoA_p2)
	print ("DoA_v1 =" , DoA_v1)
	print ("DoA_v2 =" , DoA_v2)
	print ("DoS =" , DoS)
	print ("truedof_p1 =" , truedof_p1)
	print ("truedof_p2 =" , truedof_p2)
	print ("truedof_v1 =" , truedof_v1)
	print ("truedof_v2 =" , truedof_v2)
	print ("DoE_p1 =" , DoE_p1)
	print ("DoE_p2 =" , DoE_p2)
	print ("DoE_v1 =" , DoE_v1)
	print ("DoE_v2 =" , DoE_v2)
	
