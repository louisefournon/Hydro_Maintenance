
# Bundle Method first implementation

### I. Preliminary work

import numpy as np
import random
from docplex.util.status import JobSolveStatus
import docplex.mp.model as cpx
import sys 

##############################
#      Define the oracles    #
##############################

# First toy problem oracles

def oracle1(lamb):
  x_1_star = 0
  x_2_star = 0
  eval1 = -1 + 2*lamb[0] + lamb[1] + 2*lamb[2]
  if( eval1 <= 0 ):
    x_1_star = 5
  eval2 = -1 + lamb[0] - 2*lamb[1] + 3*lamb[2]
  if( eval2 <= 0 ):
    x_2_star = 5
  f_val = -eval1*x_1_star -eval2*x_2_star + 4*lamb[0] + lamb[1] + 6*lamb[2]
  sg = [4 - 2*x_1_star - x_2_star, 1 - x_1_star + 2*x_2_star, 6 - 2*x_1_star - 3*x_2_star]
  return [f_val, sg]
  
  
# Second toy problem oracles

# Here the objective function is of the form: f(x) = max_{1 <= i <= kPieces} x^T Q_(i) x - b_i^T x

# nbVars = 10
# timeHorizon = 10
# kPieces = 5
# rnd = 10000
# 
# Q = []
# b = []
# 
# for i in range(kPieces):
#   Q.append(np.zeros((timeHorizon, timeHorizon)))
#   b.append(np.zeros(timeHorizon))
#   
#   for j in range(timeHorizon):
#     for k in range(j+1, timeHorizon):
#       Q[i][j][k] = np.exp((j+1)/(k+1))*np.cos((j+1)*(k+1))*np.sin(i+1)
#       Q[i][k][j] = Q[i][j][k]
#       
#     Q[i][j][j] = np.abs(np.sin(j+1))*((i+1)/timeHorizon) + sum( np.abs(Q[i][j][l]) for l in range(timeHorizon) if l != j)
#     b[i][j] = np.exp((j+1)/(i+1))*np.sin((i+1)*(j+1))
# 
# def oracle2(x):
#   j0 = -1
#   f_val = np.NINF
#   for j in range(kPieces):
#     pp = np.dot(x, np.dot(Q[j], x)) - np.dot(b[j], x)
#     if(f_val < pp):
#       j0 = j
#       f_val = pp
#   sg = 2*np.dot(Q[j0],x) - b[j0]
#   return [f_val, sg]

# Third problem: thermal unit

demand = np.array([2288.4, 2590.1, 3002.8, 3165.6,
 3313.9, 3454.0, 3314.7, 3073.9,
 2518.0, 2034.4, 1610.9, 1010.0,
 1613.6, 2034.7, 2515.4, 3066.5,
 3300.5, 3451.8, 3300.0, 3066.2,
 2521.4, 2047.9, 1627.2, 1029.1])

nbThermal = 9
T = 24
dt = 2

initP = np.zeros(nbThermal)
initP[0] = 700
initP[1] = 700 
initP[2] = 700
initP[3] = 150
initP[4] = 150

therm_cost = np.zeros(nbThermal)
therm_cost[0] = 30
therm_cost[1] = 35
therm_cost[2] = 37
therm_cost[3] = 45
therm_cost[4] = 47
therm_cost[5] = 60
therm_cost[6] = 100
therm_cost[7] = 110
therm_cost[8] = 150

pow_max = np.zeros(nbThermal)
pow_max[0] = 900
pow_max[1] = 900
pow_max[2] = 900
pow_max[3] = 300
pow_max[4] = 300
pow_max[5] = 200
pow_max[6] = 200
pow_max[7] = 200
pow_max[8] = 100


therm_grad = np.zeros(nbThermal);
therm_grad[0] = 100
therm_grad[1] = 100
therm_grad[2] = 100
therm_grad[3] = 30
therm_grad[4] = 30
therm_grad[5] = 20
therm_grad[6] = 20
therm_grad[7] = 20
therm_grad[8] = 10

b_inf = []
b_sup = []
for i in range(nbThermal):
  b_inf_i = -therm_grad[i]*dt*np.ones(T)
  b_sup_i = therm_grad[i]*dt*np.ones(T)
  b_inf_i[0] += initP[i]
  b_sup_i[0] += initP[i]
  b_inf.append(b_inf_i)
  b_sup.append(b_sup_i)
  

# x is the production level, lamb the dual variable and i the index of the unit

def oracleTherm(lamb, iCentrale, therm_grad, therm_cost, pow_max, initP, T, dt):
  
  # Define parameters
  A_therm = np.eye(T)
  for j in range(T - 1):
    A_therm[j+1][j] = -1
  
  b_inf = -therm_grad[iCentrale]*dt*np.ones(T)
  b_sup = therm_grad[iCentrale]*dt*np.ones(T)
  b_inf[0] += initP[iCentrale]
  b_sup[0] += initP[iCentrale]
  
  # Create model
  opt_model = cpx.Model(name="Sub Problem Thermal")
  opt_model.parameters.qpmethod = 2
  
  # Create constraints parameters
  lin_ct1 = { (i,j): A_therm[i][j] for i in range(T) for j in range(T) }
  const_ct2 = { j: b_inf[j] for j in range(T) }
  const_ct3 = { j: b_sup[j] for j in range(T) }    
    
  # Add decision variables
  
  # There is only one asset for thermal units and x_vars corresponds to the energy production  
  x_vars  = {t: opt_model.continuous_var(ub = pow_max[iCentrale]) for t in range(T)}

  # Set constraints
  sup_constraints = {j : opt_model.add_constraint(
    ct=opt_model.sum(A_therm[i][j] * x_vars[j] for j in range(T)) <= b_sup[i]) for i in range(T)}
  inf_constraints = { j : opt_model.add_constraint(
    ct=opt_model.sum(A_therm[i][j] * x_vars[j] for j in range(T)) >= b_inf[i]) for i in range(T)}
  
  # Set objective
  cost = therm_cost[iCentrale]*np.ones(T) - lamb[0]
  objective = opt_model.sum(cost[i]*x_vars[i]*dt for i in range(T))
    
  opt_model.minimize(objective)
  
  # Solve
  opt_model.solve()
  
  # We check wether or not the SP could be solved
  # opt_model.export_as_lp(path="C:/Users/Louise Fournon/Documents/ENPC/MPRO/Stages/StageEdimbourg", basename=None, hide_user_names=False)
  
  if(opt_model.get_solve_status() == JobSolveStatus.OPTIMAL_SOLUTION):
    aPower = np.array([opt_model.solution[x_vars[i]] for i in range(T)])
    obj = opt_model.objective_value
    return [obj, aPower]


#  Fourth problem : hydro units  

nMaintenance = 2

# Valley de l'Ain

hm3tom3 = 1e6
nRes_1 = 6

A_1 = np.zeros((nRes_1,nRes_1))
A_1[0][1] = 1
A_1[1][2] = 1 
A_1[2][4] = 1 
A_1[3][4] = 1
A_1[4][5] = 1

Vmin_1 = hm3tom3 * np.array([150.0, 0.96, 33.27, 0.88, 10.81, 17.6])
V0_1 = hm3tom3 * np.array([277.5, 1.011, 34.38, 1.84, 12.90, 18.1])
Vmax_1 = hm3tom3 * np.array([ 417.8, 1.5, 35.12, 4.38, 13.59, 18.2])

mxFlow_1 = np.array([72.5, 100.0, 120.0, 15.0, 60.0, 90.0])
mxPow_1 = np.array([65.0, 22.0, 20.0, 12.0, 7.0, 14.0])

wvals_1 = np.array([40.0, 35.0, 30.0, 33.0, 28.0, 25.0])
nominf_1 = np.array([21.2, 0.0, 5.0, 3.0, 2.0, 0.0])

nbTurbine_1 = 16
sigT_1 = np.array([0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5])

nbHyd = 2

# Plug in the Isere Valley

nRes_2 = 7
A_2 = np.zeros((nRes_2, nRes_2))
A_2[0][1] = 1
A_2[1][3] = 1
A_2[3][5] = 1
A_2[4][5] = 1
A_2[2][6] = 1
A_2[5][6] = 1

Vmin_2 = hm3tom3 * np.array([4.0, 3, 0.1, 0, 0, 0.1, 0.09])
V0_2 = hm3tom3 * np.array([ 5.0, 133.43, 0.5,0.05,0.05,0.3,0.1])
Vmax_2 = hm3tom3 * np.array([9.78, 223.83,1.32,0.1,0.1,0.6,0.3])

mxFlow_2 = np.array([11.0, 16.0, 10.0, 12.0, 3.3, 47.7, 25.0])
mxPow_2 = np.array([24.5, 32.0, 75.0, 75.0, 9.6, 17.7, 30.0])

wvals_2 = np.array([35, 40, 35, 30, 33, 28, 25])
nominf_2 = np.array([0, 3.9, 0.8, 0.2, 0.2, 1.5, 7])

nbTurbine_2 = 23
sigT_2 = np.array([0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 5, 5, 5, 6, 6, 6, 6])


def rho(x, iR, mxPow, mxFlow):
  return mxPow[iR]/np.log(1.0 + mxFlow[iR])*np.log(1.0 + x)

def efficiency(nRes, sigma_T, mxFlow, mxPow, wvals, nbTurbine):
  rhoEff = np.zeros( nbTurbine )
  turBnd = np.zeros( nbTurbine )
  iOff = 0
  for i in range(nRes):
    iRes = sigma_T[iOff]
    nbTur = 0
    for j in range(iOff, nbTurbine):
      if (iRes != sigma_T[j]) or (j==(nbTurbine-1)):
        if ( j == (nbTurbine-1) ):
          nbTur = j - iOff + 1
        else:
          nbTur = (j-1) - iOff + 1
        break
    deltaP = mxFlow[iRes]/nbTur;
    for j in range(iOff, iOff + nbTur):
      rhoEff[j] = (1.0/3600.0)*(rho((j-iOff +1)*deltaP, iRes, mxPow, mxFlow) - rho((j-iOff)*deltaP, iRes, mxPow, mxFlow)) / deltaP
      turBnd[j] = deltaP
    iOff += nbTur
  
    # We have to scale the watervalues to be consistent (they are in euro/MW) and not in euro/m3
    wvals[i] *= (mxPow[i]/mxFlow[i])*(1.0/3600.0)
  return [rhoEff, turBnd]
  
def upper_bound(T, nbTurbine, nRes, sigma_T, mxFlow, mxPow, wvals, Vmax, Vmin):
  
  nbVars = T*nbTurbine + 2*nRes
  x_u = np.zeros(nbVars)
  turBnd = efficiency(nRes, sigma_T, mxFlow, mxPow, wvals, nbTurbine)[1]
  for i in range(T):
      for j in range(nbTurbine):
        x_u[i*nbTurbine + j] = turBnd[j]*3600.0 # m3/h

  for i in range(nRes):
    x_u[ T*nbTurbine + i ] = Vmax[i] - Vmin[i]
    x_u[T*nbTurbine + nRes + i] = Vmax[i] - Vmin[i]
  return x_u
  
def turbine_belong(iR, sigma_T, nbTurbine):
# Figure out the turbines belonging to this reservoir (reservoir iR)
  kOff = 0
  nbTur = 0
  for k in range(nbTurbine):
    if ( sigma_T[k] == iR ): 
      # kOff is the number of the first turbine belonging to reservoir iR
      kOff = k
      break 
  nbTur = 0
  for j in range(kOff, nbTurbine):
    if ((iR != sigma_T[j]) or (j==(nbTurbine-1))):
      if ( j == ( nbTurbine - 1 ) ):
        nbTur = j - kOff + 1
      else:
        nbTur = (j-1) - kOff + 1
      # nbTur is the number of turbines belonging to reservoir iR
      break
  return [int(kOff), int(nbTur)]
  
def set_hydro_constraints(opt_model, x_vars, xi_vars, x_u, nRes, nbVars, nbTurbine, T, inflow_nom, V0, Vmin, Vmax, dt, sigT, A_connect, z):
    
  AnRes = np.zeros(nRes)
  AnOff = np.zeros(nRes)
  AnNum = np.zeros(nRes)
  a_ = np.zeros(nbVars)
  for jRes in range(nRes):
    a_ = np.zeros(nbVars)
    pair = turbine_belong(jRes, sigT, nbTurbine)
    iOff = pair[0]
    nbTur = pair[1]
    # Which are the ancestor reservoirs		
    AnRes = np.zeros( nRes )
    nbAn = 0
    for kRes in range(nRes):
      if ( A_connect[kRes][jRes] == 1 ):
        AnRes[nbAn] = kRes
        pair = turbine_belong(kRes, sigT, nbTurbine)
        AnOff[nbAn] = pair[0]
        AnNum[nbAn] = pair[1]
        nbAn += 1
    for i in range(T):
      # Remove turbined quantities at this time step			
      for j in range(iOff, iOff + nbTur):
          a_[ i*nbTurbine + j ] = -1.0*dt
      #  dd Amont Turbined stuff
      # Flow delay should be added here in a more sophisticated model
      for kAn in range(nbAn):
        for k in range(int(AnOff[kAn]), int(AnOff[kAn] + AnNum[kAn])):
          a_[i*nbTurbine + k ] = 1.0*dt
          
      # Set the maintenance constraints
      nbMainTurb = int(len(z)/nMaintenance)
      for i in range(nbMainTurb):
        for k in range(nMaintenance):
          maint_ct = { t : opt_model.add_constraint(x_vars[i*T + t] <= x_u[T*i + t]*xi_vars[nMaintenance*i + k]) for t in range(int(k*T/nMaintenance), int((k+1)*T/nMaintenance)) }
      
    
      # Add the Constraint
      # flowMatrix * x <= VmaxB
      # VmaxB = Vmax - V0 - cumulated inflows
      VmaxB = Vmax[jRes]- V0[jRes] - inflow_nom[jRes]*3600.0*dt*(i+1)		
      mxFlow_ct = opt_model.add_constraint( 
        ct=opt_model.sum(a_[i] * x_vars[i] for i in range(nbVars)) <= VmaxB) 
                  
      # flowMatrix * x >= VminB
      # VminB = Vmin - V0 - cumulated inflows
      VminB = Vmin[jRes] - V0[jRes] - inflow_nom[jRes]*3600.0*dt*(i+1)
      minFlow_ct = opt_model.add_constraint( 
        ct=opt_model.sum(a_[i] * x_vars[i] for i in range(nbVars)) >= VminB) 

    # Before moving to the next reservoir we add the water value constraints
    # Add constraints related to zF
    # At this stage a_ contains the last line of the flow equations
    a_[T*nbTurbine + nRes + jRes] = -1.0
    bz = Vmin[jRes] - V0[jRes] - inflow_nom[jRes]*3600.0*dt*T

    flow_ct1 = opt_model.add_constraint(
      ct=opt_model.sum(a_[i] * x_vars[i] for i in range(nbVars)) == bz) 
    
    # Add z0 related constraint
    bz = V0[jRes] - Vmin[jRes]
    a_ = np.zeros(nbVars)
    a_[ T*nbTurbine + jRes ] = 1.0
    
    flow_ct2 = opt_model.add_constraint( 
      ct=opt_model.sum(a_[i] * x_vars[i] for i in range(nbVars)) == bz)

# Here we consider a whole hydro valley

# The "current_point" is of the form [ x[i,t],..., V_i[0],...,V_i[T] ]


def oracleHydro(A_connect, lamb, T, dt, V0, Vmin, Vmax, nRes, nbTurbine, mxFlow, mxPow, sigT, wvals, nominf, z = None):

  nbVars1 = T*nbTurbine + 2*nRes
  nbVars2 = len(z)
  # Create model
  opt_model = cpx.Model(name="Sub Problem Hydro")
  opt_model.parameters.qpmethod = 2
  
  # Define parameters
  rhoEff = efficiency(nRes, sigT, mxFlow, mxPow, wvals, nbTurbine)[0]
  c = np.zeros(nbVars1)
  
  for i in range(T):
    for j in range(nbTurbine): 
      c[i*nbTurbine + j] = -1.0*rhoEff[j]*lamb[0][i]*dt
          
  for i in range(nRes):
    c[T*nbTurbine + i] = wvals[i]
    c[T*nbTurbine + nRes + i] = -1.0*wvals[i]	
    
  # Add decision variables 
  
  # For hydro units, the decision variables are the volume in each reservoir and the flow for each turbine
  x_u = upper_bound(T, nbTurbine, nRes, sigT, mxFlow, mxPow, wvals, Vmax, Vmin)
  if not z is None: # z should have nTurbine rows and 2 columns for each period
    xi_vars = { i : opt_model.binary_var() for i in range(nbVars2) }
    
  # We represent the variable in the following form : [x[i,t],..., z_i[0],...,z_i[T]]
  x_vars = { i : opt_model.continuous_var(ub = x_u[i]) for i in range(nbVars1) }
    
  # Set constraints
  set_hydro_constraints(opt_model, x_vars, xi_vars,  x_u, nRes, nbVars1, nbTurbine, T, nominf, V0, Vmin, Vmax, dt, sigT, A_connect, z)
  
  # Set objective
  objective = opt_model.sum(x_vars[i] * c[i] for i in range(nbVars1)) - opt_model.sum(lamb[1][i]* xi_vars[i] for i in range(len(z)))
  opt_model.minimize(objective)

  # Solve
  opt_model.solve()
  
  # We check wether or not the SP could be solved
  opt_model.export_as_lp(path="C:/Users/Louise Fournon/Documents/ENPC/MPRO/Stages/StageEdimbourg", basename=None, hide_user_names=False)
  
  if(opt_model.get_solve_status() == JobSolveStatus.OPTIMAL_SOLUTION):
    val = np.array([opt_model.solution[x_vars[i]] for i in range(nbVars1)])
    xi = np.array([opt_model.solution[xi_vars[i]] for i in range(nbVars2)])
    aPower = np.zeros(T)
    for i in range(T):
        for j in range(nbTurbine):
          aPower[i] += rhoEff[j]*val[i*nbTurbine + j]
    obj = opt_model.objective_value
    return [obj, aPower, xi]


# Lambda in argument has dimension T*(nbPbThermal + nbPbHydro)

def lagrangian(nbPbTherm, nbPbHydro, T, dt, lamb, A_connect = None, V0 = None, Vmin = None, Vmax = None, nRes = None, nbTurbine = None, mxFlow = None, mxPow = None, sigT = None, wvals = None, nominf = None, therm_grad = None, therm_cost = None, pow_max = None, initP = None, z = None):
  theta = 0
  sg1 = np.zeros(T)
  sg2 = np.zeros(len(z))
  xi = np.zeros(len(z))
  if(nbPbTherm > 0):
    for i in range(nbPbTherm):
      oracle = oracleTherm(lamb, i, therm_grad, therm_cost, pow_max, initP, T, dt)
      theta += oracle[0] # Objective
      sg1 -= oracle[1]*dt # Active power (size T)
      print("oracle[0] = ", oracle[0])
      print("Power = ", oracle[1])
  if(nbPbHydro > 0):
    for i in range(nbPbHydro):
      oracle = oracleHydro(A_connect[i], lamb, T, dt, V0[i], Vmin[i], Vmax[i], nRes[i], nbTurbine[i], mxFlow[i], mxPow[i], sigT[i], wvals[i], nominf[i], z)
      theta += oracle[0] # Objective
      sg1 -= oracle[1]*dt # Active power (size T)
      if(len(oracle) > 2):
        theta += np.dot(lamb[1], z - oracle[2])
        sg2 += z - oracle[2]
      print("oracle[0] = ", oracle[0])
      print("Power = ", oracle[1])
  theta += np.dot(demand, lamb[0])
  sg1 += demand*dt
  return [theta, sg1, sg2]



#######################################
#    Define intermediate quantities   #
#######################################   

# Returns delta_k and the j index yielding the minimum in delta

def compute_delta(function_bundle, f_low): 
  mini = min(function_bundle)
  index = np.argmin(function_bundle)
  return [index, mini - f_low]


##################################
#       Choose next iterate      #     
##################################

# Let nbRows be the number of rows of the \tilde A matrix (i.e. also the number of components of \tilde b)


def find_next_iterate(stab_center, f_lev, function_bundle, subgradient_bundle, iterates, T, ub = None, lb = None, A_tilde = None, b_1 = None, b_2 = None):
  
  nbVars1 = len(stab_center[0])
  nbVars2 = len(stab_center[1])
  # Create model
  opt_model = cpx.Model(name="Sub Problem")
  opt_model.parameters.qpmethod = 2
  
  # Ct1 = cutting planes model level constraint
  # lin_ct1 = { (j,t) : subgradient_bundle[j][t] for t in range(T) for j in range(len(iterates)) }
  # 
  # const_ct1 = { j : f_lev - function_bundle[j] + np.dot(subgradient_bundle[j][0],iterates[j][0]) + np.dot(subgradient_bundle[j][1],iterates[j][1]) for j in range(len(iterates)) }
  

  # Add decision variables
  # Actually these are the lambda variables
  # x_vars = { t : opt_model.continuous_var() for t in range(T) }
  
  # Add decision variables
  lamb1_vars = { t : opt_model.continuous_var() for t in range(nbVars1) }
  lamb2_vars = { j : opt_model.continuous_var() for j in range(nbVars2) }
  
  # Set constraints
  cp_constraints = { j : opt_model.add_constraint( 
    ct = function_bundle[j] + np.dot(lamb1_vars - iterates[j][0],subgradient_bundle[j][0]) + np.dot(lamb2_vars - iterates[j][1],subgradient_bundle[j][1])) for j in range(len(iterates))}
  
  # cp_constraints = {j : opt_model.add_constraint(
  #     ct=opt_model.sum(lin_ct1[j,t] * x_vars[t] for t in range(T)) <= const_ct1[j]) for j in range(len(iterates))}

  
  # Set objective
  # objective = opt_model.sum(1/2*(x_vars[i] - stab_center[i])**2 for i in range(T))
   
  objective = opt_model.sum((lamb1_vars[i] - stab_center[0][i])**2 for i in range(nbVars1)) + 1/2*opt_model.sum((lamb2_vars[j] - stab_center[0][j])**2  for j in range(nbVars2))
    
  opt_model.minimize(objective)
  
  #Solve
  opt_model.solve()
  print("opt_model.get_solve_status() = ", opt_model.get_solve_status())
  #We check wether or not the SP could be solved
  if(opt_model.get_solve_status() != JobSolveStatus.OPTIMAL_SOLUTION):
    isEmptyL = True
    return [isEmptyL, stab_center]
  else:
    isEmptyL = False
    lamb1 = np.array([opt_model.solution[lamb1_vars[i]] for i in range(nbVars1)])
    lamb2 = np.array([opt_model.solution[lamb2_vars[i]] for i in range(nbVars2)])
    lamb = [lamb1, lamb2]
    obj = opt_model.objective_value
    return [isEmptyL, lamb]

def cutting_planes_model(function_bundle, subgradient_bundle, iterates, z): 
  hat_W = map(lambda w, xi, zk: w + np.dot(xi,(z-zk)), function_bundle, subgradient_bundle, iterates)
  return max(hat_W)


### II. Algorithm"""

######################################
#      Format of the data            #
######################################

# We store the past evaluations of \bar W (resp. \partial \bar W) in the list
# function_bundle (resp. subgradient_bundle) 
# This way, function_bundle[i] = \bar W(z_i) and subgradient_bundle[i] = \bar \partial W(z_i)

def bundle_method(dt, T, lamb_0, nbPbTherm, nbPbHydro, A_connect, V0, Vmin, Vmax, nRes, nbTurbine, mxFlow, mxPow, sigT, wvals, nominf, therm_grad, therm_cost, pow_max, initP, z = None, ub = None, lb = None, A_tilde = None, b_1 = None, b_2 = None):

  # Data
  gamma = 0.2
  tol = 500
  f_low = -10000000
  f_lev = 1000000

  # initialization 
  k = 0
  nbVars1 = len(lamb_0[0])
  nbVars2 = len(lamb_0[1])

  iterates = [[lamb_0[0], lamb_0[1]]]

  oracle = lagrangian(nbPbTherm, nbPbHydro, T, dt, lamb_0, A_connect, V0, Vmin, Vmax, nRes, nbTurbine, mxFlow, mxPow, sigT, wvals, nominf, therm_grad, therm_cost, pow_max, initP, z)
  
  function_bundle = [-oracle[0]]
  subgradient_bundle = [-oracle[1]]
  delta = tol + 1
  best_index = 0
  sys.exit()
  while (delta > tol and k < 250):

    print("k = ", k)
    print("f_low = ", f_low)
    print("Function evaluation = ", oracle[0]) 
    print("Number of iterates = ", len(iterates))
    delta_pair = compute_delta(function_bundle, f_low)
    delta = delta_pair[1]
    
    print("delta = ", delta)
    
    best_index = delta_pair[0]
    f_lev = f_low + gamma*delta
    next_it_pair = find_next_iterate(iterates[best_index], f_lev, function_bundle, subgradient_bundle, iterates, T, z, ub = None, lb = None, A_tilde = None, b_1 = None, b_2 = None)
    
    next_it = next_it_pair[1]
    isEmptyL = next_it_pair[0]

    if (isEmptyL):
      f_low = f_lev
    else :
      
      iterates.append(next_it)
      oracle = lagrangian(nbPbTherm, nbPbHydro, T, dt, iterates[-1], A_connect, V0, Vmin, Vmax, nRes, nbTurbine, mxFlow, mxPow, sigT, wvals, nominf, therm_grad, therm_cost, pow_max, initP)

      function_bundle.append(-oracle[0])
      print("Function evaluation = ", oracle[0]) 
      subgradient_bundle.append(-oracle[1])
      sys.exit()
    k = k+1
    obj = delta + f_low
    print("f_lev = ", f_lev)

  return [iterates[best_index], obj]


nbPbTherm = nbThermal
nbPbHydro = 2

A_connect = [A_1, A_2]
V0 = [V0_1, V0_2]
Vmin = [Vmin_1, Vmin_2]
Vmax = [Vmax_1, Vmax_2]
nRes = [nRes_1, nRes_2]
nbTurbine = [nbTurbine_1, nbTurbine_2]
mxFlow = [mxFlow_1, mxFlow_2]
mxPow = [mxPow_1, mxPow_2]
sigT = [sigT_1, sigT_2]
wvals = [wvals_1, wvals_2]
nominf = [nominf_1, nominf_2]


# lamb_0 = 40*np.ones(T)
z = np.array([1,1])

lamb1_0 = np.array([29.99175809, 30.02809364, 43.97303945, 36.99451812, 45.99256514, 58.9687617 , 46.03385174, 47.06288056, 35.9554747 , 30.01260556, 34.94836093, 21.56113425, 33.4998073 , 29.99090791, 37.33781095, 47.17917975, 45.50809716, 54.39832005, 45.65966851, 46.95971906, 38.06938898, 30.00898517, 32.97738482, 27.02073265])

lamb2_0 = np.zeros(len(z))

lamb_0 = [lamb1_0, lamb2_0]
print(sys.version_info)
print(bundle_method(dt, T, lamb_0, nbPbTherm, nbPbHydro, A_connect, V0, Vmin, Vmax, nRes, nbTurbine, mxFlow, mxPow, sigT, wvals, nominf, therm_grad, therm_cost, pow_max, initP, z))