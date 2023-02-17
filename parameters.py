import numpy as np
NJ = 7
NM = 7
# Define scalars

Cp_steel       = 720.0  # J/kg K # not in use ?
Cp_slag        = 500.0  # J/kg K
Cp_dolomite    = 1.0*1500.0 # 2.0*1100.0
Cp_alumina     = 718.0 # J/kg K - https://www.engineeringtoolbox.com/specific-heat-capacity-d_391.html
Cp_mgoc        = 1500.0 #1.0*1000.0 # 2.0*1000.0 # 1000.0 # https://www.ispatguru.com/magnesia/
# Cp_mgoc        = 1300.0 # Santos et.al 2018
Cp_insulation  = 900.0
rho_steel      = 7100.0 # not in use ?
# rho_steel      = 8100.0
rho_slag       = 3500.0  # 8100.0
rho_dolomite   = 2900.0
rho_alumina    = 2500.0
rho_mgoc       = 3540.0 #  https://www.ispatguru.com/magnesia/ 
rho_insulation = 300.0
rho_carbon     = 2250.0
rho_mgo        = rho_mgoc*1.1  # FIXME: must find correct value, but must be larger than mgoc
lam_steel      = 15.0
lam_steelshell = 12.0
lam_slag       = 10.0
lam_dolomite   = 2.7
# lam_alumina    = 1.2*0.7
lam_alumina    = 2.0
# lam_mgoc       = 6.0
lam_mgoc       = 30.0  # material-properties.org says 1.31, cp  800, den = 1700
lam_insulation = 0.1 # 0.02
# lam_insulation = 0.02

Cp_wall = np.zeros(NJ)
Cp_wall[6] = Cp_steel
Cp_wall[5] = Cp_insulation
Cp_wall[4] = Cp_alumina
Cp_wall[3] = Cp_dolomite
Cp_wall[2] = Cp_mgoc
Cp_wall[1] = Cp_mgoc
Cp_wall[0] = Cp_mgoc
Cp_wall = Cp_wall/1.0
rho_wall = np.zeros(NJ)
rho_wall[6] = rho_steel
rho_wall[5] = rho_insulation
rho_wall[4] = rho_alumina
rho_wall[3] = rho_dolomite
rho_wall[2] = rho_mgoc
rho_wall[1] = rho_mgoc
rho_wall[0] = rho_mgoc
lam_wall = np.zeros(NJ)
lam_wall[6] = lam_steelshell
lam_wall[5] = lam_insulation
lam_wall[4] = lam_alumina
lam_wall[3] = lam_dolomite
lam_wall[2] = lam_mgoc
lam_wall[1] = lam_mgoc
lam_wall[0] = lam_mgoc
lam_wall = lam_wall*1.0
Cp_bottom = np.zeros(NM)
Cp_bottom[0] = Cp_steel
Cp_bottom[1] = Cp_insulation
Cp_bottom[2] = Cp_alumina
Cp_bottom[3] = Cp_dolomite
Cp_bottom[4] = Cp_mgoc
Cp_bottom[5] = Cp_mgoc
Cp_bottom[6] = Cp_mgoc
Cp_bottom = Cp_bottom/1.0
rho_bottom = np.zeros(NJ)
rho_bottom[0] = rho_steel
rho_bottom[1] = rho_insulation
rho_bottom[2] = rho_alumina
rho_bottom[3] = rho_dolomite
rho_bottom[4] = rho_mgoc
rho_bottom[5] = rho_mgoc
rho_bottom[6] = rho_mgoc
lam_bottom = np.zeros(NJ)
lam_bottom[0] = lam_steelshell
lam_bottom[1] = lam_insulation
lam_bottom[2] = lam_alumina
lam_bottom[3] = lam_dolomite
lam_bottom[4] = lam_mgoc
lam_bottom[5] = lam_mgoc
lam_bottom[6] = lam_mgoc
lamb = np.zeros((NM+1))
# lam_bottom = lam_bottom*1e-10
# lam_wall = lam_wall*1e-10

alpha_carbon = 0.05  # FIXME: Need proper value
alpha_mgo    = 1.0 - alpha_carbon