# -*- coding: utf-8 -*-
"""
In this script is implemented the model DEoS, a predictive model of explosive
detonation velocity and pressure based on first order
approximation of the detonation velocity equation. Please see:
    
Predictive model of explosive detonation parameters from an equation of state
based on detonation velocity. Fernando G. Bastante*, María Araújo* and Eduardo
Giráldez*.DOI: 10.1039/D2CP00085G Phys. Chem. Chem. Phys., 2022, 24, 8189-8195

*CINTECX, GESSMin Group, Department of Natural Resources
and Environmental Engineering (EEME), University of Vigo, Pontevedra, Spain

Data Inputs (example):
    name = "PBX9502" # the name of explosive without blanks
    rho = 1.90  # density in g/cc (must be >1g//cc)
    C = 2.30    # number of moles
    H = 2.23    # number of moles
    N = 2.21    # number of moles
    O = 2.21    # number of moles 
    W = 3.81    # weight percent of other elements (Cl/F/P/Si % must be <10)
    HR = -205.5 # enthalpy standard of formation in cal/g (1 atm, 298 K)
    filename = "PBX9502_1.90" # filename for saving results

The script can read datainputs from cli...:
    python deos.py name rho C H N O W HR outputfilename 
example:
    python deos.py PBX9502 1.90 2.30 2.23 2.21 2.21 3.81 -205.5 out_pbx9052 
the units are:
        ........  string, g/cc, mol, mol, mol, mol, %weight, cal/g, string
        
...or you can introduce them in this script (see below) 

@author: Fernando García Bastante
"""
import numpy as np
import sys

# if not datainput in cli then introduce (change) then below... 

name = "PBX9502" 
rho = 1.90  # density in g/cc (must be >1g//c)
C = 2.30   # number of moles 
H = 2.23    # number of moles 
N = 2.21    # number of moles 
O = 2.21    # number of moles 
W = 3.81    # weight percent of other elements (must be <10)
HR = -205.5 # enthalpy standard of formation in cal/g (1 atm, 298 K)
filename = name + "__" + str(rho) # filename for saving results

# reading arguments in cli
if len(sys.argv) == 10:
    name = sys.argv[1]
    rho = float(sys.argv[2])
    C = float(sys.argv[3])
    H = float(sys.argv[4])
    N = float(sys.argv[5])
    O = float(sys.argv[6])
    W = float(sys.argv[7])
    HR = float(sys.argv[8])
    filename = sys.argv[9]
    print("Using cli arguments")
    NAMES = ['script=', 'name =', 'rho (g/cc) =', 'C =', 'H =', 'N =', 'O =',\
             'weight percent other elements =', \
             'enthalpy standard (cal/g) =', 'output filename =']  
    for i, j in zip(NAMES, sys.argv):
        print(i,j)
else:
    print('no (valid) arguments in cli, using data in script for: ' + name)

#.............................................................
# constants
# C, H, N, O)
MWI = np.array([12.0111, 1.0079, 14.0067, 15.9994])
# H2O, CO2, C, O2, H2,	 N2
MWP = np.array([2* MWI[1] + MWI[3], MWI[0] + 2*MWI[3], \
                MWI[0], 2 * MWI[3], 2* MWI[1], 2* MWI[2] ])
HP = np.array([-57.8, -94.05, 10, 0, 0, 0])
COV = np.array([75, 123,	 18, 72, 43, 112])
GAMo = 1.23
F = 4184
#.............................................................
# functions
def typex(C, H, O, W):
    
    if W!= 0:
        typ = 0
        return typ

    if O - H / 2 - 2 * C >= 0:
            typ = 1
    else:
        if O - H / 2 - C > 0:
                typ = 2
        else:
            if O - H / 2 > 0:
                    typ = 3
            else:
                if C >= O:
                        typ = 4
                else:
                        typ = 5
    return typ

def jerar(C, H, N, O):

    product = np.empty([6], dtype='float')            
    H2O = np.min([O, H/2])
    CO2 = np.min([0.5*O - 0.5*H2O, C])
    H2 = 0.5 * (H-2*H2O)
    Cs = C - CO2
    O2 = (O-2*CO2-H2O)*0.5
    N2 = N/2
    product = np.array([H2O, CO2, Cs, O2, H2, N2])
    fur = np.where(((H == 0) & (C != 0) & (N != 0) & (O != 0)), 1, 0)
    return product, fur

def mw_hp(product, fur, W, HR):
    
    mwp = np.dot(product, MWP)
    mw = mwp/(1-W/100)
    q = -((np.dot(product, HP) *1000)/mw - HR)
    if q<0:
        raise ValueError ("negative detonation heat, check the data input")
    w = np.dot(product, COV)/mwp
    if fur:
        w = w * 1.07
    Do = (2*q*F*(GAMo**2-1))**0.5  
    return np.array([mw, q, w, Do])

def calc_D_P(Do, q, w, rho):
    
    Q = F * q 
    D = Do + w*rho*1000
    A = w*rho*1000 / D
    B = 0.5 * D * Do / Q               
    gcj = A + np.sqrt(np.square(1 + A) + B)
    acj = np.sqrt(1 + B / np.square(1+A))-1
    bcj = (1 + acj)/(acj*gcj)
    P = rho * np.power(D, 2)/((1 + gcj) * 1e6)
    rocj = 1000 * rho * (gcj+1)/gcj            
    return D, P, Q, gcj, acj, bcj, rocj
    
def main(rho, C, H, N, O, W, HR):
    
    typ = typex(C, H, O, W)
    product, fur = jerar(C, H, N, O)
    mw, q, w, Do = mw_hp(product, fur, W, HR)
    D, P, Q, gcj, acj, bcj, rocj = calc_D_P(Do, q, w, rho)
    return rho, mw, typ, D, P, Q/F, gcj, acj, Do, w

#.............................................................
# execution
final = main(rho, C, H, N, O, W, HR)
round_final = [round(i, 2) for i in final]
#.............................................................
#saving
import csv

header = ['name', 'rho (g/cc)', 'mw (g)', 'typ', 'D (m/s)', 'P (GPa)',\
          'Q (cal/g)', 'Gamma_cj','Jones parameter','Do (m/s)','w (km/s)']
round_final.insert(0, name) 
data = [i for i in round_final]

with open(filename + '.csv', 'a', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(header)
    writer.writerow(data)
    print ('results:')
    print(header[4:6] + [round(i,1) for i in data[4:6]])
    print("output file saved as " + filename)
    