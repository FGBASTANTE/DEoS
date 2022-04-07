# DEoS
PREDICTIVE MODEL FOR EXPLOSIVE DETONATION PARAMETERS

In this script is implemented the DEoS model, a predictive model of explosive
detonation velocity and pressure based on first order
approximation of the detonation velocity equation. Please see:
    
Predictive model of explosive detonation parameters from an equation of state
based on detonation velocity. Fernando G. Bastante*, María Araújo* and Eduardo
Giráldez*. [DOI: 10.1039/D2CP00085G Phys. Chem. Chem. Phys., 2022, 24, 8189-8195]
(https://doi.org/10.1039/D2CP00085G)

*CINTECX, GESSMin Group, Department of Natural Resources
and Environmental Engineering (EEME), University of Vigo, Pontevedra, Spain

Data Inputs (example in the script):
-    name = "PBX9502" # the name of explosive without blanks
-    rho = 1.90  # density in g/cc (must be >1g//cc)
-    C = 2.30    # number of moles
-    H = 2.23    # number of moles
-    N = 2.21    # number of moles
-    O = 2.21    # number of moles 
-    W = 3.81    # weight percent of other elements (Cl/F/P/Si % must be <10)
-    HR = -205.5 # enthalpy standard of formation in cal/g (1 atm, 298 K)
-    filename = "PBX9502_1.90" # filename for saving results

The script can read datainputs from cli...:
 -   python deos.py name rho C H N O W HR outputfilename

example:
 -   python deos.py PBX9502 1.90 2.30 2.23 2.21 2.21 3.81 -205.5 out_pbx9052 

the units are:
-        ........  string, g/cc, mol, mol, mol, mol, %weight, cal/g, string

...or you can introduce them in the script
