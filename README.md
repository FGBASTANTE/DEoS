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

Some results
+500 data

![image](https://user-images.githubusercontent.com/52360383/164480516-2752e152-155b-45fd-9b8d-4a0275497431.png)


+250 data

![image](https://user-images.githubusercontent.com/52360383/164482623-50d5769a-f7b4-4eea-8ac8-26f7d34cc66d.png)

NQ

![image](https://user-images.githubusercontent.com/52360383/164490428-9ddc84a0-23ed-4103-94a3-ef0c6cdf9ba5.png)

![image](https://user-images.githubusercontent.com/52360383/164486567-43b7c7fb-2924-4343-b2c1-4f4661ad1df7.png)

COMPOSITES RDX/TNT

![image](https://user-images.githubusercontent.com/52360383/164491432-3120a366-45ba-4e74-91c8-a1d1a4d88cc2.png)

![image](https://user-images.githubusercontent.com/52360383/164489150-8486c327-5bb3-41b2-a8ef-44145a1f7553.png)



