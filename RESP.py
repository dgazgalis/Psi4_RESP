"""
Based on the make it rain parital charges sciprt for smiles to parital charges
https://github.com/pablo-arantes/making-it-rain/blob/main/Partial_Charges.ipynb

Meant for use in conda enviroment RESP
RESP can be made from RESP.yml
conda create -n RESP --file RESP_env.txt
or 
conda env create -f RESP.yml
conda activate RESP
Output mol2 can be used for follow up calcs in a CADD package or amber
You should change the functional (defined as method) and the basis set (defined as basisSet) to what ever you need for your purposes
"""

import sys
import os
import psi4
import resp
#for d3 corrected dispersion, the package conda-forge::dftd3-python is included in the enviroment for this reason
#import dftd3
from openbabel import openbabel as ob
from rdkit import Chem
from rdkit.Chem import AllChem

#Name of file without the .sdf extentsion
molId = "mol"
#sys.argv can be used to define mol from the command line 
#Uncomment line 29 to enable
#molId = sys.argv[0]

#Parameters to optimize a small molocule and gnerate RESP charges
method = "WB97X" #Suggested functionals "WB97X-D3","M06-2X-D3", "B3LYP-D3"  range seperate dispersion corrected though this requires grimme's DFTD3 package
basisSet = "6-311+G" # Suggested basis sets for good cost to preformance with diffuse and polarizable functions "aug-cc-pVTZ", "def2-TZVPPD", "6-311++G(2d,p)" ]

#Psi4 sspecific parameters
psi4.set_num_threads(32) #number of threads to run each psi4 calc on 
psi4.set_memory('32 GB') #memory for calculation for psi4

def calcRESPCharges(mol, basisSet, method, gridPsi4 = 1):
    options = {'BASIS_ESP': basisSet, 'METHOD_ESP': method, 'RESP_A': 0.0005, 'RESP_B': 0.1, 'VDW_SCALE_FACTORS':[1.4, 1.6, 1.8, 2.0], 'VDW_POINT_DENSITY':int(gridPsi4)}
    resp_charges = resp.resp([mol], options)
    return resp_charges

inputFile = "./"+molId+".sdf"
#Convert SDF to MOL2 to preserve cordinates when reading in compound for tleap
obConversion = ob.OBConversion()
obConversion.SetInAndOutFormats("sdf", "mol2")
obmol = ob.OBMol()
obConversion.ReadFile(obmol, inputFile)
obConversion.WriteFile(obmol, molId+".mol2")

obConversion.SetInAndOutFormats("xyz", "mol2")

print('Trying:', molId)
suppl = Chem.SDMolSupplier(inputFile, removeHs=False) 
mol = next(suppl)

#Quick optimization using MMFF to put the inital mol into more sane cordinates for the more costly quantum mechanics based optimization  
#NOTE: NEED TO TEST
#AllChem.MMFFOptimizeMolecule(mol)

num_atoms = mol.GetNumAtoms()
#Assuming netural charge and singlet
#Psi4 recognizes the first "<int> <int>" pair as the charge and multiplicity, respectively
xyz_string='0 1 \n'

for counter in range(num_atoms):
    pos=mol.GetConformer().GetAtomPosition(counter)
    xyz_string = xyz_string + ("%s %12.6f %12.6f %12.6f\n" % (mol.GetAtomWithIdx(counter).GetSymbol(), pos.x, pos.y, pos.z) )

psi_mol = psi4.geometry(xyz_string)
outfile_mol2 = inputFile[:-4]+".mol2"
        
print('Running geometry optimization...')
methodNbasisSet = method+"/"+basisSet
#psi4.optimize optimizes geometry of the molocule which the the most costly step in this script
#Max iterations is set to 500. If the scf does not converge, try a 2 step protocol for QM or force field based minimisation first 
#One step
psi4.optimize(methodNbasisSet, molecule=psi_mol, maxiter=500)

#Uncomment for two step protocol and comment out line 73
"""
inital_methodNbasisSet = "B3LYP/6-311+G"
psi4.optimize(inital_methodNbasisSet, molecule=psi_mol, maxiter=500)
psi4.optimize(methodNbasisSet, molecule=psi_mol, maxiter=500)
"""
#Charge calcuation step, fitting might take some time for larger molocules 
resp_charges = calcRESPCharges(psi_mol, basisSet, method, gridPsi4 = 1)

### save coords to xyz file
psi4out_xyz = molId + '.xyz'
psi_mol.save_xyz_file(psi4out_xyz,1)

### read xyz file and write as mol2
ob_mol = ob.OBMol()
obConversion.ReadFile(ob_mol, psi4out_xyz)

### write as mol2
outfile_mol2 = "./"+molId+"_partialChgs.mol2"
obConversion.WriteFile(ob_mol, outfile_mol2)

### set new partial charges
count = 0
newChg_temp = resp_charges[1]
print("RESP Charges: ", newChg_temp)
for atom in ob.OBMolAtomIter(ob_mol):
    newChg = newChg_temp[count]
    atom.SetPartialCharge(newChg)
    count += 1

### write as mol2
outfile_mol2 = "./"+molId+"_RESP.mol2"
print("Finished. Saved compound with partial charges as mol2 file: %s" % outfile_mol2)
obConversion.WriteFile(ob_mol, outfile_mol2)