#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 09:08:15 2017

@author: nienke
"""


import rdkit
import pandas as pd
import os
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, SaltRemover
from molvs import tautomer
from rdkit.Chem.rdchem import BondDir, BondStereo, BondType


#%%
dir_RT='/users/nienkemoret/Dropbox (HMS)/CBDM-SORGER/Collaborations/LSP_data_organization/ChemInformatics_Files/Reagenttracker'
dir_chembl='/users/nienkemoret/Dropbox (HMS)/CBDM-SORGER/Collaborations/LSP_data_organization/ChemInformatics_Files/ChEMBL_24_1'

print(os.listdir(dir_RT))
print(os.listdir(dir_chembl))

#%%
os.chdir(dir_RT)
molecules_RT=pd.read_csv('smallmol_list_RT_20190617.csv')


#%%
def make_tautomer(mol):
    return(tautomer.TautomerEnumerator(mol))
    
    #%%
    enum=tautomer.TautomerEnumerator(max_tautomers=10)

#%%

tautomers_RT=[]
#smiles_tautomers=[]

for rec in molecules_RT.iterrows():
    try:
        mol=Chem.MolFromSmiles(rec[1][3])
        mol1=enum(mol)
        a=[rec[1][0],mol1]
        b=[(a[0], Chem.MolToSmiles(s)) for s in a[1]]
        #b=[rec[1][0],Chem.MolToSmiles(mol1)]#make a loop in a loop
        c=pd.DataFrame(b)
        tautomers_RT.append(c)
        print(rec[1][0])
    except Exception:
        continue 
    
    #%%
df_tautomers1=pd.concat(tautomers_RT)
    #%%
    os.chdir(dir_chembl)
    df_tautomers1.columns=['hms_id','smiles']
    df_tautomers1.to_csv('smallmol_list_RT_20190617_tautomers.csv',index=False)
    
    
    #%%
    
   # [(a[0], s) for s in a[1]]
   #[(tautomers_RT[0][0], s) for s in tautomers_RT[0][1]]    
    #%%
    mol1=Chem.MolFromSmiles('C[C@@H](C1=C(C(=O)C2=C(O1)C=CC(=C2)F)C3=CC(=CC=C3)F)N4C5=C(C(=N4)C6=CC(=C(C=C6)OC(C)C)F)C(=NC=N5)N')
    
    enum=tautomer.TautomerEnumerator()
    
    mol2=enum(mol1)
    
    print(mol1)
    print(mol2)