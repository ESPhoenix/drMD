## BASIC LIBS
import os
from os import path as p
import pandas as pd

########################## write pdb dataframe to pdb file
def df2Pdb(df, outFile,
           chain=True):
    with open(outFile,"w") as f:
        for _, row in df.iterrows():
            pdbLine = f"{row['ATOM']:<6}"
            pdbLine += f"{row['ATOM_ID']:>5}{' '*2}"
            pdbLine += f"{row['ATOM_NAME']:<4}"
            pdbLine += f"{row['RES_NAME']:<4}"
            try: 
                pdbLine += f"{row['CHAIN_ID']:<1}{' '*1}"
            except:
                pdbLine += ' '*2
            pdbLine += f"{row['RES_ID']:<7}"
            pdbLine += f"{row['X']:>8.3f}"
            pdbLine += f"{row['Y']:>8.3f}"
            pdbLine += f"{row['Z']:>8.3f}"
            pdbLine += f"{row['OCCUPANCY']:>6.2f}"
            try:
                pdbLine += f"{row['BETAFACTOR']:>6.2f}"
            except:
                pdbLine += row["BETAFACTOR"]
            try:
                pdbLine += f"{row['ELEMENT']:>12}\n"
            except:
                element = row["ATOM_NAME"][0]
                pdbLine += f"{element:>12}\n"
            f.write(pdbLine)
############################### read pdb file as pdb dataframe
def pdb2df(protPdb):
    columns = ['ATOM', 'ATOM_ID', 'ATOM_NAME', 'RES_NAME', 'CHAIN_ID', 'RES_ID', 'X', 'Y', 'Z', 'OCCUPANCY', 'BETAFACTOR', 'ELEMENT']
    data = []
    with open(protPdb, 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atom_type = line[0:6].strip()
                atom_id = int(line[6:11].strip())
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain_id = line[21:22].strip()
                if chain_id == '':
                    chain_id = None
                res_id = int(line[22:26].strip())
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                occupancy = float(line[54:60].strip())
                try:
                    temp_factor = float(line[60:66].strip())
                except:
                    temp_factor = " "*6
                element = line[76:78].strip()

                data.append([atom_type, atom_id, atom_name, res_name, chain_id, res_id, x, y, z, occupancy, temp_factor, element])

    return pd.DataFrame(data, columns=columns)
############################### read list of pdb files as dataframes, concatonate, then write back to pdb file
def mergePdbs(pdbList,outFile):
    dfList=[]
    for pdbFile in pdbList:
        df = pdb2df(pdbFile)
        dfList.append(df)
    mergedDf = pd.concat(dfList,axis=0)
    df2Pdb(df=mergedDf, outFile=outFile, chain=True)

############################### apply a a bunch of fixes to a pdb dataframe
def fix_atom_names(df): 
    # deal with unwanted apostrophies (prime)
    df.loc[:,'ATOM_NAME'] = df['ATOM_NAME'].str.replace("'", "")
    # deal with numbers at the beginnig of atom names
    df.loc[:,'ATOM_NAME'] = df['ATOM_NAME'].replace(r'^(\d+)(.+)$', r'\2\1', regex=True)
    # deal with "A" at the start of atom name
    df.loc[:,'ATOM_NAME'] = df['ATOM_NAME'].apply(lambda x: x.lstrip('A') if x.startswith('A') else x)

    ## ensure unique names
    count_series = df.groupby('ATOM_NAME').cumcount()
    df.loc[:,'ATOM_NAME'] = df['ATOM_NAME'] + "_" +count_series.astype(str)
    df.loc[:,'ATOM_NAME'] = df['ATOM_NAME'].str.replace("_0", "")
    df.loc[:,'ATOM_NAME'] = df['ATOM_NAME'].str.replace("_", "")

    return df 


