import argparse
import numpy as np
import scipy.linalg as linalg
import math
import sys
import re
import json


aa3to1=dict({'ALA':'A', 'CYS':'C', 'ASP':'D', 'GLU':'E', 'PHE':'F', 'GLY':'G',
        'HIS':'H', 'ILE':'I', 'LYS':'K', 'LEU':'L', 'MET':'M', 'ASN':'N', 'PRO':'P',
        'GLN':'Q', 'ARG':'R', 'SER':'S', 'THR':'T', 'VAL':'V', 'TRP':'W', 'TYR':'Y'})

def pdb_opener(file_name):
    if file_name.endswith('pdb.gz'):
        return gzip.open(file_name, 'rt')
    elif file_name.endswith('pdb'):
        return open(file_name, 'r')
    else:
        print("Unknown format!!")
        exit(0)

def conformation_from_file(file_name, atoms=['N', 'CA', 'C', 'O', 'CB'], residue_major=True, dtype=np.float64):
    sequence = ''
    with pdb_opener(file_name) as infile:
        for line in infile:
            # sequence
            if line.startswith('ATOM'):
                if line[12:16] == ' CA ':
                    try:
                        sequence += aa3to1[line.split()[3]]
                    except KeyError:
                        sequence += 'X'

            # coords
            if line.startswith('ATOM') and line[12:16].strip() in atoms:
                try:
                    atom_coords = np.vstack([atom_coords, np.asarray([float(line[30:38]), float(line[38:46]), float(line[46:54])], dtype=dtype)])
                except NameError:
                    atom_coords = np.array([[line[30:38], line[38:46], line[46:54]]], dtype=dtype)

        if residue_major:
            return atom_coords.reshape(-1, len(atoms), 3), sequence
        else:
            return atom_coords, sequence

def ZN_coor(file):
    with open(file) as infile:
        for line in infile:
            if line.startswith('HETATM'):
                ZN = np.array([line[30:38].split(),line[38:46].split(),line[46:54].split()], dtype=np.float64).reshape(3)
    return ZN        
        
def HIS_coor(file,coor_aa):
    with open(file) as infile:
        for line in infile:
            if line.startswith('ATOM') and line[18:26].split()[-1] == coor_aa[0] and line[17:20] == coor_aa[1]:
                if line[13:16].split()[0] == 'ND1':
                    ND1 = np.array([line[30:38].split(),line[38:46].split(),line[46:54].split()], dtype=np.float64).reshape(3)
                elif line[13:16].split()[0] == 'NE2':
                    NE2 = np.array([line[30:38].split(),line[38:46].split(),line[46:54].split()], dtype=np.float64).reshape(3)
    return ND1,NE2

def ASP_coor(file,coor_aa):
    with open(file) as infile:
        for line in infile:
            if line.startswith('ATOM') and line[18:26].split()[-1] == coor_aa[0] and line[17:20] == coor_aa[1]:
                if line[13:16].split()[0] == 'OD1':
                    OD1 = np.array([line[30:38].split(),line[38:46].split(),line[46:54].split()], dtype=np.float64).reshape(3)
                elif line[13:16].split()[0] == 'OD2':
                    OD2 = np.array([line[30:38].split(),line[38:46].split(),line[46:54].split()], dtype=np.float64).reshape(3)
    return OD1,OD2


def HIS_coor_atom(ND,NE,ZN):
    ZN_NE_dis = np.linalg.norm(NE-ZN)
    ZN_ND_dis = np.linalg.norm(ND-ZN)
    if ZN_NE_dis >= ZN_ND_dis:
        return ['ND1','CG','CD2']
    elif ZN_NE_dis < ZN_ND_dis:
        return ['NE2','CE1','ND1']
    

def ASP_coor_atom(OD1,OD2,ZN):
    ZN_OD1_dis = np.linalg.norm(OD1-ZN)
    ZN_OD2_dis = np.linalg.norm(OD2-ZN)
    if ZN_OD1_dis >= ZN_OD2_dis:
        return ['OD2','CG','CB']
    elif ZN_OD1_dis < ZN_OD2_dis:
        return ['OD1','CG','CB']

def number_to_letter(number):
    return chr(number + 64)

def letter_to_number(letter):
    return ord(letter) - 64
    
parser = argparse.ArgumentParser()
parser.add_argument('-pdb', type=str, default=None, help='name of the input pdb file')
parser.add_argument('-pdblist', type=str, default=None, help='name of the input pdb list')
parser.add_argument('-coordination_frame', type=str, default=None, help='name of the coordination frame file')
parser.add_argument('-scaffold_pdb', type=str, default=None, help='name of the scaffold pdb file')
parser.add_argument('-symmetry_number', type=int, default=None, help='the symmetry definition of the input pdb file')
parser.add_argument('-metal_num', type=int, default=None, help='the number of zincs associated with each chain')
parser.add_argument('-output_path', type=str, default=None, help="the path of the output file")
parser.add_argument('-output_prefix', type=str, default=None, help="the name of the output file")
parser.add_argument('-match_path', type=str, default=None, help="the path of the Match output pdb")
args = parser.parse_args()

path = args.pdb
pdblist = args.pdblist
ZN_num = args.metal_num
output_path = args.output_path
match_path = args.match_path
good_pdb = []

if path != None:
    good_pdb.append(path)
elif pdblist != None:
    with open(pdblist) as f:
        for line in f:
            good_pdb.append(line.strip())

if args.output_prefix == None:
    output_name = ''
else:
    output_name = args.output_prefix
            
for ii in good_pdb:
    _,sequence = conformation_from_file(args.scaffold_pdb,atoms=['CA'])
    scaffold_pose_seq_num = len(sequence)
    symm_num = args.symmetry_number
    zn_cluster = re.findall(r'Zn_Cluster[0-9]{4}',ii)
    with open(args.coordination_frame) as f:
        json_load = json.load(f)
    zn_pair_list = []
    zn_pair_list2 = []
    for pair in range(0,len(zn_cluster)):
        if pair % 2 == 0:
            zn_pair_list.append(json_load[zn_cluster[pair]]['center_match'])
        else:
            zn_pair_list2.append(json_load[zn_cluster[pair]]['center_match'])

    zn_res_list = []
    for jj in zn_pair_list:
        zn_res = [aa[1:] for aa in re.compile(r'[HD][0-9]+').findall(jj) if aa != '']
        match_file_path = match_path+'/'+jj
        ZN = ZN_coor(match_file_path)
        
        res = []
        for resi,resn in enumerate(jj.split('_')[-2]):
            if resn == 'H':
                atom_1,atom_2 = HIS_coor(match_file_path,[zn_res[resi],'HIS'])
                res_list = HIS_coor_atom(atom_1,atom_2,ZN)
                res.append(res_list)
            elif resn == 'D':
                atom_1,atom_2 = ASP_coor(match_file_path,[zn_res[resi],'ASP'])
                res_list = ASP_coor_atom(atom_1,atom_2,ZN)
                res.append(res_list)
        zn_res_list.append([zn_res,res])        
        
    zn_res_list2 = []
    for jj in zn_pair_list2:
        zn_res = [aa[1:] for aa in re.compile(r'[HD][0-9]+').findall(jj) if aa != '']
        match_file_path = match_path+'/'+jj
        ZN = ZN_coor(match_file_path)
        
        res = []
        for resi,resn in enumerate(jj.split('_')[-2]):
            if resn == 'H':
                atom_1,atom_2 = HIS_coor(match_file_path,[zn_res[resi],'HIS'])
                res_list = HIS_coor_atom(atom_1,atom_2,ZN)
                res.append(res_list)
            elif resn == 'D':
                atom_1,atom_2 = ASP_coor(match_file_path,[zn_res[resi],'ASP'])
                res_list = ASP_coor_atom(atom_1,atom_2,ZN)
                res.append(res_list)
        zn_res_list2.append([zn_res,res])
        
    with open(ii) as infile:
        total_ZN_coor_old = []
        for line in infile:
            if line.startswith('HETATM'):
                total_ZN_coor_old.append([line[30:38].split(),line[38:46].split(),line[46:54].split()])
    total_ZN_coor_old = np.array(total_ZN_coor_old, dtype=np.float64).reshape(-1,3)
    
    total_ZN_coor = []
    if symm_num > 2:
        for c_num in range(len(total_ZN_coor_old)//symm_num):
            total_ZN_coor.append(total_ZN_coor_old[symm_num*c_num])
    elif symm_num == 2:
        total_ZN_coor = total_ZN_coor_old
        
    with open(ii) as f:
        old_lines = [line for line in f]
            
    distance_cst = []
    angle1_cst = []
    angle2_cst = []
    dihedral_cst = []
    pdb_file = []
 
    if symm_num == 2:
        pass
        
    else:
        for symm in range(0,symm_num):   
            for aa in range(0,ZN_num):
                ZN = str((scaffold_pose_seq_num+ZN_num)*(symm+1)-ZN_num+aa+1)
                res_1 = [str(int(res)+(scaffold_pose_seq_num+ZN_num)*symm) for res in zn_res_list[aa][0]]
                new_symm = (symm+symm_num-1)%symm_num
                res_2 = [str(int(res)+(scaffold_pose_seq_num+ZN_num)*new_symm) for res in zn_res_list2[aa][0]]
                for ri,ra in enumerate(res_1):
                    distance_cst.append('AtomPair  '+ zn_res_list[aa][1][ri][0] +'  '+ ra +'  ZN  '+ ZN + ' BOUNDED 1.9 2.3 0.2 zinc_coordination')
                    angle1_cst.append('Angle  '+ zn_res_list[aa][1][ri][1] +'  '+ ra +'  '+ zn_res_list[aa][1][ri][0] +'  '+ ra +'  ZN  '+ ZN +' CIRCULARHARMONIC 2.181 0.2') 
                    if ri == 0:
                        angle2_cst.append('Angle  '+ zn_res_list[aa][1][ri][0] +'  '+ ra +'  ZN  '+ ZN +'  '+ zn_res_list[aa][1][(ri+1)%len(res_1)][0] +'  '+ res_1[(ri+1)%len(res_1)] +' CIRCULARHARMONIC 1.911 0.2')
                    if zn_res_list[aa][1][ri][0] in ['NE2', 'ND1']:
                        dihedral_cst.append('Dihedral  '+ zn_res_list[aa][1][ri][2] +'  '+ ra +'  '+ zn_res_list[aa][1][ri][1] +'  '+ ra +'  '+ zn_res_list[aa][1][ri][0] +'  '+ ra +'  ZN  '+ ZN +' CIRCULARHARMONIC 3.141 0.2')
                    elif zn_res_list[aa][1][ri][0] in ['OD1', 'OD2']:
                        dihedral_cst.append('Dihedral  '+ zn_res_list[aa][1][ri][2] +'  '+ ra +'  '+ zn_res_list[aa][1][ri][1] +'  '+ ra +'  '+ zn_res_list[aa][1][ri][0] +'  '+ ra +'  ZN  '+ ZN +' CIRCULARHARMONIC 3.141 0.3')
                for ri,ra in enumerate(res_2):
                    distance_cst.append('AtomPair  '+ zn_res_list2[aa][1][ri][0] +'  '+ ra +'  ZN  '+ ZN + ' BOUNDED 1.9 2.3 0.2 zinc_coordination')
                    angle1_cst.append('Angle  '+ zn_res_list2[aa][1][ri][1] +'  '+ ra +'  '+ zn_res_list2[aa][1][ri][0] +'  '+ ra +'  ZN  '+ ZN +' CIRCULARHARMONIC 2.181 0.2')
                    if ri == 0:
                        angle2_cst.append('Angle  '+ zn_res_list2[aa][1][ri][0] +'  '+ ra +'  ZN  '+ ZN +'  '+ zn_res_list2[aa][1][(ri+1)%len(res_2)][0] +'  '+ res_2[(ri+1)%len(res_2)] +' CIRCULARHARMONIC 1.911 0.2')
                    if zn_res_list2[aa][1][ri][0] in ['NE2', 'ND1']:
                        dihedral_cst.append('Dihedral  '+ zn_res_list2[aa][1][ri][2] +'  '+ ra +'  '+ zn_res_list2[aa][1][ri][1] +'  '+ ra +'  '+ zn_res_list2[aa][1][ri][0] +'  '+ ra +'  ZN  '+ ZN +' CIRCULARHARMONIC 3.141 0.2')
                    elif zn_res_list2[aa][1][ri][0] in ['OD1', 'OD2']:
                        dihedral_cst.append('Dihedral  '+ zn_res_list2[aa][1][ri][2] +'  '+ ra +'  '+ zn_res_list2[aa][1][ri][1] +'  '+ ra +'  '+ zn_res_list2[aa][1][ri][0] +'  '+ ra +'  ZN  '+ ZN +' CIRCULARHARMONIC 3.141 0.3')
                for ri in range(len(res_1)):
                    angle2_cst.append('Angle  '+ zn_res_list[aa][1][ri][0] +'  '+ res_1[ri] +'  ZN  '+ ZN +'  '+ zn_res_list2[aa][1][ri][0] +'  '+ res_2[ri] +' CIRCULARHARMONIC 1.911 0.2')
                        
        zn_res_list = zn_res_list + zn_res_list2                
        for each in zn_res_list:
            for jj in each[0]:
                pdb_file.append('REMARK PDBinfo-LABEL:   '+str(jj)+' Zn_res')
 
        chain_name = []
        for num in range(ZN_num):
            chain_name.append(number_to_letter(symm_num + num + 1)) 

        new_lines = []   
        for line in old_lines:
            if line.startswith('ATOM') and line[21] == 'A':
                new_lines.append(line)

            if line.startswith('HETATM') and line[17:20].strip() == 'ZNX' and line[21] in chain_name:
                line = line[:21]+str(number_to_letter(letter_to_number(line[21])-symm_num+1))+line[22:]
                new_lines.append(line)


    cst_list = distance_cst + angle1_cst + angle2_cst + dihedral_cst
    with open(output_path+'/'+output_name+ii.split('/')[-1].split('.pdb')[0]+'.cst','w') as w:
        w.writelines('\n'.join(cst_list))

    with open(output_path+'/'+output_name+ii.split('/')[-1].split('.pdb')[0]+'.pdb','w') as a:
        for line in new_lines:
            a.write(line)
        a.writelines('\n'.join(pdb_file))          

        
 


