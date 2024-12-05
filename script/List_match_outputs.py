import os
import argparse
import numpy as np
import glob

def HIS_coor(file,coor_aa):
    with open(file) as infile:
        His_coor = []
        for line in infile:
            if line.startswith('ATOM') and line[18:26].split()[-1] == coor_aa[0] and line[17:20] == coor_aa[1]:
                if line[13:15].split()[0] == 'ND':
                    ND = np.array([line[30:38].split(),line[38:46].split(),line[46:54].split()], dtype=np.float64).reshape(3)
                    His_coor.append(['ND',ND])
                elif line[13:15].split()[0] == 'NE':
                    NE = np.array([line[30:38].split(),line[38:46].split(),line[46:54].split()], dtype=np.float64).reshape(3)
                    His_coor.append(['NE',NE])

    return His_coor

def ASP_coor(file,coor_aa):
    with open(file) as infile:
        ASP_coor = []
        for line in infile:
            if line.startswith('ATOM') and line[18:26].split()[-1] == coor_aa[0] and line[17:20] == coor_aa[1]:
                if line[13:16].split()[0] == 'OD1':
                    OD1 = np.array([line[30:38].split(),line[38:46].split(),line[46:54].split()], dtype=np.float64).reshape(3)
                    ASP_coor.append(['OD1',OD1])
                elif line[13:16].split()[0] == 'OD2':
                    OD2 = np.array([line[30:38].split(),line[38:46].split(),line[46:54].split()], dtype=np.float64).reshape(3)
                    ASP_coor.append(['OD2',OD2])
                    
    return ASP_coor

def HIS_coor_atom(HIS_coor,ZN):
    ND = HIS_coor[0][1]
    NE = HIS_coor[1][1]
    ZN_NE_dis = np.linalg.norm(NE-ZN)
    ZN_ND_dis = np.linalg.norm(ND-ZN)
    if ZN_NE_dis >= ZN_ND_dis:
        return HIS_coor[0]
    elif ZN_NE_dis < ZN_ND_dis:
        return HIS_coor[1]


def ASP_coor_atom(ASP_coor,ZN):
    OD1 = ASP_coor[0][1]
    OD2 = ASP_coor[1][1]
    ZN_OD1_dis = np.linalg.norm(OD1-ZN)
    ZN_OD2_dis = np.linalg.norm(OD2-ZN)

    if ZN_OD1_dis >= ZN_OD2_dis:
        return ASP_coor[1]
    elif ZN_OD1_dis < ZN_OD2_dis:
        return ASP_coor[0]


parser = argparse.ArgumentParser()
parser.add_argument('-input_list', type=str, default=None, help="the list of the Match output pdb")
parser.add_argument('-output_path', type=str, default=None, help="the path of the output file")
parser.add_argument('-output_name', type=str, default=None, help="the name of the output file")
args = parser.parse_args()

if args.output_name == None:
    output_name = 'output'
else:
    output_name = args.output_name

with open(args.input_list) as f:
    pdb_list = [line.strip() for line in f]

    
lines = []
first_line = 'description residue_type index atom1 atom1_x atom1_y atom1_z atom2 atom2_x atom2_y atom2_z ZN ZN_x ZN_y ZN_z\n'
lines.append(first_line)
for ii in pdb_list:
    if ii[-3:] != 'pdb':
        continue
    pdb_name = ii.split('/')[-1]
    coor_res = pdb_name.split('_')[3]
    index = pdb_name.split('_')[4]
    with open(ii) as infile:
        file_coordinate_aa = []
        for line in infile:
            if line.startswith('REMARK') and line[51:54] in ['HIS','ASP']:
                file_coordinate_aa.append([line[54:59].split()[-1],line[51:54]])
            if line.startswith('HETATM'):
                ZN = np.array([line[30:38].split(),line[38:46].split(),line[46:54].split()], dtype=np.float64).reshape(3)


    coor_aa1 = file_coordinate_aa[0]
    if coor_aa1[1] == 'HIS':
        His_coor1 = HIS_coor(ii,coor_aa1)       
    elif coor_aa1[1] == 'ASP':
        ASP_coor1 = ASP_coor(ii,coor_aa1)

    coor_aa2 = file_coordinate_aa[1]
    if coor_aa2[1] == 'HIS':
        His_coor2 = HIS_coor(ii,coor_aa2)
    elif coor_aa2[1] == 'ASP':
        ASP_coor2 = ASP_coor(ii,coor_aa2)

    if coor_aa1[1] == 'HIS':
        His_coor1_ZN = HIS_coor_atom(His_coor1,ZN)
        atom1 = His_coor1_ZN[0]
        atom1_x = His_coor1_ZN[1][0]
        atom1_y = His_coor1_ZN[1][1]
        atom1_z = His_coor1_ZN[1][2]
    elif coor_aa1[1] == 'ASP':
        ASP_coor1_ZN = ASP_coor_atom(ASP_coor1,ZN)
        atom1 = ASP_coor1_ZN[0]
        atom1_x = ASP_coor1_ZN[1][0]
        atom1_y = ASP_coor1_ZN[1][1]
        atom1_z = ASP_coor1_ZN[1][2]
    if coor_aa2[1] == 'HIS':
        His_coor2_ZN = HIS_coor_atom(His_coor2,ZN)
        atom2 = His_coor2_ZN[0]
        atom2_x = His_coor2_ZN[1][0]
        atom2_y = His_coor2_ZN[1][1]
        atom2_z = His_coor2_ZN[1][2]
    elif coor_aa2[1] == 'ASP':
        ASP_coor2_ZN = ASP_coor_atom(ASP_coor2,ZN)
        atom2 = ASP_coor2_ZN[0]
        atom2_x = ASP_coor2_ZN[1][0]
        atom2_y = ASP_coor2_ZN[1][1]
        atom2_z = ASP_coor2_ZN[1][2]
    ZN_x = ZN[0]
    ZN_y = ZN[1]
    ZN_z = ZN[2]
    line = pdb_name+' '+coor_res+' '+index+' '+atom1+' '+str(atom1_x)+' '+str(atom1_y)+' '\
           +str(atom1_z)+' '+atom2+' '+str(atom2_x)+' '+str(atom2_y)+' '+str(atom2_z)\
           +' ZN '+str(ZN_x)+' '+str(ZN_y)+' '+str(ZN_z)+'\n'
    lines.append(line)

w = open(args.output_path+'/'+output_name+'.list','w')
w.writelines(lines)
w.close()
