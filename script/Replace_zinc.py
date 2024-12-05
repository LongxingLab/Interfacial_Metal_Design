import argparse
import numpy as np
import scipy.linalg as linalg
import sys
import re
import json
from pyrosetta import *
from pyrosetta.rosetta import *
sys.path.append('script/')
import npose_util as nu



def polyA( pdb ):
    init()
    pose = pose_from_file( pdb )
    rts = rosetta.core.chemical.ChemicalManager.get_instance().residue_type_set("fa_standard")
    rtype = rts.name_map( 'ALA' )
    res = rosetta.core.conformation.ResidueFactory.create_residue( rtype )
    for ir in range(1, pose.size() + 1):
        name = pose.residue(ir).name3()
        if name in ['GLY', 'PRO', 'CYS']:
            continue
        pose.replace_residue( ir, res, True )
    return pose

def HIS_coor(file,coor_aa):
    with open(file) as infile:
        for line in infile:
            if line.startswith('ATOM') and line[18:26].split()[-1] == coor_aa[0] and line[17:20] == coor_aa[1]:
                if line[13:15].split()[0] == 'CG':
                    CG = np.array([line[30:38].split(),line[38:46].split(),line[46:54].split()], dtype=np.float64).reshape(3)
                elif line[13:15].split()[0] == 'ND':
                    ND = np.array([line[30:38].split(),line[38:46].split(),line[46:54].split()], dtype=np.float64).reshape(3)
                elif line[13:15].split()[0] == 'CD':
                    CD = np.array([line[30:38].split(),line[38:46].split(),line[46:54].split()], dtype=np.float64).reshape(3)
                elif line[13:15].split()[0] == 'NE':
                    NE = np.array([line[30:38].split(),line[38:46].split(),line[46:54].split()], dtype=np.float64).reshape(3)
                elif line[13:15].split()[0] == 'CE':
                    CE = np.array([line[30:38].split(),line[38:46].split(),line[46:54].split()], dtype=np.float64).reshape(3)

    return CG,ND,CD,NE,CE

def ASP_coor(file,coor_aa):
    with open(file) as infile:
        for line in infile:
            if line.startswith('ATOM') and line[18:26].split()[-1] == coor_aa[0] and line[17:20] == coor_aa[1]:
                if line[13:16].split()[0] == 'OD1':
                    OD1 = np.array([line[30:38].split(),line[38:46].split(),line[46:54].split()], dtype=np.float64).reshape(3)
                elif line[13:16].split()[0] == 'OD2':
                    OD2 = np.array([line[30:38].split(),line[38:46].split(),line[46:54].split()], dtype=np.float64).reshape(3)
                elif line[13:16].split()[0] == 'CG':
                    CG = np.array([line[30:38].split(),line[38:46].split(),line[46:54].split()], dtype=np.float64).reshape(3)
                elif line[13:16].split()[0] == 'CB':
                    CB = np.array([line[30:38].split(),line[38:46].split(),line[46:54].split()], dtype=np.float64).reshape(3)
    return OD1,OD2,CG,CB

def HIS_coor_atom(CG,ND,CD,NE,CE,ZN):
    ZN_NE_dis = np.linalg.norm(NE-ZN)
    ZN_ND_dis = np.linalg.norm(ND-ZN)
    if ZN_NE_dis >= ZN_ND_dis:
        ZN_dis = ZN_ND_dis
        ZN_U1 = ND
        ZN_U2 = CG
        ZN_U3 = CD
        ZN_U4 = CE
    elif ZN_NE_dis < ZN_ND_dis:
        ZN_dis = ZN_NE_dis
        ZN_U1 = NE
        ZN_U2 = CE
        ZN_U3 = ND
        ZN_U4 = CD
    return ZN_dis,ZN_U1,ZN_U2,ZN_U3,ZN_U4

def ASP_coor_atom(OD1,OD2,CG,CB,ZN):
    ZN_OD1_dis = np.linalg.norm(OD1-ZN)
    ZN_OD2_dis = np.linalg.norm(OD2-ZN)
    ZN_U2 = CG
    ZN_U3 = CB
    if ZN_OD1_dis >= ZN_OD2_dis:
        ZN_dis = ZN_OD2_dis
        ZN_U1 = OD2
        ZN_U4 = OD1
    elif ZN_OD1_dis < ZN_OD2_dis:
        ZN_dis = ZN_OD1_dis
        ZN_U1 = OD1
        ZN_U4 = OD2
    return ZN_dis,ZN_U1,ZN_U2,ZN_U3,ZN_U4

def ang(U2,U1,D1):
    U2U1 = U2 - U1
    D1U1 = D1 - U1
    U2U1_norm = np.linalg.norm(U2U1)
    D1U1_norm = np.linalg.norm(D1U1)
    U2U1_dot_D1U1 = U2U1.dot(D1U1)
    cos = np.arccos(U2U1_dot_D1U1/(U2U1_norm * D1U1_norm))
    ang_U2D1 = np.rad2deg(cos)
    return ang_U2D1

def tor(a, b, c, d):
    b0 = -1.0*(b-a)
    b1 = c - b
    b2 = d - c
    b1 = b1/np.linalg.norm(b1, axis=-1)

    v = b0 - np.sum(b0*b1, axis=-1)*b1
    w = b2 - np.sum(b2*b1, axis=-1)*b1

    x = np.sum(v*w, axis=-1)
    y = np.sum(np.cross(b1, v)*w, axis=-1)
    angle = np.rad2deg(np.arctan2(y, x))
    if angle < 0:
        angle = -angle
    return angle

def pair_zn_coorsdinate(pdbfile):
    with open(pdbfile) as infile:
        file_coordinate_aa = []
        for line in infile:
            if line.startswith('REMARK') and line[51:54] in ['HIS','ASP']:
                file_coordinate_aa.append([line[54:59].split()[-1],line[51:54]])
            if line.startswith('HETATM'):
                ZN = np.array([line[30:38].split(),line[38:46].split(),line[46:54].split()], dtype=np.float64).reshape(3)
    return file_coordinate_aa,ZN

def replace_res_HD(file_coordinate_aa,pairfile,wefold_pose_for_replace,scaffold_seq_num,symm_num): 
    HIS1_num = int(file_coordinate_aa[0][0])
    ASP2_num = int(file_coordinate_aa[1][0])
    pose_zinc = pose_from_file(pairfile)
    HIS_coor_zinc1 = pose_zinc.residue(1)
    ASP_coor_zinc2 = pose_zinc.residue(2)
    HIS_replace = pose_from_sequence('AX['+str(pose_zinc.residue_type_ptr(1)).split(':')[0]+']A').residue(2)
    ASP_replace = pose_from_sequence('ADA').residue(2)
    wefold_pose_for_replace.replace_residue(HIS1_num+scaffold_seq_num*symm_num,HIS_replace,True)
    for jj in range(1, HIS_coor_zinc1.nchi()+1):
        wefold_pose_for_replace.set_chi(jj,HIS1_num+scaffold_seq_num*symm_num,HIS_coor_zinc1.chi(jj))
    wefold_pose_for_replace.replace_residue(ASP2_num+scaffold_seq_num*symm_num,ASP_replace,True)
    for jj in range(1, ASP_coor_zinc2.nchi()+1):
        wefold_pose_for_replace.set_chi(jj,ASP2_num+scaffold_seq_num*symm_num,ASP_coor_zinc2.chi(jj))
    return wefold_pose_for_replace

def replace_res_HH(file_coordinate_aa,pairfile,wefold_pose_for_replace,scaffold_seq_num,symm_num): 
    HIS1_num = int(file_coordinate_aa[0][0])
    HIS2_num = int(file_coordinate_aa[1][0])
    pose_zinc = pose_from_file(pairfile)
    HIS_coor_zinc1 = pose_zinc.residue(1)
    HIS_coor_zinc2 = pose_zinc.residue(2)
    HIS_replace1 = pose_from_sequence('AX['+str(pose_zinc.residue_type_ptr(1)).split(':')[0]+']A').residue(2)
    HIS_replace2 = pose_from_sequence('AX['+str(pose_zinc.residue_type_ptr(2)).split(':')[0]+']A').residue(2)
    wefold_pose_for_replace.replace_residue(HIS1_num+scaffold_seq_num*symm_num,HIS_replace1,True)
    for jj in range(1, HIS_coor_zinc1.nchi()+1):
        wefold_pose_for_replace.set_chi(jj,HIS1_num+scaffold_seq_num*symm_num,HIS_coor_zinc1.chi(jj))
    wefold_pose_for_replace.replace_residue(HIS2_num+scaffold_seq_num*symm_num,HIS_replace2,True)
    for jj in range(1, HIS_coor_zinc2.nchi()+1):
        wefold_pose_for_replace.set_chi(jj,HIS2_num+scaffold_seq_num*symm_num,HIS_coor_zinc2.chi(jj))
    return wefold_pose_for_replace

def wefold_bb_stubs(coords):
    assert(coords.ndim == 3) # resi : atoms : xyz
    n = coords[:, 0, :3]
    ca = coords[:, 1, :3]
    c = coords[:, 2, :3]

    stub = np.zeros((len(n),4,4),dtype=coords.dtype)
    stub[:,3,3] = 1
    e1 = c[:, :3] - ca[:, :3]
    e1 /= np.linalg.norm(e1, axis=1)[:, None]
    e3 = np.cross(e1, n[:, :3] - ca[:, :3])
    e3 /= np.linalg.norm(e3, axis=1)[:, None]
    e2 = np.cross(e3, e1)
    e2 /= np.linalg.norm(e2, axis=1)[:, None]
    stub[:, :3, 0] = e1
    stub[:, :3, 1] = e2
    stub[:, :3, 2] = e3
    stub[:, :3, 3] = ca[:, :3]

    assert np.allclose(np.linalg.det(stub), 1)

    return stub

def create_stubs(pdb):
    npose = nu.npose_from_file(pdb)
    pose =  np.array(nu.extract_atoms(npose,[nu.N,nu.CA,nu.C])[:,:3],dtype=np.float32).reshape([-1,3,3])
    stubs = wefold_bb_stubs(pose)
    return stubs


parser = argparse.ArgumentParser()
parser.add_argument('-pdb', type=str, default=None, help='name of the input pdb file')
parser.add_argument('-pdblist', type=str, default=None, help='name of the input pdb list')
parser.add_argument('-scaffold_pdb', type=str, default=None, help='name of the scaffold pdb file')
parser.add_argument('-metal_pdb', type=str, default=None, help='name of the metal pdb file')
parser.add_argument('-metal_params', type=str, default=None, help='name of the metal params file')
parser.add_argument('-match_path', type=str, default=None, help="the path of the Match output pdb")
parser.add_argument('-coordination_frames', type=str, default=None, help='name of the coordination frame file')
parser.add_argument('-output_path', type=str, default=None, help="the path of the output file")
parser.add_argument('-output_prefix', type=str, default=None, help="the name of the output file")
args = parser.parse_args()
init('-extra_res_fa '+args.metal_params+' -ignore_unrecognized_res -mute all')

path = args.pdb
pdblist = args.pdblist

if args.output_prefix == None:
    output_name = ''
else:
    output_name = args.output_prefix
    
pdb_list = []

if path != None:
    pdb_list.append(path)
elif pdblist != None:
    with open(pdblist) as f:
        for line in f:
            pdb_list.append(line.strip())

for ii in pdb_list:
    scaffold_name = ii.split('/')[-1].split('___')[0]
    scaffold_pose = pose_from_file(args.scaffold_pdb)
    scaffold_pose_seq_num = len(scaffold_pose.sequence())
    wefold_pose = pose_from_file(ii)
    wefold_pose_seq_num = len(wefold_pose.sequence())
    symm = int(wefold_pose_seq_num/scaffold_pose_seq_num)
    xyz = pyrosetta.rosetta.numeric.xyzVector_double_t()
    zn_pose = pose_from_file(args.metal_pdb)
    zn_cluster = re.findall(r'Zn_Cluster[0-9]{4}',ii)
    zn_pair = []
    for i in range(0,int(len(zn_cluster)/2)):
        zn_pair.append([zn_cluster[i*2],zn_cluster[i*2+1]])        
    with open(args.coordination_frames) as f:
        json_load = json.load(f)

    zn_score_good = []
    for pair in zn_pair:
        zn_score_good.append([json_load[pair[0]]['center_match'],json_load[pair[1]]['center_match']])
    
    pdbfile = args.match_path
    for cluster_pair in zn_score_good:
        for one_side in cluster_pair:
            pair_coor , _ = pair_zn_coorsdinate(pdbfile+one_side)
            residue_type = [pair_coor[0][1],pair_coor[1][1]]
            residue_type.sort()
            if residue_type[0] == 'ASP' and  residue_type[1] == 'HIS':
                for symm_num in range(0,symm):
                    wefold_pose = replace_res_HD(pair_coor,pdbfile+one_side,wefold_pose,scaffold_pose_seq_num,symm_num) 
            elif residue_type[0] == 'HIS' and  residue_type[1] == 'HIS':
                for symm_num in range(0,symm):
                    wefold_pose = replace_res_HH(pair_coor,pdbfile+one_side,wefold_pose,scaffold_pose_seq_num,symm_num)

                
    new_pose_name = args.output_path+'/'+output_name+ii.split('/')[-1].split('.pdb')[0]+'.pdb'
    wefold_pose.dump_pdb(new_pose_name)    

    if symm == 2:
        for ii in zn_score_good:
            wefold_stubs = create_stubs(new_pose_name)
            zn_num1 = [int(aa) for aa in re.split(r'[HD]',re.search(r'[HD][0-9]+[HD][0-9]+',ii[0]).group()) if aa != ''][0]
            zn_num2 = [int(aa) for aa in re.split(r'[HD]',re.search(r'[HD][0-9]+[HD][0-9]+',ii[1]).group()) if aa != ''][0]
            zn_stubs_1 = create_stubs(pdbfile+ii[0])
            zn_stubs_2 = create_stubs(pdbfile+ii[1])
            _ , zn1 = pair_zn_coorsdinate(pdbfile+ii[0])
            _ , zn2 = pair_zn_coorsdinate(pdbfile+ii[1])
            zn_coor1 = np.hstack((np.array(zn1)[None, :],np.ones([1])[None, :]))
            zn_coor2 = np.hstack((np.array(zn2)[None, :],np.ones([1])[None, :]))
            zn_coor_new_1 = wefold_stubs[zn_num1-1]@ (np.linalg.inv(zn_stubs_1[0])@zn_coor1.T)
            zn_coor_new_2 = wefold_stubs[scaffold_pose_seq_num+zn_num2-1]@ (np.linalg.inv(zn_stubs_2[0])@zn_coor2.T)
            zn_coor_new = ((zn_coor_new_2+zn_coor_new_1)/2).reshape((-1,4))[0][:3]
            wefold_pose.append_pose_by_jump(zn_pose,1)
            wefold_pose.set_xyz(AtomID(wefold_pose.residue(len(wefold_pose.sequence())).atom_index("ZN"), len(wefold_pose.sequence())),xyz.assign(zn_coor_new[0],zn_coor_new[1],zn_coor_new[2]))

    else:
        for pair in zn_score_good:
            for symm_num in range(0,symm):
                wefold_stubs = create_stubs(new_pose_name)
                zn_num1 = [int(aa) for aa in re.split(r'[HD]',re.search(r'[HD][0-9]+[HD][0-9]+',pair[0]).group()) if aa != ''][0]
                zn_num2 = [int(aa) for aa in re.split(r'[HD]',re.search(r'[HD][0-9]+[HD][0-9]+',pair[1]).group()) if aa != ''][0]
                zn_stubs_1 = create_stubs(pdbfile+pair[0])
                zn_stubs_2 = create_stubs(pdbfile+pair[1])
                _ , zn1 = pair_zn_coorsdinate(pdbfile+pair[0])
                _ , zn2 = pair_zn_coorsdinate(pdbfile+pair[1])
                zn_coor1 = np.hstack((np.array(zn1)[None, :],np.ones([1])[None, :]))
                zn_coor2 = np.hstack((np.array(zn2)[None, :],np.ones([1])[None, :]))
                zn_coor_new_1 = wefold_stubs[scaffold_pose_seq_num*symm_num+zn_num1-1]@ (np.linalg.inv(zn_stubs_1[0])@zn_coor1.T)
                zn_coor_new_2 = wefold_stubs[scaffold_pose_seq_num*((symm_num+1)%symm)+zn_num2-1]@ (np.linalg.inv(zn_stubs_2[0])@zn_coor2.T)
                zn_coor_new = ((zn_coor_new_2+zn_coor_new_1)/2).reshape((-1,4))[0][:3]
                wefold_pose.append_pose_by_jump(zn_pose,1)
                wefold_pose.set_xyz(AtomID(wefold_pose.residue(len(wefold_pose.sequence())).atom_index("ZN"), len(wefold_pose.sequence())),xyz.assign(zn_coor_new[0],zn_coor_new[1],zn_coor_new[2]))

    wefold_pose.dump_pdb(new_pose_name)


        