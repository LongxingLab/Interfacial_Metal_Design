import numpy as np
import copy
import argparse
from collections import defaultdict
import sys
import re

description = 'Get the sequence of specified chain for pdb files ...'
parser = argparse.ArgumentParser( description = description )
parser.add_argument('-pdb', type=str, help='name of the input pdb file')
parser.add_argument('-symmetry', type=str, required=True, help='the symmetry definition of the input pdb file')
parser.add_argument('-num_zincs', type=int, required=True, help='the number of zincs associated with each chain')
parser.add_argument('-not_design_threshold', type=float, default=10.0, help='above not design')
parser.add_argument('-full_pose', default=False, action='store_true', help='for full pose, C2')
parser.add_argument('-must_design_threshold', type=float, default=6.0, help='above not design')
parser.add_argument('-pos_file', type=str, help='pos file of scaffold')
parser.add_argument('-sep', default=',', type=str, help='a list file of all pdb files')
parser.add_argument('-output', type=str, help="the name of the output file")
args = parser.parse_args()

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

def conformation_from_pdblines(pdblines, atoms=['N', 'CA', 'C', 'O', 'CB'], residue_major=True, dtype=np.float64):

    sequence = ''
    for line in pdblines:
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

def psudo_CB(N, CA=None, C=None):

    # assume the first 
    if CA is None:
        N, CA, C = N[...,0,:], N[...,1,:], N[...,2,:]
    #
    # calc the position of the CB atom
    b = CA - N
    c = C - CA
    a = np.cross(b, c)
    CB = -0.58273431*a + 0.56802827*b - 0.54067466*c + CA

    return CB        
        
with open(args.pos_file) as f:
    lines = [line.strip() for line in f]

designable_res = set([int(ii) for ii in lines[1].split(':')[-1].split()])

pdb = args.pdb

results = []

m = re.match(r'C([1-9])$', args.symmetry)
if m == None:
    print("Symmetry definition error!")
    exit(0)
else:
    sym_num = int(m.group(1))

#sym_num = 10



base_name = pdb.split('/')[-1]

if base_name.endswith('.pdb.gz'):
    tag = base_name[:-7]
elif base_name.endswith('.pdb'):
    tag = base_name[:-4]
else:
    exit(0)

scaffold_tag = tag.split('___')[0]

with open(pdb) as f:
    pdblines = [line for line in f if line.startswith('ATOM')]
coords, _ = conformation_from_pdblines(pdblines,atoms = ['N', 'CA', 'C', 'O'],residue_major = True)

pseudo_CB = psudo_CB(coords)
monomer_len = len(pseudo_CB) // sym_num 
chainA = pseudo_CB[:monomer_len]
chain_other = pseudo_CB[monomer_len:]

dist_squared = np.sum(np.square(chainA[:,None,:] - chain_other[None,...]), axis=-1)

not_design_res = np.where( np.all(dist_squared > args.not_design_threshold * args.not_design_threshold, axis=-1) )[0] + 1
must_design_res = np.where( np.any(dist_squared < args.must_design_threshold * args.must_design_threshold, axis=-1) )[0] + 1

not_design_res = set(not_design_res)
must_design_res = set(must_design_res)
designable_res = designable_res

designable_res = designable_res - not_design_res
designable_res = designable_res | must_design_res

designable_res = list(designable_res)

if args.full_pose:
    temp_designable_res = np.array(designable_res)
    for ii in range(sym_num-1):
        temp_designable_res += monomer_len + args.num_zincs
        designable_res += list(temp_designable_res)

designable_res.sort()
designable_res = args.sep.join([str(ii) for ii in designable_res])

results.append(tag + ' ' + designable_res)

if args.output is None:
    for ii in results:
        print(ii)
else:
    with open(args.output, 'w') as f:
        f.write('\n'.join(results))
        f.write('\n')
