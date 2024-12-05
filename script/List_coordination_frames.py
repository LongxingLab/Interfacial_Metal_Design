import pandas as pd
import numpy as np
import json
import argparse
import os
from scipy.spatial.distance import squareform
import re

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

def wefold_bb_stubs(coords):
    assert(coords.ndim == 3) # resi : atoms : xyz
    n = coords[:, 0, :3]
    ca = coords[:, 1, :3]
    c = coords[:, 2, :3]
    stubs = np.zeros((len(n),4,4),dtype=coords.dtype)
    stubs[:,3,3] = 1
    e1 = c[:, :3] - ca[:, :3]
    e1 /= np.linalg.norm(e1, axis=1)[:, None]
    e3 = np.cross(e1, n[:, :3] - ca[:, :3])
    e3 /= np.linalg.norm(e3, axis=1)[:, None]
    e2 = np.cross(e3, e1)
    e2 /= np.linalg.norm(e2, axis=1)[:, None]
    stubs[:, :3, 0] = e1
    stubs[:, :3, 1] = e2
    stubs[:, :3, 2] = e3
    stubs[:, :3, 3] = ca[:, :3]
    assert np.allclose(np.linalg.det(stubs), 1)
    return stubs

parser = argparse.ArgumentParser()
parser.add_argument('-scaffold_pdb', type=str, default=None, help='name of the scaffold pdb file')
parser.add_argument('-matchlist', type=str, default=None, help="name of the Match pdb list")
parser.add_argument('-output_path', type=str, default=None, help="the path of the output file")
parser.add_argument('-output_name', type=str, default=None, help="the name of the output file")
args = parser.parse_args()

pdblist_file = args.matchlist
if args.output_name == None:
    output_name = 'output'
else:
    output_name = args.output_name

with open(pdblist_file) as f:
    match_lines = [line.strip() for line in f]
_ = match_lines.pop(0) 

coords,_ = conformation_from_file(args.scaffold_pdb,atoms=['N','CA','C','O'],residue_major=True)
res1_wefold_stubs = wefold_bb_stubs(coords)[0]

ref_frame = np.linalg.inv(res1_wefold_stubs)

df = pd.read_csv(pdblist_file, sep='\s+')

atom1 = df[['atom1_x', 'atom1_y', 'atom1_z']].to_numpy()
atom2 = df[['atom2_x', 'atom2_y', 'atom2_z']].to_numpy()
zn    = df[['ZN_x', 'ZN_y', 'ZN_z']].to_numpy()

vec1 = atom1 - zn
vec1 /= np.linalg.norm(vec1, axis=-1)[...,None]
vec2 = atom2 - zn
vec2 /= np.linalg.norm(vec2, axis=-1)[...,None]
x_axis = vec1 + vec2
x_axis /= np.linalg.norm(x_axis, axis=-1)[...,None]

z_axis = np.cross(vec1, vec2)
z_axis /= np.linalg.norm(z_axis, axis=-1)[...,None]

y_axis = np.cross(z_axis, x_axis)

xforms = np.zeros((len(zn),4,4))
xforms[:,3,3]  = 1
xforms[:,:3,0] = x_axis
xforms[:,:3,1] = y_axis
xforms[:,:3,2] = z_axis
xforms[:,:3,3] = zn

def xform_magnitude_sq_fast( trans_err2, traces, lever2 ):

    # trans_part = rts[...,:3,3]
    # err_trans2 = np.sum(np.square(trans_part), axis=-1)

    # rot_part = rts[...,:3,:3]
    # traces = np.trace(rot_part,axis1=-1,axis2=-2)
    cos_theta = ( traces - 1 ) / 2

    # We clip to 0 here so that negative cos_theta gets lever as error
    clipped_cos = np.clip( cos_theta, 0, 1)

    err_rot2 = ( 1 - np.square(clipped_cos) ) * lever2

    # err = np.sqrt( err_trans2 + err_rot2 )
    err =  trans_err2 + err_rot2 

    return err


def mm2(inv_xform, xforms, traces, trans_err2):

    a = inv_xform
    b = xforms

    # leaving out the 4th term because we know it's 0
    traces[:] = np.sum( a[0,:3] * b[:,:3,0], axis=-1 )
    traces += np.sum( a[1,:3] * b[:,:3,1], axis=-1 )
    traces += np.sum( a[2,:3] * b[:,:3,2], axis=-1 )

    # we know the 4th term here has a 1 in b
    trans_err2[:] = np.square(np.sum( a[0,:3] * b[:,:3,3], axis=-1) + a[0,3])
    trans_err2 += np.square(np.sum( a[1,:3] * b[:,:3,3], axis=-1) + a[1,3])
    trans_err2 += np.square(np.sum( a[2,:3] * b[:,:3,3], axis=-1) + a[2,3])

def xform_dist_matrix(xforms, lever, inverse_xforms = None):
    
    num_xforms = len(xforms)
    lever2 = lever*lever
    
    if ( xforms.dtype != np.float32 ):
        xforms = xforms.astype(np.float32)
        
    xforms_sym = np.copy(xforms)
    
    xforms_sym[:,:3,1] = -1 * xforms_sym[:,:3,1]
    xforms_sym[:,:3,2] = -1 * xforms_sym[:,:3,2]

    if ( inverse_xforms is None ):
        inverse_xforms = np.linalg.inv(xforms)
    else:
        if ( inverse_xforms.dtype != inverse_xforms ):
            inverse_xforms = inverse_xforms.astype(np.float32)
            
    distance_list = []
    for idx in range(num_xforms-1):
        
        traces_1 = np.zeros(num_xforms-1-idx, dtype=np.float32)
        trans_err2_1 = np.zeros(num_xforms-1-idx, dtype=np.float32)

        traces_2 = np.zeros(num_xforms-1-idx, dtype=np.float32)
        trans_err2_2 = np.zeros(num_xforms-1-idx, dtype=np.float32)
        
        mm2(inverse_xforms[idx], xforms[idx+1:], traces_1, trans_err2_1)
        distances_1 = xform_magnitude_sq_fast( trans_err2_1, traces_1, lever2 )
        
        mm2(inverse_xforms[idx], xforms_sym[idx+1:], traces_2, trans_err2_2)
        distances_2 = xform_magnitude_sq_fast( trans_err2_2, traces_2, lever2 )

        distances = np.min(np.c_[distances_1,distances_2], axis=-1)
        
        distance_list.append(distances)
        
    return squareform(np.hstack(distance_list))

def xform_to_str(xform):
    return '_'.join([str(ii) for ii in xform[:3,:].flatten(order='F')])

def clustering_report(xforms,
                      dist_matrix,
                      centers,
                      members,
                      labels):

    output = {}
    for idx, curr_center in enumerate(centers):
        cluster_members = members[idx]

        #local_matrix = dist_matrix[cluster_members][:,cluster_members]
        dist_center = dist_matrix[curr_center,cluster_members]
        
        mean = np.mean(dist_center)
        largest_dist = np.max(dist_center)

        center_name = labels[curr_center]
        info = {}
        info["avg_dist"] = float(mean)
        info["max_dist"] = float(largest_dist)
        info['center_xform'] = xform_to_str(xforms[curr_center])
        info['center_match'] = center_name.split()[0]
        cluster_labels = []
        for mem_idx in cluster_members:
            cluster_labels.append(labels[mem_idx].split()[0])
        info["descriptions"] = cluster_labels
        info['residue_index'] = list(map(int,re.split(r'[HD]',center_name.split()[0].split('_')[3])[1:]))
        info['cluster_token'] = 1
        info['matching_clusters'] = [1]
        output["Zn_Cluster{:04d}".format(idx)] = info

    return output

def compute_cluster_centers(dist_matrix,
                            members):
    print("Reassigning centers to true centers of the cluster")
    centers = []
    min_distances = np.zeros(dist_matrix.shape[0])
    for cluster_members in members:
        local_matrix = dist_matrix[cluster_members][:,cluster_members]
        local_idx = np.argmin(local_matrix.sum(axis=0))
        global_idx = cluster_members[local_idx]
        centers.append( global_idx )
        min_distances[cluster_members] = local_matrix[local_idx]

    return centers, min_distances


def fast_clustering(dist_matrix,
                  convergence_dist=0.1,
                  recompute_center=True,
                  min_nclusters=-1,
                  max_nclusters=-1):
    length = dist_matrix.shape[0]
    min_distances = np.zeros(length)
    min_distances.fill(np.finfo(float).max)
    assignments = np.zeros(length, np.int32)
    centers = []
    members = []
    curr_center = 0
    while np.max(min_distances) > convergence_dist or (min_nclusters != -1 and len(center_idx) < min_nclusters):
        need_update = dist_matrix[curr_center] < min_distances
        assignments[need_update] = curr_center
        min_distances[need_update] = dist_matrix[curr_center][need_update]
        centers.append(curr_center)
        curr_center = np.argmax(min_distances)
 
        if max_nclusters != -1 and len(centers) > max_nclusters:
            break
    for idx, curr_center in enumerate(centers):
        cluster_members = np.nonzero(assignments == curr_center)[0]
        members.append(cluster_members)
    print("Total culsters identified using cutoff value {} is {}".format(convergence_dist, len(centers)))
    if recompute_center:
        centers, min_distances = compute_cluster_centers(dist_matrix, members)
 
    return centers, members, min_distances


dist_matrix = xform_dist_matrix(xforms, 4.0)
centers, members, min_distances = fast_clustering(dist_matrix, 2.0)
reports = clustering_report(ref_frame @ xforms, dist_matrix, centers, members, match_lines)
with open(args.output_path+'/'+output_name+'.dat', 'w') as f:
    f.write(json.dumps(reports, indent=4))









