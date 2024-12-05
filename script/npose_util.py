#!/usr/bin/env python

import os
import sys
import math
import pandas as pd
import numpy as np
import warnings
import gzip
import struct
import itertools

from collections import OrderedDict

if ( hasattr(os, 'ATOM_NAMES') ):
    assert( hasattr(os, 'PDB_ORDER') )

    ATOM_NAMES = os.ATOM_NAMES
    PDB_ORDER = os.PDB_ORDER
else:
    ATOM_NAMES=['N', 'CA', 'CB', 'C', 'O']
    PDB_ORDER = ['N', 'CA', 'C', 'O', 'CB']

_byte_atom_names = []
_atom_names = []
for i, atom_name in enumerate(ATOM_NAMES):
    long_name = " " + atom_name + "       "
    _atom_names.append(long_name[:4])
    _byte_atom_names.append(atom_name.encode())

    globals()[atom_name] = i

R = len(ATOM_NAMES)

def gzopen(name, mode="rt"):
    if (name.endswith(".gz")):
        return gzip.open(name, mode)
    else:
        return open(name, mode)


_space = " ".encode()[0]

def space_strip(string):
    start = 0
    while(start < len(string) and string[start] == _space):
        start += 1
    end = len(string)
    while( end > 0 and string[end-1] == _space):
        end -= 1
    return string[start:end]


def byte_startswith( haystack, needle ):
    if ( len(haystack) < len(needle) ):
        return False
    for i in range(len(needle)):
        if ( haystack[i] != needle[i] ):
            return False
    return True


def byte_equals( haystack, needle ):
    if ( len(haystack) != len(needle) ):
        return False
    for i in range(len(needle)):
        if ( haystack[i] != needle[i] ):
            return False
    return True


def getline( bytess, start ):
    cur = start
    while (cur < len(bytess) and bytess[cur] != 10 ):
        cur += 1
    cur += 1
    return bytess[start:cur], cur

# ord(" ") == 32
# ord("-") == 45
# ord(".") == 46
# ord("0") == 48
# ord("9") == 57


def stof(string):
    multiplier = 0

    parsed = False
    negative = 1
    start = 0
    end = len(string) - 1
    for i in range(len(string)):
        char = string[i]
        # print(char)
        if ( not parsed ):
            if ( char == 32 ): # " "
                start = i + 1
                continue
            if ( char == 45 ): # "-"
                start = i + 1
                negative = -1
                continue
        if ( char == 32 ): # " "
            break
        if ( char == 46 ): # "."
            multiplier = np.float64(1)
            parsed = True
            continue
        if ( char >= 48 and char <= 57 ): # 0 9
            parsed = True
            multiplier /= np.float64(10)
            end = i
            continue
        print("Float parse error! Unrecognized character: ", char)
        assert(False)

    if ( not parsed ):
        print("Float parse error!")
        assert(False)

    result = np.float64(0)

    if ( multiplier == 0 ):
        multiplier = 1
    for i in range(end, start-1, -1):
        char = string[i]
        if ( char == 46 ): # "."
            continue
        value = np.float64(char - 48) # 0

        result += value * multiplier
        multiplier *= np.float64(10)

    result *= negative

    return np.float32(result)


def parse_aa(name3):

    a = name3[0]
    b = name3[1]
    c = name3[2]

    if ( a <= 73 ): # I
        if ( a <= 67 ): # C
            if ( a == 65 ): # A
                if ( b == 83 ):
                    if ( c == 78 ):
                        return 78
                    if ( c == 80 ):
                        return 68
                    return 88
                if ( b == 76 and c == 65 ):
                    return 65
                if ( b == 82 and c == 71 ):
                    return 82
            if ( a == 67 and b == 89 and c == 83 ): # C
                    return 67
        else:
            if ( a == 71 ): # G
                if ( b != 76 ):
                    return 88
                if ( c == 85 ):
                    return 69
                if ( c == 89 ):
                    return 71
                if ( c == 78 ):
                    return 81
            if ( a == 72 and b == 73 and c == 83 ): # H
                    return 72
            if ( a == 73 and b == 76 and c == 69 ): # I
                    return 73
    else: 
        if ( a <= 80 ): # P
            if ( a == 76 ): # L
                if ( b == 69 and c == 85 ):
                    return 76
                if ( b == 89 and c == 83 ):
                    return 75
            if ( a == 77 ): # M
                if ( b == 69 and c == 84 ):
                    return 77
            if ( a == 80 ): # P
                if ( b == 72 and c == 69 ):
                    return 70
                if ( b == 82 and c == 79 ):
                    return 80
        else:
            if ( a == 83 and b == 69 and c == 82 ): # S
                    return 83
            if ( a == 84 ): # T
                if ( c == 82 ):
                    if ( b == 72 ):
                        return 84
                    if ( b == 89 ):
                        return 89
                    return 88
                if ( b == 82 and c == 80 ):
                    return 87
            if ( a == 86 and b == 65 and c == 76 ): # V
                    return 86
    return 88


# _atom = "ATOM".encode()
_null_line = "ATOMCB00000000000000000000000000000000000000000000000000000000"#.encode()
_null_line_size = len(_null_line)

# _CB = "CB".encode()
# _empty_bytes = "".encode()

# Switches to next residue whenever
#  Line doesn't start with atom
#  Line isn't long enough
#  Res/resnum/chain changes

def read_npose_from_data( data, null_line_atom_names, NCACCBR, scratch_residues, scratch_chains, scratch_aa):


    _null_line = null_line_atom_names[:_null_line_size]
    array_byte_atom_names = null_line_atom_names[_null_line_size:]

    _CB = _null_line[4:6]
    _atom = _null_line[:4]
    _empty_bytes = _null_line[:0]

    N = NCACCBR[0]
    CA = NCACCBR[1]
    C = NCACCBR[2]
    CB = NCACCBR[3]
    R = NCACCBR[4]

    byte_atom_names = []
    for i in range(len(array_byte_atom_names)//4):
        aname = array_byte_atom_names[i*4:i*4+4]
        byte_atom_names.append(space_strip(aname))

    seqpos = 0
    res_ident = _empty_bytes
    res_has_n_atoms = 0
    next_res = False
    scratch_residues[0].fill(0)

    # lines.append(_nulline)

    cursor = 0
    keep_going = True
    while keep_going:
        line, cursor = getline(data, cursor)
        if ( cursor >= len(data) ):
            keep_going = False
            line = _null_line

        # print(iline)

        if ( not byte_startswith(line, _atom) or len(line) < 54 ):
            next_res = True
            res_ident = _empty_bytes
            continue

        ident = line[17:27]
        if ( not byte_equals( ident, res_ident ) ):
            next_res = True

        if ( next_res ):
            if ( res_has_n_atoms > 0 ):

                res = scratch_residues[seqpos]
                if ( res_has_n_atoms != R ):
                    missing = np.where(res[:,3] == 0)[0]

                    # We only know how to fix missing CB
                    first_missing = byte_atom_names[missing[0]]
                    if ( len(missing) > 1 or not byte_equals(first_missing,  _CB) ):
                        print("Error! missing atoms:")
                        for i in range(len(missing)):
                            print(byte_atom_names[missing[i]])
                        print("in residue:")
                        print(res_ident)
                        assert(False)

                    # Fixing CB
                    xform = get_stub_from_n_ca_c(res[N,:3], res[CA,:3], res[C,:3])
                    res[CB] = get_CB_from_xform( xform )


                seqpos += 1
                #If we run out of scratch, double its size
                if ( seqpos == len(scratch_residues) ):
                    old_size = len(scratch_residues)
                    new_scratch = np.zeros((old_size*2, R, 4), np.float32)
                    for i in range(old_size):
                        new_scratch[i] = scratch_residues[i]
                    scratch_residues = new_scratch

                    new_scratch2 = np.zeros((old_size*2), np.byte)
                    for i in range(old_size):
                        new_scratch2[i] = scratch_chains[i]
                    scratch_chains = new_scratch2

                    new_scratch2 = np.zeros((old_size*2), np.byte)
                    for i in range(old_size):
                        new_scratch2[i] = scratch_aa[i]
                    scratch_aa = new_scratch2


                scratch_residues[seqpos].fill(0)

            res_ident = ident
            res_has_n_atoms = 0
            next_res = False


        # avoid parsing stuff we know we don't need
        if ( res_has_n_atoms == R ):
            continue

        atom_name = space_strip(line[12:16])

        # figure out which atom we have
        atomi = -1
        for i in range( R ):
            if ( byte_equals( atom_name, byte_atom_names[i] ) ):
                atomi = i
                break
        if ( atomi == -1 ):
            continue

        res = scratch_residues[seqpos]
        if ( res[atomi,3] != 0 ):
            print("Error! duplicate atom:")
            print( atom_name )
            print("in residue:" )
            print( res_ident )
            assert(False)

        res_has_n_atoms += 1

        res[atomi,0] = stof(line[30:38])
        res[atomi,1] = stof(line[38:46])
        res[atomi,2] = stof(line[46:54])
        res[atomi,3] = 1

        if ( res_has_n_atoms == 1 ):
            scratch_chains[seqpos] = line[21]
            scratch_aa[seqpos] = parse_aa(line[17:20])

    to_ret = np.zeros((seqpos, R, 4))
    for i in range(seqpos):
        to_ret[i] = scratch_residues[i]
    
    return to_ret.reshape(-1, 4), scratch_residues, scratch_chains, scratch_aa


g_scratch_residues = np.zeros((1000,R,4), np.float32)
g_scratch_chains = np.zeros((1000), np.byte)
g_scratch_aa = np.zeros((1000), np.byte)
# for i in range(1000):
#     g_scratch_residues.append(np.zeros((R,4), np.float32))

_array_byte_atom_names = list(" "*len(_byte_atom_names)*4)
for i in range(len(_byte_atom_names)):
    name = _atom_names[i]
    for j in range(len(name)):
        _array_byte_atom_names[i*4+j] = name[j]
_array_atom_names = "".join(_array_byte_atom_names)

_null_line_atom_names = (_null_line + _array_atom_names).encode()

NCACCBR = np.zeros(5, np.int64)
if ( "N" in locals() ):
    NCACCBR[0] = N
if ( "CA" in locals() ):
    NCACCBR[1] = CA
if ( "C" in locals() ):
    NCACCBR[2] = C
if ( "CB" in locals() ):
    NCACCBR[3] = CB
NCACCBR[4] = R

def npose_from_file(fname, chains=False, aa=False):
    return npose_from_file_fast(fname, chains, aa)

def npose_from_file_fast(fname, chains=False, aa=False):
    with gzopen(fname, "rb") as f:
        data = f.read()
    return npose_from_bytes(data, chains, aa)

def npose_from_bytes(data, chains=False, aa=False):

    global g_scratch_residues
    global g_scratch_chains
    global g_scratch_aa

    npose, scratch, scratch_chains, scratch_aa = read_npose_from_data( data, _null_line_atom_names, NCACCBR, g_scratch_residues, 
                                                                                                g_scratch_chains, g_scratch_aa)
    np.around(npose, 3, npose)  # get rid of random numba noise

    g_scratch_residues = scratch
    g_scratch_chains = scratch_chains
    g_scratch_aa = scratch_aa

    output = [npose]
    if ( chains ):
        output.append(bytes(scratch_chains[:nsize(npose)]).decode("ascii"))
    if ( aa ):
        output.append(bytes(scratch_aa[:nsize(npose)]).decode("ascii"))
    if ( len(output) == 1 ):
        return output[0]
    return output

def cross(vec1, vec2):
    result = np.zeros(3)
    a1, a2, a3 = vec1[0], vec1[1], vec1[2]
    b1, b2, b3 = vec2[0], vec2[1], vec2[2]
    result[0] = a2 * b3 - a3 * b2
    result[1] = a3 * b1 - a1 * b3
    result[2] = a1 * b2 - a2 * b1
    return result

def get_stub_from_n_ca_c(n, ca, c):
    e1 = ca - n
    e1 /= np.linalg.norm(e1)

    e3 = cross( e1, c - n )
    e3 /= np.linalg.norm(e3)

    e2 = cross( e3, e1 )

    stub = np.identity(4)
    stub[:3,0] = e1
    stub[:3,1] = e2
    stub[:3,2] = e3
    stub[:3,3] = ca

    return stub

def extract_atoms(npose, atoms):
    return npose.reshape(-1, R, 4)[...,atoms,:].reshape(-1,4)

def get_CB_from_xform(xform):
    CB_pos = np.array([0.53474492, -0.76147505, -1.21079691, 1.0])
    return xform @ CB_pos
