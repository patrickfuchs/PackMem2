#-*- coding: utf-8 -*-
"""Functions to create, analyse list variables."""
# R. Gautier 2015
# M. Zygadlo 2025

import numpy as np

from core import matrix as m

def get_glyc_lipids(lipid_list, RESNAME_GLYC):
    """
    Get the glycerol atom name(s) for each lipid in the bilayer.

    --------------------
    INPUT
    lipid_list:  list
        Contains the name of every lipid(s) in the bilayer given by the user
    RESNAME_GLYC: dictionary
        Contains as keys to lipid name and as keys the glycerol atom name
    
    --------------------
    OUTPUT
    str
        The glycerol atom name for each lipid given   
    """
    atom_mb =  []
    # Find the glycerol atom for all lipids selected
    for lip in lipid_list:
        atom_mb.append(RESNAME_GLYC[lip])
    # Delete the multi occurrences
    atom_mb = set(atom_mb)
    # Cast to string
    atom_mb = ' '.join(atom_mb)
    return atom_mb

def min_max_mean(data):
    """
    Find the minimum, maximum and mean of an array.

    --------------------
    INPUT
    data: numpy array

    --------------------
    OUTPUT
    float
        The minimum of the data
    float
        The maximum of the data
    float
        The mean of the data
    """
    mini = min(data)
    maxi = max(data)
    mean = np.mean(data)
    return mini, maxi, mean

def create_array(v1, v2, step):
    """
    Create numpy array from v1 to v2 by step step.

    --------------------
    INPUT
    v1: float
    v2: float
    step: float
    --------------------
    OUTPUT
    numpy array
        Containing floats from v1 to v2 by step step
    """
    return np.arange(v1, v2, step)

def create_arrayZ(residues, array_resids, resname_glyc, dist_suppl_Z, z_extr, up=True):
    """
    Create a dictionary of the Z position for upper and lower leaflet.

    --------------------
    INPUT
    residues: MDAnalysis ResidueGroup
        Contains the names of all the residues selected
    array_resids: numpy array
        Contains all the residue numbers selected
    resname_glyc: dictionary
        Contains the residue name as key and the reference atom as value
    dist_suppl_Z: float
        Supplementary distance from the z coord
    z_extr: float
        Either maximum or minimum value of z coord
    up: boolean
        If we are on the upper leaflet
    --------------------
    OUTPUT
    dictionary
        For each resid, it contains floats ranging from z_coord to zmin
        or from zmax to z_coord by 1.0 steps
    """
    leaflet_arrayZ = {}
    for resid in array_resids:
        # Select the residue with the selected res ID
        selec_residue = residues[residues.resids == resid][0]
        # Get the central atom of the residue
        glyc = get_glyc_lipids([selec_residue.resname], resname_glyc)
        # Get the atom
        atom = selec_residue.atoms[selec_residue.atoms.names == glyc][0]
        # Get the Z coord
        z_coord = atom.position[2]
        if up:
            tmp = create_array(round(z_coord - dist_suppl_Z, 2),
                                        round(z_extr +1.0, 2), m.SIZE)
        else:
            tmp = create_array(round(z_coord + dist_suppl_Z, 2),
                                            round(z_extr - 1.0, 2), -m.SIZE)
        # Reverse it
        tmp = np.flip(tmp)
        leaflet_arrayZ[resid] = tmp
    return leaflet_arrayZ
