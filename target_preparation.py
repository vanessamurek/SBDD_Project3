import sys
from Bio.PDB import *


class LocationSelector(Select):
    def accept_atom(self, atom):
        """ For atoms with alternate locations, selects first alternative and removes other locations. """
        if (not atom.is_disordered()) or atom.get_altloc() == 'A':
            atom.set_altloc(' ')
            return True
        else:
            return False

def get_model_from_pdb(pdb_file_path, pdb_id):
    """ Creates a structure object using PDB file path and ID and returns the first model. """
    structure_object = PDBParser(PERMISSIVE=True)
    structure = structure_object.get_structure(pdb_id, pdb_file_path)
    return structure[0]


def remove_het_atoms(model, reference_ligand):
    """ Removes all HET atoms from the model, except for the specified ligand. """
    residue_to_remove = []
    chain_to_remove = []
    for chain in model:
        for residue in chain:
            if residue.id[0] != ' ' and residue.id[0] != reference_ligand:
                residue_to_remove.append((chain.id, residue.id))
        if len(chain) == 0:
            chain_to_remove.append(chain.id)
    for residue in residue_to_remove:
        model[residue[0]].detach_child(residue[1])
    for chain in chain_to_remove:
        model.detach_child(chain)
    return model


if __name__ == "__main__":
    # Specify file path to PDB structure as first argument
    pdb_structure_file_path = sys.argv[1]

    # Define PDB ID and get first model from PDB structure
    pdb_id = '3B7E'
    model = get_model_from_pdb(pdb_structure_file_path, pdb_id)

    # Define residue ID for reference ligand
    reference_ligand = 'H_ZMR'

    # Remove solvents, water molecules and other HET atoms, except for reference ligand from model
    model_without_het_atoms = remove_het_atoms(model, reference_ligand)

    # Create a PDB structure of each of the two chains, handling alternative locations
    io = PDBIO()
    io.set_structure(model_without_het_atoms['A'])
    io.save("3b7e_chain_a_clean.pdb", select=LocationSelector())

    io.set_structure(model_without_het_atoms['B'])
    io.save("3b7e_chain_b_clean.pdb", select=LocationSelector())
