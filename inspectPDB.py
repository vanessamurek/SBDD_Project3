import sys
from Bio.PDB import *
from pprint import pprint

# TODO: inspect structure -> How many chains? Duplicates? Solvents? Water?
# ConstructionWarning: 4 mal -> Chain A / B is discontinuous at line 6679, 6743, 6792, 7272 (markiert Übergänge zwischen den chains)
# 1 model (	the model serial number when a single coordinate entry contains multiple structures), 2 chains (A und B)
# Residues rausgefiltert und gezählt

# TODO: Identify reference ligand -> Must be included in docking
# TODO: Remove components to be excluded (solvents, water)
# TODO: In case of alternative records, choose the one with higher occupancy
# TODO: Multiple chains? Treat them seperately, create one PDB structure for each chain
from pandas.io.pytables import Selection

if __name__ == "__main__":
    pdb_structure_file_name = sys.argv[1]
    structure_object = PDBParser(PERMISSIVE=True)
    structure = structure_object.get_structure("3B7E", pdb_structure_file_name)

    model = structure[0]

    # TODO: Implement check for occupancy values: Only 1.00 and 0.50?
    #testcom
    
    class NMROutputSelector2(Select):
        # Inherit methods from Select class
        def accept_atom(self, atom):
            if (not atom.is_disordered()) or atom.get_altloc() == 'A':
                atom.set_altloc(' ')  # Eliminate alt location ID before output.
                return True
            else:  # Alt location was not one to be output.
                return False


    residue_to_remove = []
    chain_to_remove = []
    for chain in model:
        for residue in chain:
            if residue.id[0] != ' ' and residue.id[0] != 'H_ZMR':
                residue_to_remove.append((chain.id, residue.id))
        if len(chain) == 0:
            chain_to_remove.append(chain.id)

    for residue in residue_to_remove:
        model[residue[0]].detach_child(residue[1])

    for chain in chain_to_remove:
        model.detach_child(chain)

    io = PDBIO()
    io.set_structure(structure[0]['A'])
    io.save("3b7e_chain_a.pdb", select=NMROutputSelector2())

    io.set_structure(structure[0]['B'])
    io.save("3b7e_chain_b.pdb", select=NMROutputSelector2())
