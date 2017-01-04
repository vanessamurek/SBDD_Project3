import sys
import random
import itertools


def parse_sdf_to_ligand_list(sdf_file_path):
    """ Parses SDF file and returns a list of all the molecules contained within. """
    sdf = open(sdf_file_path, 'r')
    all_molecules = []
    current_molecule = ''
    for line in sdf:
        current_molecule += line
        if '$$$$' in line:
            all_molecules.append(current_molecule)
            current_molecule = ''
    return all_molecules

def extract_random_ligands(complete_ligand_list, number_of_ligands):
    """ Takes ligand list and extracts desired number of random ligands, constructing a new ligand list. """
    new_ligand_list = []
    for _ in itertools.repeat(None, number_of_ligands):
        current_ligand = random.choice(complete_ligand_list)
        complete_ligand_list.remove(current_ligand)
        new_ligand_list.append(current_ligand)
    return new_ligand_list

def write_sdf(ligand_list):
    """ Writes ligand list to SDF file format. """
    new_sdf = open('random_ligands.sdf', 'w')
    for ligand in ligand_list:
        new_sdf.write("%s" % ligand)

if __name__ == "__main__":
    # Specify file path to SDF file containing decoys as first argument
    sdf_structure_file_path = sys.argv[1]

    # Parse ligands to list from SDF
    all_ligands = parse_sdf_to_ligand_list(sdf_structure_file_path)

    # Specify desired number of decoys to extract
    number_decoys_to_extract = 250

    # Filter desired number of random ligands from list
    extracted_ligands = extract_random_ligands(all_ligands, number_decoys_to_extract)

    # Write extracted ligands to file
    write_sdf(extracted_ligands)




