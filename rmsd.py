from __future__ import print_function
from sdf_reader import SDFReader
import numpy.linalg
import sys
__author__ = 'Vanessa Murek'

def error(*params):
    print(*params, file=sys.stderr)

def rmsd(molecule_a, molecule_b):
    atoms_in_molecule = len(molecule_a)
    build_term = 0
    for atom in range(0, (atoms_in_molecule)):
        atom_pos_a = numpy.array((molecule_a[atom].position.x, molecule_a[atom].position.y, molecule_a[atom].position.z))
        atom_pos_b = numpy.array((molecule_b[atom].position.x, molecule_b[atom].position.y, molecule_b[atom].position.z))
        eucl_dist = numpy.linalg.norm(atom_pos_a-atom_pos_b)
        build_term += eucl_dist**2
    rmsd = numpy.sqrt((1/atoms_in_molecule)*build_term)
    return rmsd

if __name__ == "__main__":
    # check if we have everything we need
    if len(sys.argv) < 2:
        error("Not enough arguments provided.")
        error("Usage:")
        error("  python rmsd.py <input-sdf>")
        sys.exit(-1)
    # the first given argument is the sdf file that we will parse and use
    file_path_reference = sys.argv[1]
    file_path_autodock_a = sys.argv[2]
    file_path_caddsuite_a = sys.argv[3]
    file_path_autodock_b = sys.argv[4]
    file_path_caddsuite_b = sys.argv[5]

    molecules_reference = SDFReader.load(file_path_reference)
    molecules_autodock = SDFReader.load(file_path_autodock_a)
    molecules_caddsuite = SDFReader.load(file_path_caddsuite_a)
    molecules_autodock_b = SDFReader.load(file_path_autodock_b)
    molecules_caddsuite_b = SDFReader.load(file_path_caddsuite_b)


    print(str(molecules_reference[0].name))
    print("Reference against autodock chain A:")
    print("RMSD between molecule " + str(molecules_reference[0].name).split(",")[0] + " and molecule " +
                  str(molecules_autodock[0].name).split(",")[0] + " is " +
                  str(rmsd(molecules_reference[0], molecules_autodock[0])) + " Angstrom")

    print("RMSD between molecule " + str(molecules_reference[1].name).split(",")[0] + " and molecule " +
                  str(molecules_autodock[1].name).split(",")[0] + " is " +
                  str(rmsd(molecules_reference[1], molecules_autodock[1])) + " Angstrom")

    print("Reference against caddsuite chain A:")
    print("RMSD between molecule " + str(molecules_reference[0].name).split(",")[0] + " and molecule " +
                  str(molecules_caddsuite[0].name).split(",")[0] + " is " +
                  str(rmsd(molecules_reference[0], molecules_caddsuite[0])) + " Angstrom")

    print("RMSD between molecule " + str(molecules_reference[1].name).split(",")[0] + " and molecule " +
                  str(molecules_caddsuite[1].name).split(",")[0] + " is " +
                  str(rmsd(molecules_reference[1], molecules_caddsuite[1])) + " Angstrom")

    print("Reference against autodock chain B:")
    print("RMSD between molecule " + str(molecules_reference[0].name).split(",")[0] + " and molecule " +
                  str(molecules_autodock_b[0].name).split(",")[0] + " is " +
                  str(rmsd(molecules_reference[0], molecules_autodock_b[0])) + " Angstrom")

    print("RMSD between molecule " + str(molecules_reference[1].name).split(",")[0] + " and molecule " +
                  str(molecules_autodock_b[1].name).split(",")[0] + " is " +
                  str(rmsd(molecules_reference[1], molecules_autodock_b[1])) + " Angstrom")

    print("Reference against caddsuite chain B:")
    print("RMSD between molecule " + str(molecules_reference[0].name).split(",")[0] + " and molecule " +
                  str(molecules_caddsuite_b[0].name).split(",")[0] + " is " +
                  str(rmsd(molecules_reference[0], molecules_caddsuite_b[0])) + " Angstrom")

    print("RMSD between molecule " + str(molecules_reference[1].name).split(",")[0] + " and molecule " +
                  str(molecules_caddsuite_b[1].name).split(",")[0] + " is " +
                  str(rmsd(molecules_reference[1], molecules_caddsuite_b[1])) + " Angstrom")

    #print("Loaded %d molecules from %s." % (len(molecules), file_name))

