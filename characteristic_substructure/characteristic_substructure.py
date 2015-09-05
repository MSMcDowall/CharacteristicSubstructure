from .algorithm import CSAlgorithm
import argparse

# Main module for the application
# main method acts as entry point for application


def main():
    """
    Creates the algorithm object with the correct parameters and calls methods based on the command line arguments.

    Uses the argparse module to parse the command line arguments that are given
    The file name for the SMILES file is mandatory, the flags and threshold specification are optional
    :return: None
    """
    parser = argparse.ArgumentParser(description="compare the structure of molecules in a SMILES file")
    parser.add_argument("smiles_file", help="a file containing SMILES strings")
    parser.add_argument("threshold", nargs='?', default=0.8, type=float,
                        help="the relative frequency of structures in the molecules")
    parser.add_argument("-c", "--characteristic", action="store_true",
                        help="display the characteristic substructure only")
    parser.add_argument("-r", "--representative", action="store_true",
                        help="display all of the representative substructures")
    args = parser.parse_args()
    cs = CSAlgorithm(smiles_file=args.smiles_file)
    if args.threshold:
        cs.threshold = args.threshold

    c_structure = None
    all_structures = None
    if args.characteristic:
        c_structure = cs.find_characteristic_substructure()
    if args.representative:
        all_structures = cs.find_all_representative_structures()
    if not args.characteristic and not args.representative:
        c_structure = cs.find_characteristic_substructure()
        all_structures = cs.find_all_representative_structures()

    # Display the structures in the command line
    if c_structure:
        print 'Characteristic Substructure'
        print c_structure.adjacency_dictionary.keys()
    if all_structures:
        print 'Representative Structures'
        print all_structures

