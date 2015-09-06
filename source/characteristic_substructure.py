from algorithm import CSAlgorithm
import argparse
from datetime import datetime


def main():
    """
    Creates the algorithm object with the correct parameters and calls methods based on the command line arguments.

    Uses the argparse module to parse the command line arguments that are given
    The file name for the SMILES file is mandatory, the flags and threshold specification are optional
    :return: None
    """
    a = datetime.now()
    parser = argparse.ArgumentParser(description="compare the structure of molecules in a SMILES file")
    parser.add_argument("smiles_file", help="a file containing SMILES strings")
    parser.add_argument("threshold", nargs='?', default=0.8, type=float,
                        help="the relative frequency of structures in the molecules")
    parser.add_argument("-c", "--characteristic", action="store_true",
                        help="display the characteristic substructure only")
    parser.add_argument("-r", "--representative", action="store_true",
                        help="display all of the representative substructures")
    args = parser.parse_args()

    cs = CSAlgorithm()
    if args.threshold:
        cs.threshold = args.threshold

    smiles_set = data_input(args.smiles_file)
    paths = cs.find_graphs_paths(smiles_set)

    c_structure = None
    all_structures = None
    if args.characteristic:
        c_structure = cs.find_characteristic_substructure(paths)
    if args.representative:
        all_structures = cs.find_all_representative_structures(paths)
    if not args.characteristic and not args.representative:
        c_structure = cs.find_characteristic_substructure(paths)
        all_structures = cs.find_all_representative_structures(paths)
    b = datetime.now()
    # Display the structures in the command line
    if c_structure:
        print 'Characteristic Substructure'
        print c_structure.adjacency_dictionary.keys()
        data_output(cs.create_cs_results(), args.smiles_file[:-4] + str(args.threshold) +
                    '_CharacteristicSubstructure.txt')
    if all_structures:
        print 'Representative Structures'
        print all_structures
        data_output(cs.structures_output(all_structures),
                    args.smiles_file[:-4] + str(args.threshold) + '_RepresentativeStructures.txt')
    c = b-a
    print c.total_seconds()


def data_input(file_name):
    """
    Takes in a txt file, reads each line and stores the result in a list

    :return: a list of each of the lines of the text file
    """
    smiles_set = []
    reader = open(file_name, mode='rb')
    for line in reader:
        words = line.split(" ")
        smiles_set.append(words[0])
    reader.close()
    return smiles_set


def data_output(string_list, file_name):
    """
    Writes a list of strings to a file with the given filename

    :param string_list: list of display strings
    :return: None
    """
    display_string = ''.join(string_list)
    writer = open(file_name, mode='wb')
    writer.write(display_string)
    writer.close()


if __name__ == '__main__':
    main()