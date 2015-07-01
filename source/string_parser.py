import re

from source import molecule, tokenizer

smiles_string = ''
molecule_graph = molecule.Molecule(smiles_string)
_tokens = tokenizer.tokenize_smiles(smiles_string)
_previous_atom = None
_break_points = {}
#_branch_root = #stack

# Regular expressions used in testing
letters = re.compile('([A-Za-z][a-z]?)')
numbers = re.compile('([0-9])')

# Read in string from file
def create_atom(token, isotope=False, hydrogen=False, charge=False):
    global _previous_atom
    lowercase = re.compile('([a-z])')
    if re.match(lowercase, token):
        token = token.capitalize()
        return molecule_graph.add_aromatic_atom(token, isotope, hydrogen, charge)
    else:
        return molecule_graph.add_atom(token, isotope, hydrogen, charge)

def square_atom(index):
    global _tokens
    isotope = False
    if re.match('[0-9]', _tokens[index+1]):
        isotope = _tokens[index+1]
        if re.match('[0-9]', _tokens[index+2]):
            isotope = _tokens[index+1]+_tokens[index+2]
            index =+ 2
        else:
            index =+ 1
        isotope = int(isotope)
    atom = _tokens[index]
    index =+ 1
    hydrogen_count = False
    if re.match('H', _tokens[index]):
        if re.match(numbers, _tokens[index+1]) and re.match(numbers, _tokens[index+2]):
            hydrogen_count = int(_tokens[index+1]+_tokens[index+2])

        elif re.match(numbers, _tokens[index+1]):
            hydrogen_count = int(_tokens[index+1])
        elif _tokens[index+1] == 'H':
            hydrogen_count = 2
        else:
            hydrogen_count = 1




for token in _tokens:
    if re.match(letters, token):
        atom = create_atom(token)
        if _previous_atom:
            molecule_graph.add_single_bond(_previous_atom, atom)
        _previous_atom = atom

    elif token == '[':
        pass







	
