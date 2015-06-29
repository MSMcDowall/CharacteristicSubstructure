import molecule
import tokenizer
import re


SMILES_string = ''
_previous_atom = None
_break_points = {}
#_branch_root = #stack

# Read in string from file
_tokens = tokenizer.tokenize_SMILES(SMILES_string)
_graph = molecule.Molecule(SMILES_string)


def create_atom(token, isotope=False, hydrogen=False, charge=False):
    global _previous_atom
    aromatic = False
    if re.match('[a-z][a-z]?', token):
        vertex = _graph.add_aromatic_atom(token)
        # capitalize first letter of token
    else:
        vertex = _graph.add_vertex(token)
    if _previous_atom:
        _graph.add_edge(_previous_atom, vertex)
    _previous_atom = vertex

def square_atom(index):
    isotope = False
    if re.match('[0-9]', tokens[index+1]):
        isotope = _tokens[index+1]
        if re.match('[0-9]', _tokens[index+2]):
            isotope = _tokens[index+1]+_tokens[index+2]
            index = index + 2
        else:
            index =+ 1
        # Change isotope to integer
    atom = _tokens[index]
    if re.match('H', tokens[index+1]
		
	


	
