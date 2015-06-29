import re


def tokenize_smiles(molecule_string):
    # Regular expressions to break the SMILES strings component parts
    separate_symbols = re.compile('([-=#\$:/\\ \(\).\[\]]|[A-Za-z][a-z]?|[0-9])')
    # Separate the molecule string into individual letters
    symbols = separate_symbols.findall(molecule_string)

    # Separate any lowercase letters which are aromatic to make atomic symbol
    tokens = []
    aromatic_pattern = re.compile('b|c|n|o|p|s|se|as')
    letters_pattern = re.compile('([A-Za-z])')
    for token in symbols:
        if len(token)>1 and re.match(aromatic_pattern, token[1]):
            separated = re.findall(letters_pattern, token)
            tokens.append(separated[0])
            tokens.append(separated[1])
        else:
            tokens.append(token)
    return tokens

string = 'C12CN=[12C]He1C(c)Oo'
print tokenize_smiles(string)
