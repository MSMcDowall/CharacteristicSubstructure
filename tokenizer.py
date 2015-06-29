import re


def tokenize_SMILES(molecule_string):
	#probelem getting the regular expression to check the aromatic set first and then [A-Z][a-z]?
	# |\d|b|c|n|o|p|s|se|as|[A-Z][a-z]?
	# Regular expressions to break the SMILES strings component parts
	seperate_letters = re.compile('([-=#\$:/\\ \(\).\[\]]|[A-Za-z]|[0-9])')
	aromatic = re.compile('b|c|n|o|p|s|se|as')
	lowercase = re.compile('[a-z]')
	# Seperate the molecule string into individual letters
	letters = seperate_letters.findall(molecule_string)
	print letters
	# concatenate capital letter with small letter (not aromatic) to make atomic symbol
	tokens = []
	for token in letters:
		if re.match(lowercase, token) and not re.match(aromatic, token):
			position = letters.index(token)
			tokens.append(letters[position-1]+letters[position])
		else:
			tokens.append(token)
	return tokens
	
string = 'C12CN[12C]He'
print tokenize_SMILES(string)
