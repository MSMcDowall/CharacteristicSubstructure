import re
import molecule

def square_atom(scanner, token):
    return token

def organic(scanner, token):
    return token

def ring(scanner, token):
    return token

def branch_start(scanner, token):
    return token

def branch_end(scanner, token):
    return token

def bond(scanner, token):
    return token

smiles_pattern = re.compile(r"""((?P<square_atom>
                                  (?P<open_bracket>\[)
                                  (?P<isotope>\d+)?
                                  (?P<element>[A-Z][a-z]?
                                        |(?P<aromatic>b|c|n|o|p|s|se|as))
                                  (?P<chiral>@|@@)
                                  (?P<hydrogen>[H])?
                                  (?P<hcount>\d+)?
                                  (?P<charge>[+]|[-])?
                                  (?P<chargecount>\d+)?
                                  (?P<posdouble>[+][+])?
                                  (?P<negdouble>[-][-])?
                                  (?P<close_bracket>\]))
                            |(?P<organic>B|C|N|O|S|P|F|Cl|Br|I|
                                (?P<oaromatic>b|c|n|o|s|p))
                            |(?P<bond>[-=#$:])
                            |(?P<ring>\d)
                            |(?P<branch_start>\()
                            |(?P<branch_end>\))
                            |(?P<break>\.))""", re.X)
smiles = '[12bH2+2]CN=C.O(CF)'
tokens = re.finditer(smiles_pattern, smiles)
molecule = molecule.Molecule(smiles)
for a in tokens:
    d = a.groupdict()
    if d['square_atom']:
        print 'This is a square bracket'
        if d['aromatic']:
            print 'This is aromatic'
    elif d['organic']:
        print 'This is organic'
    elif d['bond']:
        print 'This is a bond'
    elif d['ring']:
        print 'This is a ring'
    elif d['branch_start']:
        print 'This it the branch start'
    elif d['branch_end']:
        print 'This is the branch end'
    elif d['break']:
        print 'This is a break'