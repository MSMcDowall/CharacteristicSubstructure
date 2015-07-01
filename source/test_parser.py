import re
import molecule

def square_atom(token, molecule_object):
    mole = molecule_object
    d = token.groupdict()
    isotope = None
    hcount = None
    charge = None
    if d['isotope']:
        isotope = d['isotope']
    element = d['element']
    if d['hydrogen']:
        hcount = 1
    if d['hcount']:
        hcount = d['hcount']
    if d['charge']:
        charge = d['charge']
    elif d['posdouble']:
        charge = '+2'
    elif d['negdouble']:
        charge = '-2'
    if d['aromatic']:
        mole.add_aromatic_atom(element, isotope, hcount, charge)
    else:
        mole.add_atom(element, isotope, hcount, charge)

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
                                  (?P<chiral>@+)
                                  ((?P<hydrogen>[H])?
                                  (?P<hcount>\d+)?)
                                  (((?P<charge>[+]|[-])?
                                  (?P<chargecount>\d+)?)|
                                  (?P<posdouble>[+][+])?|
                                  (?P<negdouble>[-][-])?)
                                  (?P<class>\d+)?
                                  (?P<close_bracket>\]))
                            |(?P<organic>B|C|N|O|S|P|F|Cl|Br|I|
                                (?P<oaromatic>b|c|n|o|s|p))
                            |(?P<bond>[-=#$:])
                            |(?P<ring>\d)
                            |(?P<branch_start>\()
                            |(?P<branch_end>\))
                            |(?P<break>\.))
                            |(?P<cis_trans>\\|/)""", re.X)

smiles = 'C[12bH2+2]CN=C.O(CF)'
tokens = re.finditer(smiles_pattern, smiles)
mol = molecule.Molecule(smiles)
for a in tokens:
    d = a.groupdict()
    if d['square_atom']:
        square_atom(a, mol)
        print mol.vertices
        print 'This is a square bracket'
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