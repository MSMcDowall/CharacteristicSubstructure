import re
import molecule

smiles_string_pattern = re.compile(r"""(?P<square_atom>
                                            (?P<open_bracket>\[)
                                                (?P<isotope>\d+)?
                                                (?P<element>[A-Z][a-z]?
                                                    |(?P<aromatic>b|c|n|o|p|s|se|as))
                                                (?P<chiral>@|@@)?
                                                (?P<hydrogen>[H])?
                                                (?P<hcount>\d+)?
                                                (?P<charge>[+]|[-])?
                                                (?P<chargecount>\d+)?
                                                (?P<posdouble>[+][+])?
                                                (?P<negdouble>[-][-])?
                                            (?P<close_bracket>\]))
                                        |(?P<organic>B|C|N|O|S|P|F|Cl|Br|I
                                            |(?P<oaromatic>b|c|n|o|s|p))
                                        |(?P<bond>
                                            (?P<single>-)|
                                            (?P<double>=)|
                                            (?P<triple>\#)|
                                            (?P<quadruple>$)|
                                            (?P<arom_bond>:))
                                        |(?P<ring>\d)
                                        |(?P<branch_start>\()
                                        |(?P<branch_end>\))
                                        |(?P<dot>\.)
                                        |(?P<cis_trans>\\|/)
                                    """, re.VERBOSE)
_previous_atom = None
_previous_bond = None   # If the previous token is a bond save token here
_break_points = {}      # Uses the number of the break point as a key and the vertex
_branch_root = []       # Used as stack with the append() and pop() methods


def square_atom(token, mole):
    global _previous_atom
    d = token.groupdict()
    hcount = None
    charge = None
    atom = None
    if d['hydrogen'] is not None:
        hcount = '1'
        if d['hcount'] is not None:
            hcount = d['hcount']
    if d['charge'] is not None:
        charge = d['charge']
        if d['chargecount'] is not None:
            charge = d['charge'] + d['chargecount']
    elif d['posdouble'] is not None:
        charge = '+2'
    elif d['negdouble'] is not None:
        charge = '-2'
    if d['aromatic'] is not None:
        atom = mole.add_aromatic_atom(d['element'], d['isotope'], hcount, charge)
        if _previous_atom is not None:
            mole.add_aromatic_bond(_previous_atom, atom)
    elif d['aromatic'] is None:
        atom = mole.add_atom(d['element'], d['isotope'], hcount, charge)
        add_bond(atom, mole)
    _previous_atom = atom

def organic(token, mole):
    global _previous_atom
    d = token.groupdict()
    if d['oaromatic'] is not None:
        atom = mole.add_aromatic_atom(d['organic'])
        if _previous_atom is not None:
            mole.add_aromatic_bond(_previous_atom, atom)
    elif d['oaromatic'] is None:
        atom = mole.add_atom(d['organic'])
        if _previous_atom is not None:
            add_bond(atom, mole)
    _previous_atom = atom
    print _previous_atom

def add_bond(atom, mole):
    global _previous_bond
    new_bond = None
    if _previous_atom is not None and _previous_bond is None:
        new_bond = mole.add_single_bond(_previous_atom, atom)
    if _previous_atom is not None and _previous_bond is not None:
        d = _previous_bond.groupdict()
        if d['single'] is not None:
            new_bond = mole.add_single_bond(_previous_atom, atom)
        if d['double'] is not None:
            new_bond = mole.add_double_bond(_previous_atom, atom)
        if d['triple'] is not None:
            new_bond = mole.add_triple_bond(_previous_atom, atom)
        if d['quadruple'] is not None:
            new_bond = mole.add_quadruple_bond(_previous_atom, atom)
        if d['arom_bond'] is not None:
            new_bond = mole.add_aromatic_bond(_previous_atom, atom)
        _previous_bond = None
    return new_bond

def ring(token, mole):
    global _break_points, _previous_atom, _previous_bond
    number = token.groupdict()['ring']
    _previous_atom.ring_break = True
    if number not in _break_points:                  # The number has not been encountered yet (open ring)
        _break_points[number] = [_previous_atom, _previous_bond]
    elif number in _break_points:                    # The number has been encountered before (close ring)
        print 'seen before'
        ring_atom = _break_points[number][0]
        ring_bond = _break_points[number][1]
        if ring_bond is None and _previous_bond is None:    # No bond symbol has been specified
            e = mole.add_single_bond(_previous_atom, ring_atom)
            print e
            print 'singley'
        elif ring_bond is not None:                         # A bond symbol was specified at the ring opening
            _previous_bond = ring_bond
            e = add_bond(ring_atom, mole)
            print e
            print 'ring bond'
        elif _previous_bond is not None:                    # A bond symbol was specified at the ring closing
            e = add_bond(ring_atom, mole)
            print e
            print 'previous'

def branch_start():
    global _branch_root
    _branch_root.append(_previous_atom)

def branch_end():
    global _previous_atom, _branch_root
    _previous_atom = _branch_root.pop()

def bond(token):
    global _previous_bond
    _previous_bond = token

def dot():
    global _previous_atom
    _previous_atom = None

def parse_smiles(smiles):
    global _previous_atom, _previous_bond, _break_points, _branch_root
    tokens = re.finditer(smiles_string_pattern, smiles)
    mol = molecule.Molecule(smiles)
    _previous_atom = None
    _previous_bond = None
    _break_points = {}
    _branch_root = []
    for a in tokens:
        d = a.groupdict()
        if d['organic'] is not None:
            organic(a, mol)
        elif d['square_atom'] is not None:
            square_atom(a, mol)
        elif d['bond'] is not None:
            bond(a)
        elif d['ring'] is not None:
            ring(a, mol)
        elif d['branch_start'] is not None:
            branch_start()
        elif d['branch_end'] is not None:
            branch_end()
        elif d['dot'] is not None:
            dot()
    return mol

if __name__ == '__main__':
    #complicated = 'O=C7N2c1ccccc1[C@@]64[C@@H]2[C@@H]3[C@@H](OC/C=C5\[C@@H]3C[C@@H]6N(CC4)C5)C7'
    mol = parse_smiles('[12C@H2--][N@@H+3]')
    for m in mol.vertices:
        print str(m) + ': '
        print mol.dictionary_string(m)
