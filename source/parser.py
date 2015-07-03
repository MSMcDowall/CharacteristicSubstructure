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
_break_points = {}      # Uses the number of the break point as a key and the vertex
_branch_root = []       # Used as stack with the append() and pop() methods
_previous_bond = None  # Change to weight of bond if the previous token is a bond

def square_atom(token, mole):
    global _previous_atom
    d = token.groupdict()
    hcount = None
    charge = None
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
        atom = mole.add_aromatic_atom(d['element'], d['isotope'], hcount, charge)
        if _previous_atom:
            mole.add_aromatic_bond(_previous_atom, atom)
    else:
        atom = mole.add_atom(d['element'], d['isotope'], hcount, charge)
        add_bond(atom, mole)
    _previous_atom = atom

def organic(token, mole):
    global _previous_atom
    d = token.groupdict()
    if d['oaromatic']:
        atom = mole.add_aromatic_atom(d['organic'])
        if _previous_atom:
            mole.add_aromatic_bond(_previous_atom, atom)
    else:
        atom = mole.add_atom(d['organic'])
        add_bond(atom, mole)
    _previous_atom = atom

def add_bond(atom, mole):
    global _previous_bond
    if _previous_atom and not _previous_bond:
        mole.add_single_bond(_previous_atom, atom)
    if _previous_atom and _previous_bond:
        d = _previous_bond.groupdict()
        if d['single']:
            mole.add_single_bond(_previous_atom, atom)
        if d['double']:
            mole.add_double_bond(_previous_atom, atom)
        if d['triple']:
            mole.add_triple_bond(_previous_atom, atom)
        if d['quadruple']:
            mole.add_quadruple_bond(_previous_atom, atom)
        if d['arom_bond']:
            mole.add_aromatic_bond(_previous_atom, atom)
        _previous_bond = None

def ring(token, mole):
    global _break_points, _previous_atom, _previous_bond
    _previous_atom.ring_break = True
    if token not in _break_points.keys():           # The number has not been encountered yet (open ring)
        _break_points['token'] = [_previous_atom, _previous_bond]
    elif token in _break_points:                    # The number has been encountered before (close ring)
        ring_atom = _break_points['token'][0]
        ring_bond = _break_points['token'][1]
        if not ring_bond and not _previous_bond:    # No bond symbol has been specified
            mole.add_single_bond(_previous_atom, ring_atom)
        elif ring_bond:                             # A bond symbol was specified at the ring opening
            _previous_bond = ring_bond
            add_bond(ring_bond, mole)
        elif _previous_bond:                        # A bond symbol was specified at the ring closing
            add_bond(ring_bond,mole)

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
    tokens = re.finditer(smiles_string_pattern, smiles)
    mol = molecule.Molecule(smiles)
    for a in tokens:
        d = a.groupdict()
        if d['organic']:
            organic(a, mol)
        elif d['square_atom']:
            square_atom(a, mol)
        elif d['bond']:
            bond(a)
        elif d['ring']:
            ring(a, mol)
        elif d['branch_start']:
            branch_start()
        elif d['branch_end']:
            branch_end()
        elif d['dot']:
            dot()
    return mol
complicated = 'O=C7N2c1ccccc1[C@@]64[C@@H]2[C@@H]3[C@@H](OC/C=C5\[C@@H]3C[C@@H]6N(CC4)C5)C7'
mol = parse_smiles('')
for m in mol.vertices:
    print str(m) + ': '
    print mol.dictionary_string(m)
