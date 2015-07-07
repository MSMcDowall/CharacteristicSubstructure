import re
import molecule

class Parser(object):
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

    def __init__(self):
        self._previous_atom = None
        self._previous_bond = None   # If the previous token is a bond save token here
        self._break_points = {}      # Uses the number of the break point as a key and the vertex
        self._branch_root = []       # Used as stack with the append() and pop() methods

    def square_atom(self, token, mole):
        d = token.groupdict()

        if d['hydrogen'] is not None:
            hcount = '1'
            if d['hcount'] is not None:
                hcount = d['hcount']
        else:
            hcount = None

        if d['charge'] is not None:
            charge = d['charge']
            if d['chargecount'] is not None:
                charge = d['charge'] + d['chargecount']
        elif d['posdouble'] is not None:
            charge = '+2'
        elif d['negdouble'] is not None:
            charge = '-2'
        else:
            charge = None

        if d['aromatic'] is not None:
            atom = mole.add_aromatic_atom(d['element'], d['isotope'], hcount, charge)
            self.add_bond_to_aromatic_atom(atom, mole)
        else:
            atom = mole.add_atom(d['element'], d['isotope'], hcount, charge)
            self.add_bond(atom, mole)
        self._previous_atom = atom

    def organic(self, token, mole):
        d = token.groupdict()
        if d['oaromatic'] is not None:
            atom = mole.add_aromatic_atom(d['organic'])
            self.add_bond_to_aromatic_atom(atom, mole)
        else:
            atom = mole.add_atom(d['organic'])
            self.add_bond(atom, mole)
        self._previous_atom = atom
        print self._previous_atom

    def add_bond(self, atom, mole):
        if self._previous_atom is not None:
            if self._previous_bond is None:
                new_bond = mole.add_single_bond(self._previous_atom, atom)
            if self._previous_bond is not None:
                d = self._previous_bond.groupdict()
                if d['single'] is not None:
                    new_bond = mole.add_single_bond(self._previous_atom, atom)
                if d['double'] is not None:
                    new_bond = mole.add_double_bond(self._previous_atom, atom)
                if d['triple'] is not None:
                    new_bond = mole.add_triple_bond(self._previous_atom, atom)
                if d['quadruple'] is not None:
                    new_bond = mole.add_quadruple_bond(self._previous_atom, atom)
                if d['arom_bond'] is not None:
                    new_bond = mole.add_aromatic_bond(self._previous_atom, atom)
            _previous_bond = None
        else:
            new_bond = None
        return new_bond

    def add_bond_to_aromatic_atom(self, atom, mole):
        if self._previous_atom is molecule.AromaticAtom:
            mole.add_aromatic_bond(self._previous_atom, atom)
        elif self._previous_atom is not None and self._previous_atom is not molecule.AromaticAtom:
            mole.add_single_bond(self._previous_atom, atom)

    def ring(self, token, mole):
        number = token.groupdict()['ring']
        self._previous_atom.ring_break = True
        if number not in self._break_points:                  # The number has not been encountered yet (open ring)
            self._break_points[number] = [self._previous_atom, self._previous_bond]
        elif number in self._break_points:                    # The number has been encountered before (close ring)
            ring_atom = self._break_points[number][0]
            ring_bond = self._break_points[number][1]
            if ring_bond is None and self._previous_bond is None:    # No bond symbol has been specified
                print 'single'
                e = mole.add_single_bond(self._previous_atom, ring_atom)
                print e
            elif ring_bond is not None:                         # A bond symbol was specified at the ring opening
                self._previous_bond = ring_bond
                e = self.add_bond(ring_atom, mole)
            elif self._previous_bond is not None:                    # A bond symbol was specified at the ring closing
                e = self.add_bond(ring_atom, mole)

    def branch_start(self):
        self._branch_root.append(self._previous_atom)

    def branch_end(self):
        self._previous_atom = self._branch_root.pop()

    def bond(self, token):
        self._previous_bond = token

    def dot(self):
        self._previous_atom = None

    def parse_smiles(self, smiles):
        tokens = re.finditer(self.smiles_string_pattern, smiles)
        mol = molecule.Molecule(smiles)
        for a in tokens:
            d = a.groupdict()
            if d['organic'] is not None:
                self.organic(a, mol)
            elif d['square_atom'] is not None:
                self.square_atom(a, mol)
            elif d['bond'] is not None:
                self.bond(a)
            elif d['ring'] is not None:
                self.ring(a, mol)
            elif d['branch_start'] is not None:
                self.branch_start()
            elif d['branch_end'] is not None:
                self.branch_end()
            elif d['dot'] is not None:
                self.dot()
        return mol

if __name__ == '__main__':
    #complicated = 'O=C7N2c1ccccc1[C@@]64[C@@H]2[C@@H]3[C@@H](OC/C=C5\[C@@H]3C[C@@H]6N(CC4)C5)C7'
    mole = Parser().parse_smiles('C1SON1')
    for m in mole.vertices:
        print str(m) + ': '
        print m.position
        print mole.dictionary_string(m)
