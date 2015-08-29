import re
from molecule import Molecule


class Parser(object):
    """
    Contains the functions and variables required to parse a SMILES string into a molecule object

    Follows the OpenSMILES specification [opensmiles.org]
    Regular expressions modified from Frowns smiles_parsers found at frowns.sourceforge.net (Brian Kelley 2001-2002)
    """
    # Regular expression detailing all chemical symbols
    element_symbols = r"C[laroudsemf]?|Os?|N[eaibdpos]?|S[icernbmg]?|P[drmtboau]?|"  \
                      r"H[eofgas]?|c|n|o|s|p|A[lrsgutcm]|B[eraik]?|Dy|E[urs]|F[erm]?|"  \
                      r"G[aed]|I[nr]?|Kr?|L[iaur]|M[gnodt]|R[buhenaf]|T[icebmalh]|" \
                      r"U|V|W|Xe|Yb?|Z[nr]|\*"

    # Regular expression which splits SMILES string into groups
    smiles_string_pattern = re.compile(r"""(?P<square_atom>
                                            (?P<open_bracket>\[)
                                                (?P<isotope>\d+)?
                                                (?P<element>""" + element_symbols + r"""
                                                    |(?P<aromatic>b|c|n|o|p|s|se|as))
                                                (?P<chiral>@|@@)?
                                                (?P<hydrogen>[H])?
                                                (?P<hcount>\d+)?
                                                (?P<charge>[+]|[-])?
                                                (?P<chargecount>\d+)?
                                                (?P<posdouble>[+][+])?
                                                (?P<negdouble>[-][-])?
                                            (?P<close_bracket>\]))
                                        |(?P<organic>Br|Cl|B|C|N|O|S|P|F|I
                                            |(?P<oaromatic>b|c|n|o|s|p))
                                        |(?P<bond>
                                            (?P<single>-)|
                                            (?P<double>=)|
                                            (?P<triple>\#)|
                                            (?P<quadruple>\$)|
                                            (?P<arom_bond>\:))
                                        |(?P<ring>\d)
                                        |(?P<branch_start>\()
                                        |(?P<branch_end>\))
                                        |(?P<dot>\.)
                                        |(?P<cis_trans>\\|/)
                                    """, re.VERBOSE)

    def __init__(self):
        self._previous_atom = None   # Holds the previous atom object which was encountered
        self._previous_bond = None   # If the previous token is a bond save token here
        self._break_points = {}      # Uses the number of the break point as a key and the vertex
        self._branch_root = []       # Used as stack with the append() and pop() methods

    def parse_smiles(self, smiles):
        """
        Takes in a SMILES string and creates a molecule object

        :param smiles: the SMILES string that is to be parsed
        :return: a molecule object which has the structured described in the SMILES string
        """
        tokens = re.finditer(self.smiles_string_pattern, smiles)
        molecule = Molecule(smiles)
        for a in tokens:
            d = a.groupdict()
            if d['organic'] is not None:
                self._organic(a, molecule)
            elif d['square_atom'] is not None:
                self._square_atom(a, molecule)
            elif d['bond'] is not None:
                self._bond(a)
            elif d['ring'] is not None:
                self._ring(a, molecule)
            elif d['branch_start'] is not None:
                self._branch_start()
            elif d['branch_end'] is not None:
                self._branch_end()
            elif d['dot'] is not None:
                self._dot()
        return molecule

    def _square_atom(self, token, molecule):
        """
        Adds an atom to the molecule which has specified hydrogen atoms, charge and weight (isotope)

        :param token: the token from the regular expression
        :param molecule: the molecule object which is being created
        :return: None
        """
        d = token.groupdict()

        # Checks if an explicit hydrogen (and number of hydrogens) has been specified
        if d['hydrogen'] is not None:
            hcount = '1'
            if d['hcount'] is not None:
                hcount = d['hcount']
        else:
            hcount = None

        # Checks if an explicit charge has been specified
        if d['charge'] is not None:
            charge = d['charge']
            if d['chargecount'] is not None:
                charge = d['charge'] + d['chargecount']
        # Accepts the now depreciated ++ and -- to signify +2 and -2 respectively
        elif d['posdouble'] is not None:
            charge = '+2'
        elif d['negdouble'] is not None:
            charge = '-2'
        else:
            charge = None

        # Add atom including attributes to molecule
        # Call the method to add a bond between this atom and the previous atom
        # The bond method chosen depends on if the atom is aromatic (lowercase) or not
        if d['aromatic'] is not None:
            atom = molecule.add_aromatic_atom(d['element'], d['isotope'], hcount, charge)
            self._add_bond_to_aromatic_atom(atom, molecule)
        else:
            atom = molecule.add_atom(d['element'], d['isotope'], hcount, charge)
            self._add_bond(atom, molecule)
        self._previous_atom = atom

    def _organic(self, token, molecule):
        """
        Adds an organic atom [B,C,N,O,S,P,F,C,Br,I] which has implied hydrogen atoms and standard charge and weight

        :param token: the token from the regular expression
        :param molecule: the molecule object which is being created
        :return: None
        """
        d = token.groupdict()
        # Add atom to molecule and call the method to add a bond between this atom and the previous atom
        # The bond method chosen depends on if the atom is aromatic (lowercase) or not
        if d['oaromatic'] is not None:
            atom = molecule.add_aromatic_atom(d['organic'])
            self._add_bond_to_aromatic_atom(atom, molecule)
        else:
            atom = molecule.add_atom(d['organic'])
            self._add_bond(atom, molecule)
        # Before moving to the next part of the SMILES.txt string set this atom as the previous atom
        self._previous_atom = atom

    def _add_bond(self, atom, molecule):
        """
        Adds a bond to the molecule and if a bond symbol has been encountered it tests the type

        :param token: the token from the regular expression
        :param molecule: the molecule object which is being created
        :return: None
        """
        if self._previous_atom is not None:
            if self._previous_bond is None:
                new_bond = molecule.add_single_bond(self._previous_atom, atom)
            if self._previous_bond is not None:
                d = self._previous_bond.groupdict()
                if d['single'] is not None:
                    new_bond = molecule.add_single_bond(self._previous_atom, atom)
                if d['double'] is not None:
                    new_bond = molecule.add_double_bond(self._previous_atom, atom)
                if d['triple'] is not None:
                    new_bond = molecule.add_triple_bond(self._previous_atom, atom)
                if d['quadruple'] is not None:
                    new_bond = molecule.add_quadruple_bond(self._previous_atom, atom)
                if d['arom_bond'] is not None:
                    new_bond = molecule.add_aromatic_bond(self._previous_atom, atom)
            self._previous_bond = None
        else:
            new_bond = None
        return new_bond

    def _add_bond_to_aromatic_atom(self, atom, molecule):
        """
        Tests if a single or an aromatic bond should be added to an aromatic atom

        :param token: the token from the regular expression
        :param molecule: the molecule object which is being created
        :return: None
        """
        if self._previous_atom is not None:
            if self._previous_atom.aromatic:
                molecule.add_aromatic_bond(self._previous_atom, atom)
            elif self._previous_atom is not None and not self._previous_atom.aromatic:
                molecule.add_single_bond(self._previous_atom, atom)

    def _ring(self, token, molecule):
        """
        Called when a number is encountered that is not in square brackets

        :param token: the token from the regular expression
        :param molecule: the molecule object which is being created
        :return: None
        """
        number = token.groupdict()['ring']
        self._previous_atom.ring_break = True
        if number not in self._break_points:                  # The number has not been encountered yet (open ring)
            self._break_points[number] = [self._previous_atom, self._previous_bond]
        elif number in self._break_points:                    # The number has been encountered before (close ring)
            ring_atom = self._break_points[number][0]
            ring_bond = self._break_points[number][1]
            if ring_bond is None and self._previous_bond is None:    # No bond symbol has been specified
                if ring_atom.aromatic:
                    self._add_bond_to_aromatic_atom(ring_atom, molecule)
                else:
                    molecule.add_single_bond(self._previous_atom, ring_atom)
            elif ring_bond is not None:                         # A bond symbol was specified at the ring opening
                self._previous_bond = ring_bond
                e = self._add_bond(ring_atom, molecule)
            elif self._previous_bond is not None:                    # A bond symbol was specified at the ring closing
                e = self._add_bond(ring_atom, molecule)

    def _branch_start(self):
        """
        Called when opening parentheses are encountered so a branch has been started

        When a branch is started the previous atom is stored as branch root
        Using a stack ensures that nested branches are closed in the correct order
        :return: None
        """
        self._branch_root.append(self._previous_atom)

    def _branch_end(self):
        """
        Called when closing parentheses are encountered so a branch has ended

        After a branch return back to the root of the branch
        :return: None
        """
        self._previous_atom = self._branch_root.pop()

    def _bond(self, token):
        """
        Called when a bond symbol is encountered

        When a bond symbol is encountered the match token is stored in the previous bond
        Once the subsequent atom has been made a bond can be made depending on the properties of the token
        :param token: the token from the regular expression
        :return: None
        """
        self._previous_bond = token

    def _dot(self):
        """
        Called when a dot is encountered which means that there is no bond between the atoms on either side of it

        By removing the reference to the previous atom it ensures no bond is created
        :return: None
        """
        self._previous_atom = None

if __name__ == '__main__':
    Parser().parse_smiles('C1CCC1')
