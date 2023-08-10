"""Mol2 File Parser

This script allows the user to parse a mol2 file and extract chemical, structural,
and general information about a molecule.

This tool accepts .mol2 files.

There are no dependencies or requirements to install. Only Python requires.

This file should be imported as a module for use. 

"""

import re
import operator
import collections


class Mol2Parser:

    """
    Mol2Parser is a tool to extract information from mol2 files.
    
    ...
    
    Attributes
    ----------
    ...
        
    Methods
    -------
    _section1_extract(pattern=r"(@<TRIPOS>MOLECULE).*(@<TRIPOS>ATOM)")
        Extract information from mol2 files.
        
    _section2_extract(pattern=r"(@<TRIPOS>ATOM).*(@<TRIPOS>BOND)")
        Extract information from mol2 files.
        
    _section3_extract(pattern=r"(@<TRIPOS>BOND).*")
        Extract information from mol2 files.
        
    get_information(kind="name")
        Print general information about molecule
        
    get_molecule(kind="atom_name")
        Return general information about molecule's atoms
        
    get_bond(kind='atoms_bond')
        Return general information about molecule's bonds
        
    parse()
    
        Start parsing the mol2 file.
    """

    def __init__(self, mol2_file: str) -> None:

        """
        Parameters
        ----------
        mol2_file: str
            Name of the mol2 file as a string. 
        """

        self.mol2_file = mol2_file

    def _section1_extract(
        self, mol2_str: str, pattern: str = r"(@<TRIPOS>MOLECULE).*(@<TRIPOS>ATOM)"
    ) -> dict:

        """
        Extract a selection of string from mol2 files (private method). 
        
        Parameters
        ----------
        mol2_str: str
            Content of the opened mol2 file as a string.
            
        pattern: str
            Pattern to match to mol2_str (default is "(@<TRIPOS>MOLECULE).*(@<TRIPOS>ATOM)")
        
        Return
        ------
        
        informations: dict
            Return a dictionary of all extracted information
        """

        matches = re.finditer(pattern, mol2_str, re.MULTILINE | re.DOTALL)
        match_text = next(matches).group()
        lines = match_text.splitlines()

        info_list = ["name", "general", "type", "charge", "status_bits", "comment"]
        self.informations = dict(zip(info_list, lines[1:-1]))
        self.informations["general"] = dict(
            zip(
                [
                    "atoms_count",
                    "bonds_count",
                    "substructure_counts",
                    "features_count",
                    "sets_count",
                ],
                self.informations["general"].split(),
            )
        )

        return self.informations

    def _section2_extract(
        self, mol2_str: str, pattern: str = r"(@<TRIPOS>ATOM).*(@<TRIPOS>BOND)"
    ) -> dict:

        """
        Extract a selection of string from mol2 files (private method). 
        
        Parameters
        ----------
        mol2_str: str
            Content of the opened mol2 file as a string.
            
        pattern: str
            Pattern to match to mol2_str (default is "(@<TRIPOS>ATOM).*(@<TRIPOS>BOND)")
        
        Return
        ------
        
        informations: dict
            Return a dictionary of all extracted information
        """

        matches = re.finditer(pattern, mol2_str, re.MULTILINE | re.DOTALL)
        match_text = next(matches).group()
        lines = match_text.splitlines()

        self.index_info = collections.defaultdict(list)
        atom_name_info = {}
        coords_info = {}
        atom_type_info = {}
        subset_id_info = {}
        subset_name_info = {}
        charge_info = {}

        for line in lines[1:-1]:

            line = line.split()
            index = line[0]

            self.index_info[line[1][0]].append(index)
            atom_name_info[index] = line[1]
            coords_info[index] = line[2:5]
            atom_type_info[index] = line[5]
            subset_id_info[index] = line[6]
            subset_name_info[index] = line[7]
            charge_info[index] = line[8]

        self.molecule_info = dict(
            zip(
                [
                    "element",
                    "atom_name",
                    "coords",
                    "atom_type",
                    "subset_id",
                    "subset_name",
                    "charge",
                ],
                [
                    self.index_info,
                    atom_name_info,
                    coords_info,
                    atom_type_info,
                    subset_id_info,
                    subset_name_info,
                    charge_info,
                ],
            )
        )
        return self.molecule_info

    def _section3_extract(
        self, mol2_str: str, pattern: str = r"(@<TRIPOS>BOND).*"
    ) -> dict:

        """
        Extract a selection of string from mol2 files (private method). 
        
        Parameters
        ----------
        mol2_str: str
            Content of the opened mol2 file as a string.
            
        pattern: str
            Pattern to match to mol2_str (default is "(@<TRIPOS>BOND).*")
        
        Return
        ------
        
        informations: dict
            Return a dictionary of all extracted information
        """

        matches = re.finditer(pattern, mol2_str, re.MULTILINE | re.DOTALL)
        match_text = next(matches).group()
        lines = match_text.splitlines()

        atom_bond_info = collections.defaultdict(list)
        bond_type_info = collections.defaultdict(list)
        self.bonding_info = {"atoms_bond": atom_bond_info, "bonds_type": bond_type_info}

        for line in lines[1:]:

            line = line.split()
            if line[0] == "@<TRIPOS>SUBSTRUCTURE":
                break
            atom_bond_info[line[1]].append(line[2])
            bond_type_info[line[1]].append(line[3])

        return self.bonding_info

    def get_information(self, kind: str = "name") -> str:

        """
        Return general information about the molecule: name, charge, type, 
        general structure, status_bits, and comment
        
        Parameter
        ---------
        
        kind: str
            
            Select desired variable to display (default is "name"). Possible value:
            "name", "general", "type", "charge", "status_bits", "comment".
            
        Return
        ------
        variable: str
            The value of the selected kind.
            
        Raises
        ------
        ValueError
            If selected kind is not in the possible values.
            
        """

        valid_kinds = ["name", "general", "type", "charge", "status_bits", "comment"]

        if kind not in valid_kinds:

            raise ValueError("Wrong kind")

        else:

            variable = self.informations.get(kind, "Missed or not mentioned Value")
            return variable

    def get_molecule(self, element: str, kind: str = "atom_name") -> dict:

        """
        Return chemical and structural information about the molecule: atom name,
        atom coordinates, atom type, subset id, subset name, atomic charge
        
        Parameter
        ---------
        element: str
            Element symbol. Like "C".
            
        kind: str
            
            Select desired variable to display (default is "atom_name"). Possible value:
            'atom_name', 'coords', 'atom_type', 'subset_id', 'subset_name', 'charge'.
            
        Return
        ------
        variable: dict
            The value of the selected kind for the desired element in the dictionary.
            
        Raises
        ------
        ValueError
            If selected kind is not in the possible values.
            
        """

        valid_keys = [
            "atom_name",
            "coords",
            "atom_type",
            "subset_id",
            "subset_name",
            "charge",
        ]

        if kind not in valid_keys:

            raise ValueError("Wrong kind")

        else:

            f = operator.itemgetter(*self.index_info[element])
            variable = list(f(self.molecule_info[kind]))
            return variable

    def get_bond(self, kind: str = "atoms_bond") -> dict:

        """
        Return bonding information about the molecule: connected atoms and their
        bond order
        
        Parameter
        ---------
        kind: str
            
            Select desired variable to display (default is "atomes_bond"). Possible value:
            'atoms_bond', 'bonds_type'.
            
        Return
        ------
        variable: dict
            The value of the selected kind in the dictionary.
            
        Raises
        ------
        ValueError
            If selected kind is not in the possible values.
            
        """

        valid_keys = ["atoms_bond", "bonds_type"]

        if kind not in valid_keys:

            raise ValueError("Wrong kind")

        else:
            variable = self.bonding_info[kind]
            return variable

    def parse(self) -> None:

        """
        Parsing the mol2 file to extract required information.
        
        """

        with open(self.mol2_file) as file:
            self.mol2_str = file.read()

        self.information = self._section1_extract(self.mol2_str)
        self.molecule_info = self._section2_extract(self.mol2_str)
        self.bonding_info = self._section3_extract(self.mol2_str)
