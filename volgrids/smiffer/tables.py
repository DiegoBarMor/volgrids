planar_prot = {
    "ARG" : "CZ NE NH1 NH2",
    "HIS" : "CD2 CE1 CG ND1 NE2",
    "PHE" : "CD1 CD2 CE1 CE2 CG CZ",
    "TRP" : "CD1 CD2 CE2 CE3 CG CH2 CZ2 CZ3 NE1",
    "TYR" : "CD1 CD2 CE1 CE2 CG CZ",
}

planar_rna = {
    "U" : "N1 C2 N3 C4 C5 C6",
    "C" : "N1 C2 N3 C4 C5 C6",
    "A" : "N1 C2 N3 C4 C5 C6 N7 C8 N9",
    "G" : "N1 C2 N3 C4 C5 C6 N7 C8 N9",
}

ww_scale = {
    "ASP" : -1.23, "GLU" : -2.02, "SER" : -0.13, "THR" : -0.14, "ASN" : -0.42,
    "GLN" : -0.58, "CYS" :  0.24, "PRO" : -0.45, "ALA" : -0.17, "VAL" : -0.07,
    "ILE" :  0.31, "LEU" :  0.56, "MET" :  0.23, "PHE" :  1.13, "TYR" :  0.94,
    "TRP" :  1.85, "HIS" : -0.96, "LYS" : -0.99, "ARG" : -0.81, "GLY" : -0.01,
    "U"   :  1.13, "C"   :  1.13, "A"   :  1.85, "G"   :  1.85, # RNA uses analogous protein values
}

nucleic_backbone_phosphate = ["O5'", "P", "OP1", "OP2", "O3'"]
nucleic_backbone_sugar = ["C1'", "C2'", "C3'", "C4'", "C5'", "O2'", "O4'"]
nucleic_bases = {
    "U" : ["N1", "C2", "N3", "C4", "C5", "C6", "O2", "O4"],
    "C" : ["N1", "C2", "N3", "C4", "C5", "C6", "O2", "N4"],
    "A" : ["N1", "C2", "N3", "C4", "C5", "C6", "N7", "C8", "N9", "N6"],
    "G" : ["N1", "C2", "N3", "C4", "C5", "C6", "N7", "C8", "N9", "N2", "O6"],
}

prot_hba = {
    "ALA": [ ### https://chemistry.stackexchange.com/questions/42085/can-an-amide-nitrogen-be-a-hydrogen-bond-acceptor
        [("C",), "O"]
    ],
    "ARG": [
        [("C",), "O"],
    ],
    "ASN": [
        [("C",), "O"],
        [("CG",), "OD1"]
    ],
    "ASP": [
        [("C",), "O"],
        [("CG",), "OD1"],
        [("CG",), "OD2"]
    ],
    "CYS": [
        [("C",), "O"],
        [("CB",), "SG"]
    ],
    "GLU": [
        [("C",), "O"],
        [("CD",), "OE1"],
        [("CD",), "OE2"]
    ],
    "GLN": [
        [("C",), "O"],
        [("CD",), "OE1"]
    ],
    "GLY": [
        [("C",), "O"]
    ],
    "HIS": [
        [("C",), "O"],
        [("CE1", "CG",), "ND1"], # pseudo-antecedent
    ],
    "ILE": [
        [("C",), "O"]
    ],
    "LEU": [
        [("C",), "O"]
    ],
    "LYS": [
        [("C",), "O"],
    ],
    "MET": [
        [("C",), "O"],
        [("CG",), "SD"]
    ],
    "PHE": [
        [("C",), "O"]
    ],
    "PRO": [
        [("C",), "O"]
    ],
    "SER": [
        [("C",), "O"],
        [("CB",), "OG"]
    ],
    "THR": [
        [("C",), "O"],
        [("CB",), "OG1"]
    ],
    "TRP": [
        [("C",), "O"],
    ],
    "TYR": [
        [("C",), "O"],
        [("CZ",), "OH"]
    ],
    "VAL": [
        [("C",), "O"]
    ]
}

prot_hbd = {
    "ALA": [
        [("CA",), "N"]
    ],
    "ARG": [
        [("CA",), "N"],
        [("CD", "CZ",), "NE"], # pseudo-antecedent
        [("CZ",), "NH1"],
        [("CZ",), "NH2"]
    ],
    "ASN": [
        [("CA",), "N"],
        [("CG",), "ND2"]
    ],
    "ASP": [
        [("CA",), "N"]
    ],
    "CYS": [
        [("CA",), "N"],
        [("CB",), "SG"]
    ],
    "GLU": [
        [("CA",), "N"]
    ],
    "GLN": [
        [("CA",), "N"],
        [("CD",), "NE2"]
    ],
    "GLY": [
        [("CA",), "N"]
    ],
    "HIS": [
        [("CA",), "N"],
        [("CD2", "CE1",), "NE2"] # pseudo-antecedent
    ],
    "ILE": [
        [("CA",), "N"]
    ],
    "LEU": [
        [("CA",), "N"]
    ],
    "LYS": [
        [("CA",), "N"],
        [("CE",), "NZ"]
    ],
    "MET": [
        [("CA",), "N"]
    ],
    "PHE": [
        [("CA",), "N"]
    ],
    "PRO": [],
    "SER": [
        [("CA",), "N"],
        [("CB",), "OG"]
    ],
    "THR": [
        [("CA",), "N"],
        [("CB",), "OG1"]
    ],
    "TRP": [
        [("CA",), "N"],
        [("CD1", "CE2",), "NE1"] # pseudo-antecedent
    ],
    "TYR": [
        [("CA",), "N"],
        [("CZ",), "OH"]
    ],
    "VAL": [
        [("CA",), "N"]
    ]
}

rna_hba = { # https://onlinelibrary.wiley.com/iucr/itc/Fa/ch22o2v0001/sec22o2o4.pdf
    "U": [
        [("C2'",), "O2'"],
        [("C3'", "P",), "O3'"], # pseudo-antecedent (special case)
        [("C1'", "C4'",), "O4'"], # pseudo-antecedent
        [("C5'", "P",), "O5'"], # pseudo-antecedent
        [("P",), "OP1"],
        [("P",), "OP2"],
        [("C2",), "O2"],
        [("C4",), "O4"]
    ],
    "C": [
        [("C2'",), "O2'"],
        [("C3'", "P",), "O3'"], # pseudo-antecedent (special case)
        [("C1'", "C4'",), "O4'"], # pseudo-antecedent
        [("C5'", "P",), "O5'"], # pseudo-antecedent
        [("P",), "OP1"],
        [("P",), "OP2"],
        [("C2",), "O2"],
        [("C2", "C4",), "N3"], # pseudo-antecedent
    ],
    "A": [
        [("C2'",), "O2'"],
        [("C3'", "P",), "O3'"], # pseudo-antecedent (special case)
        [("C1'", "C4'",), "O4'"], # pseudo-antecedent
        [("C5'", "P",), "O5'"], # pseudo-antecedent
        [("P",), "OP1"],
        [("P",), "OP2"],
        [("C2", "C6",), "N1"], # pseudo-antecedent
        [("C2", "C4",), "N3"], # pseudo-antecedent
        [("C5", "C8",), "N7"], # pseudo-antecedent
    ],
    "G": [
        [("C2'",), "O2'"],
        [("C3'", "P",), "O3'"], # pseudo-antecedent (special case)
        [("C1'", "C4'",), "O4'"], # pseudo-antecedent
        [("C5'", "P",), "O5'"], # pseudo-antecedent
        [("P",), "OP1"],
        [("P",), "OP2"],
        [("C2", "C4",), "N3"], # pseudo-antecedent
        [("C5", "C8",), "N7"], # pseudo-antecedent
        [("C6",), "O6"]
    ]
}

rna_hbd = {
    "U": [
        [("C2'",), "O2'"],
        [("C3'",), "O3'"], # (special case)
        [("C5'",), "O5'"], # (special case)
        [("C2", "C4",), "N3"] # pseudo-antecedent
    ],
    "C": [
        [("C2'",), "O2'"],
        [("C3'",), "O3'"], # (special case)
        [("C5'",), "O5'"], # (special case)
        [("C4",), "N4"]
    ],
    "A": [
        [("C2'",), "O2'"],
        [("C3'",), "O3'"], # (special case)
        [("C5'",), "O5'"], # (special case)
        [("C6",), "N6"]
    ],
    "G": [
        [("C2'",), "O2'"],
        [("C3'",), "O3'"], # (special case)
        [("C5'",), "O5'"], # (special case)
        [("C2", "C6",), "N1"], # pseudo-antecedent
        [("C2",), "N2"]
    ]
}
