
from weblogolib import *
from weblogolib.colorscheme import *

from corebio.seq import (Seq,SeqList)
from corebio.seq import (Alphabet, Seq, SeqList, unambiguous_dna_alphabet,
                         unambiguous_rna_alphabet, unambiguous_protein_alphabet,
                         dna_alphabet,rna_alphabet)

rna_custom = Alphabet("ACGUX", zip('acgux', 'ACGUX'))

nucleotide = ColorScheme(
        [
            SymbolColor("G", "orange"),
            SymbolColor("TU", "red"),
            SymbolColor("C", "blue"),
            SymbolColor("A", "green"),
            SymbolColor("X", "black")
        ],
)

def dna(string):
    """Create an alphabetic sequence representing a stretch of DNA."""
    return Seq(string, alphabet=dna_alphabet)

def rna(string):
    """Create an alphabetic sequence representing a stretch of RNA."""
    return Seq(string, alphabet=rna_custom)

def create_logo(sequences,filename):

    # assume RNA for now
    sequences = SeqList([rna(seq) for seq in sequences],alphabet=rna_custom)

    # clLogoData = LogoData(length=len(sequences[0]),alphabet=rna_custom)
    data = LogoData.from_seqs(sequences)

    # specify changes to the logo
    options = LogoOptions()
    options.color_scheme = nucleotide
    options.stacks_per_line = 1000

    # specify changes to the format
    frmt = LogoFormat(data,options)

    # create eps file and save
    eps = eps_formatter(data,frmt)
    file = open(filename, "w")
    file.write(eps)
    file.close()

def test_logo():

    sequences = [
    "AUGCUAGUUUCGCXAGUAAACGGG",
    "CGGGGAUUUAUGCXUUUAAUACUG",
    "ACUGCGAGAGGGCXGAUUUUAGCG",
    "AUGGGGGGGGGAUXUUUUUUAUGG",
    "AUGGCUAUUUUCGXGAUUAAAGCG"
    ]

    create_logo(sequences,filename='test.eps')


if __name__ == "__main__":
    pass