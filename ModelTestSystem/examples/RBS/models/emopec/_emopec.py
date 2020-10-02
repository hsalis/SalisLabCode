
"""

EMOPEC main functions.

Copyright 2015 Michael Klausen, Mads Bonde, Morten Sommer, all rights reserved.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

Please cite:

  Mads T. Bonde, Margit Pedersen, Michael S. Klausen, Sheila I. Jensen, Tune
  Wulff, Scott Harrison, Alex T. Nielsen, Markus J. Herrgard and Morten O.A.
  Sommer, A novel algorithm based on the comprehensive characterization of the
  Shine-Dalgarno sequence for improved predictable tuning of protein
  expression, Nature Methods (2015)

"""

from __future__ import division

from emopec._sequences import SEQS, SEQS_MAX
from emopec._predicted_sequences import PREDICTED_SEQS


def spacing_penalty(rbs_calc_spacing, optimal_spacing=5,
                    push=(12.2, 2.5, 2.0, 3.0),
                    pull=(0.048, 0.24, 0.0)):
    """Calculate spacing penalty using RBS Calculator formula.

    :param int optimal_spacing: Number of nucleotides which is the opt. spacing.
    :param tuple push: list of c1, c2, c3, c4 if spacing is less than opt.
    :param tuple pull: list of c1, c2, c3 if spacing in more than opt.

    Note that spacing is in RBS calculator spacing, which is one base smaller
    than EMOPEC spacing.

        >>> spacing_penalty(5)
        0.0
        >>> spacing_penalty(3)
        1.525
        >>> spacing_penalty(7)
        0.6719999999999999

    """
    import math
    if rbs_calc_spacing == optimal_spacing:
        return 0.
    elif rbs_calc_spacing < optimal_spacing:
        ds = rbs_calc_spacing - optimal_spacing
        c1, c2, c3, c4 = push
        return c1 / (1.0 + math.exp(c2 * (ds + c3)))**c4
    else:
        ds = rbs_calc_spacing - optimal_spacing
        c1, c2, c3 = pull
        return c1 * ds**2 + c2 * ds + c3


def get_expression(sd_seq, sd_dist=6, penalty=0.235, use_predicted=True):
    """Get the expression level including penalty.

    :param str sd_seq:         SD sequence to get expression for.
    :param int sd_dist:        Length of spacing between SD and CDS.
    :param float penalty:      Penalty constant.
    :param bool use_predicted: Whether or not to use predicted sequence values.

        >>> get_expression('AGGAGA', 6)
        9.41932973140501
        >>> get_expression('AGGAGA', 4)
        4.127100076035202

    Using predicted sequences::

        >>> get_expression('AAAAAA', 6)
        1.4114860599086005
        >>> get_expression('AAAAAA', 6, use_predicted=False)
        0.0

    """
    sd_seq   = sd_seq.upper().replace('U', 'T')
    if sd_seq in SEQS:
        raw_expr = SEQS[sd_seq]
    elif use_predicted:
        raw_expr = PREDICTED_SEQS[sd_seq]
    else:
        return 0.

    #rbs_calc_spacing
    rcs = sd_dist - 1
    return 10.0**(raw_expr - spacing_penalty(rcs) * penalty)


def predict_spacing(leader, penalty=0.235, max_spacing=12):
    """Predict spacing based on a spacing penalty and EMOPEC signal.

        >>> predict_spacing('ATACAAGTCGCTTAAGGCTTGCCAACGAACCATTGCCGCC')
        ('ATACAAGTCGCTTAAGGCTTGCCAAC', 'GAACCA', 'TTGCCGCC', 1.6219460336049762)

    """
    _leader = leader.upper().replace('U', 'T')
    expr = 0
    spcg = 0
    for i in range(max_spacing):
        sd = _leader[-7-i:-1-i]
        _expr = get_expression(sd, i+1, penalty=penalty)
        if _expr > expr:
            expr = _expr
            spcg = i
    #Get the fractions
    upstream = leader[:-7-spcg]
    sd       = leader[-7-spcg:-1-spcg]
    spacing  = leader[-1-spcg:]
    return upstream, sd, spacing, expr


import codecs

IUPAC_ENCODING = {
    45: 0, 65: 1,  66: 14, 67: 8,  68: 7, 71: 4,  72: 11, 75: 6,
    77: 9, 78: 15, 82: 5,  83: 12, 84: 2, 86: 13, 87: 3,  89: 10
}


def make_library(up_seq, sd_seq, sp_seq, cd_seq, n=5, target=1.0, current=None,
                 maxdev=10, dg_seq=None):
    """Create a library based on the EMOPEC model.

    :param str up_seq: Sequence upstream of SD sequence.
    :param str sd_seq: 6-nucleotide SD sequence.
    :param str sp_seq: Spacing between SD and CDS.
    :param str cd_seq: Coding sequence.
    :param int n:      Library size
    :param target:     Target in relative expression (interval 0-1).
    :type  target:      float
    :param current:    Relative expression level to start the library (0-1).
                       Default is expression level of current sequence.
    :type  current:     float
    :param int maxdev: How many sequences to consider as new SD candidates.
    :param str dg_seq: IUPAC degenerate 5'-UTR sequence with constraints.

        >>> lib = make_library('ACAAGTCGCTTAAGGCTTGCCAAC', 'GAACCA', 'TTGCCGCC',
        ...                    'ATGAAGTTTATCATTAAATTGTTCCCGGAAATCACCATCAAAAGCC')
        >>> for vals in lib:
        ...     print('{}: {:.4f} {:.2f}kcal/mol {:.0%}'.format(*vals))
        AGACAT: 2.6361 0.30kcal/mol 40%
        TAGAGT: 3.6126 0.70kcal/mol 55%
        CGAGTG: 4.6015 2.50kcal/mol 70%
        AGAGTA: 5.6119 0.00kcal/mol 86%
        AGGAGT: 6.3893 1.20kcal/mol 98%

    With constraints::

        >>> l2 = make_library('ACAAGTCGCTTAAGGCTTGCCAAC', 'GAACCA', 'TTGCCGCC',
        ...                   'ATGAAGTTTATCATTAAATTGTTCCCGGAAATCACCATCAAAAGCC',
        ...                   dg_seq='NNNNNNNNNNVRRKWNNNNNNNN')
        >>> for vals in l2:
        ...     print('{}: {:.4f} {:.2f}kcal/mol {:.0%}'.format(*vals))
        AGAAAT: 2.6914 0.40kcal/mol 41%
        TAGAGT: 3.6126 0.70kcal/mol 55%
        TGGAAG: 4.4340 2.10kcal/mol 68%
        GAGAGT: 5.4777 1.90kcal/mol 84%
        AGGAGT: 6.3893 1.20kcal/mol 98%

    """
    spc = len(sp_seq)
    max_expr = get_expression(SEQS_MAX, spc)

    #Get the target level in fluorescence units
    target = max_expr * target

    #Get the initial expression level
    if current is None:
        current = get_expression(sd_seq, spc)
    else:
        current = max_expr * current

    #Normalize SD sequence to uppercase string
    sd_seq = sd_seq.upper().replace('U', 'T')

    #Leader cutoff
    lc = 35 - 6 - len(sp_seq) #6nt for SD
    #Unformatted prototype leader
    pl = (up_seq[-lc:] + '{}' + sp_seq + cd_seq[:35]).format

    import emopec.viennarna as rna

    #Run the RBS Calculator
    dG_mRNA = rna.mfe(pl(sd_seq))

    okseq = lambda s: True
    if dg_seq:
        sd_constraints = dg_seq[-7-spc:-1-spc].upper().replace('U', 'T')
        if len(sd_constraints) < 6:
            sd_constraints = (6 - len(sd_constraints)) * 'N' + sd_constraints
        enc_constraints, _ = codecs.charmap_encode(sd_constraints, 'replace', IUPAC_ENCODING)
        enc_constraints = bytearray(enc_constraints)
        def okseq(s):
            enc_sd, _ = codecs.charmap_encode(s, 'replace', IUPAC_ENCODING)
            return all(a & b for a, b in zip(bytearray(enc_sd), enc_constraints))

    all_seqs = [(s, get_expression(s, spc)) for s in SEQS if okseq(s)]

    stepsize = (target - current) / n
    steps = sorted([target - stepsize * i for i in range(n)])
    j = 0
    lib = []
    used = {sd_seq, }
    for ct in steps:
        candidates = sorted(all_seqs, key=lambda s: abs(s[1] - ct))[:maxdev]

        best = (None, None, 99999, 0.)
        for cand_sd, expr in candidates: #all_seqs[j:j+maxdev]:
            #Make sure there are no dupes in the library
            if cand_sd in used:
                continue

            dG_mRNA_mut = rna.mfe(pl(cand_sd))

            #Minimize delta delta G mRNA
            ddG = abs(dG_mRNA_mut - dG_mRNA)
            if not best[0] or ddG < best[2]:
                expr_pct = expr / max_expr
                best = (cand_sd, expr, ddG, expr_pct)

        if best[0]:
            lib.append(best)
            used.add(best[0])

    return sorted(lib, key=lambda l: l[1])
