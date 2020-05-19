import random
import functools
import collections

iupac_count_table = {
    'A': 1, 'C': 1, 'G': 1, 'T': 1,
    'R': 2, 'Y': 2, 'S': 2, 'W': 2,
    'K': 2, 'M': 2, 'B': 3, 'D': 3,
    'H': 3, 'V': 3, 'N': 4, 'U': 1
}

def is_homolog_legal(seq, homology):
    '''
    Check is homology is satisfiable
    '''
    if homology > len(seq):
        return False
    return True

def is_seq_constr_legal(seq):
    '''
    Check if all characters used is standard IUPAC code.
    On failure returns (False, a list of illegal chars used).
    '''
    global iupac_count_table
    chars    = set(seq)
    alphabet = set(iupac_count_table.keys())
    if chars <= alphabet:
        return (True, None)
    else:
        return (False, sorted(chars-alphabet))

def get_k_mer_count(sub_seq_list):
    '''
    Get the possible number of candidate sequences.
    '''
    global iupac_count_table
    product = 1
    for nt in sub_seq_list:
        product *= iupac_count_table[nt]
    return product

def compress_locs(locs, homology):
    '''
    Compress contigs in locs.
    '''
    if not locs:
        return None
    else:
        compressed = []
        x = locs[0]
        y = x+homology
        for pos in xrange(1, len(locs)):
            if locs[pos]-locs[pos-1] == 1:
                y = locs[pos]+homology
            else:
                compressed.append((x, y))
                x = locs[pos]
                y = x+homology
        compressed.append((x, y))
    return compressed

def is_seq_constr_sufficient(seq, homology, toolbox_size):
    '''
    Check if a desired size toolbox can be generated.
    On failure returns (False, a list of start:end tuples with constricted motifs).
    '''
    seq_list = list(seq)
    sufficiency_status = True
    insufficiency_locs = []
    for i in xrange(len(seq)-homology+1):
        if get_k_mer_count(sub_seq_list=seq_list[i:i+homology]) < toolbox_size:
            sufficiency_status = False
            insufficiency_locs.append(i)
    return (sufficiency_status, compress_locs(locs=insufficiency_locs, homology=homology))

def get_computable_form(struct):
    '''
    Get computable components from an RNA secondary structure.
    '''
    opened  = []
    closed  = []
    pairs   = []
    invalid = []
    for pos in xrange(len(struct)):
        pairing = struct[pos]
        if pairing == '(':
            opened.append(pos)
        elif pairing == ')':
            closed.append(pos)
            try:
                pairs.append((opened.pop(), closed.pop()))
            except:
                pass
        elif pairing not in ['x', '.']:
            invalid.append(pos)
    return pairs, opened, closed, invalid

def is_structure_valid(struct):
    '''
    Check if the given RNA structure is legal and balanced.
    On failure returns (False, three lists of indices: unclosed, unopened and invalid charas/parens)
    '''
    pairs, opened, closed, invalid = get_computable_form(struct)
    if not opened:
        opened  = None
    if not closed:
        closed  = None
    if not invalid:
        invalid = None
    if opened or closed or invalid:
        return (False, opened, closed, invalid)
    else:
        return (True, None, None, None)

def hamming(a, b):
    '''
    Get hamming distance between two strings.
    '''
    hdist = 0
    for i in xrange(len(a)):
        if a[i] != b[i]:
            hdist += 1
    return hdist

def is_contiguous(item1, item2):
    '''
    Check if two pairs are contiguous.
    '''
    if (item1[0]-item2[0] <= 1) and (item2[1]-item1[1] <= 1):
        return True
    return False

def is_structure_not_conflict(struct, homology):
    '''
    Q. Does the structure not enforce internal repeats?
    A. Stack evalation of stack of pairs.
    On failure returns (False, a list of locations with stems >= homology)
    '''
    pairs, opened, closed, invalid = get_computable_form(struct)
    if pairs:
        bad_contigs = []
        stretch     = 1
        pairs       = pairs[::-1]
        start       = pairs.pop()
        end         = None
        item1       = start
        item2       = end
        while pairs:
            end   = pairs.pop()
            item2 = end
            if is_contiguous(item1, item2):
                item1 = item2
                item2 = None
                stretch += 1
            else:
                if stretch >= homology:
                    bad_contigs.append(((item1[0], start[0]), (start[1], item1[1])))
                stretch = 1
                start = item2
                item1 = start
        if stretch >= homology:
            bad_contigs.append(((end[0], start[0]), (start[1], end[1])))
        if bad_contigs:
            return (False, bad_contigs)
    return (True, None)
    

def main():
    print 'Testing Sequence Constraint Validation Checkers ... ',
    legal_alphabet  = set('ACGTRYSWKMBDHVN')
    seq_legal = ''.join(random.sample(legal_alphabet, 1)[0] for _ in xrange(100))
    assert is_seq_constr_legal(seq=seq_legal)[0] == True

    illegal_alphabet = set('ABCDEFGHIJKLMNOPQRSTUVWXYZ')
    seq_illegal = ''.join(random.sample(illegal_alphabet, 1)[0] for _ in xrange(100))
    assert is_seq_constr_legal(seq=seq_illegal)[0] == False
    print 'OK'

    print 'Testing Sequence Constraint Sufficiency Checkers ... ',
    seq      = 'N'*20 + 'TTGACA' + 'N'*17 + 'TATAAT' + 'N'*6 + 'CCN' + 'N'*20
    homology = 11
    assert is_seq_constr_sufficient(seq, homology, toolbox_size=1000)[0] == True
    assert is_seq_constr_sufficient(seq, homology, toolbox_size=2000)[0] == False
    print 'OK'

    print 'Testing Struture Constraint Validation Checkers ...',
    struct = '....(((.xx.)))...()...(())...'
    assert is_structure_valid(struct)[0] == True
    struct = '....(((.xx.))).).()...(())...'
    assert is_structure_valid(struct)[0] == False
    struct = '.(..(((.xx.)))...()...(())...'
    assert is_structure_valid(struct)[0] == False
    struct = '....(((.xx.))).X.()...(())...'
    assert is_structure_valid(struct)[0] == False
    struct = '.B..(((.xx.))).).().(.(())...'
    assert is_structure_valid(struct)[0] == False
    print 'OK'

    print 'Testing Struture Constraint Homology Conflict Checkers ...',
    struct = '(((((((((((......)))))))))))...(((((((((((......)))))))))))'
    homology = 10
    assert is_structure_not_conflict(struct, homology)[0] == False
    struct = '(((((((((((......))))).))))))...(((((((((((......)))))))))))'
    homology = 10
    assert is_structure_not_conflict(struct, homology)[0] == False
    struct = '(((((((((((......)))).)))))))...(((((.((((((......)))))))))))'
    homology = 10
    assert is_structure_not_conflict(struct, homology)[0] == True
    print 'OK'

if __name__ == '__main__':
    main()

