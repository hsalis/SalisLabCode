from string       import maketrans
from Bio.SeqUtils import MeltingTemp
from time         import time

class Synthesis(object):

    # Setup synthesis assessment variables
    comp_table = maketrans('ATGC', 'TACG')
    runs_tuple = (('CCCCCCCCC',            'GGGGGGGGG'),
                  ('AAAAAAAAAAAAA',        'TTTTTTTTTTTTT'),
                  ('TCTTCTTCTTCTTCTTCT',   'TCGTCGTCGTCGTCGTCG',   'TCCTCCTCCTCCTCCTCC',   'GGAGGAGGAGGAGGAGGA',   'GGTGGTGGTGGTGGTGGT',
                   'GGCGGCGGCGGCGGCGGC',   'AATAATAATAATAATAAT',   'AAGAAGAAGAAGAAGAAG',   'AACAACAACAACAACAAC',   'CCACCACCACCACCACCA',
                   'CCTCCTCCTCCTCCTCCT',   'CCGCCGCCGCCGCCGCCG',   'ATAATAATAATAATAATA',   'ATTATTATTATTATTATT',   'ATGATGATGATGATGATG',
                   'ATCATCATCATCATCATC',   'AGAAGAAGAAGAAGAAGA',   'AGTAGTAGTAGTAGTAGT',   'AGGAGGAGGAGGAGGAGG',   'AGCAGCAGCAGCAGCAGC',
                   'ACAACAACAACAACAACA',   'ACTACTACTACTACTACT',   'ACGACGACGACGACGACG',   'ACCACCACCACCACCACC',   'TAATAATAATAATAATAA',
                   'TATTATTATTATTATTAT',   'TAGTAGTAGTAGTAGTAG',   'TACTACTACTACTACTAC',   'TTATTATTATTATTATTA',   'TTGTTGTTGTTGTTGTTG',
                   'TTCTTCTTCTTCTTCTTC',   'TGATGATGATGATGATGA',   'TGTTGTTGTTGTTGTTGT',   'TGGTGGTGGTGGTGGTGG',   'TGCTGCTGCTGCTGCTGC',
                   'TCATCATCATCATCATCA',   'GAAGAAGAAGAAGAAGAA',   'GATGATGATGATGATGAT',   'GAGGAGGAGGAGGAGGAG',   'GACGACGACGACGACGAC',
                   'GTAGTAGTAGTAGTAGTA',   'GTTGTTGTTGTTGTTGTT',   'GTGGTGGTGGTGGTGGTG',   'GTCGTCGTCGTCGTCGTC',   'GCAGCAGCAGCAGCAGCA',
                   'GCTGCTGCTGCTGCTGCT',   'GCGGCGGCGGCGGCGGCG',   'GCCGCCGCCGCCGCCGCC',   'CAACAACAACAACAACAA',   'CATCATCATCATCATCAT',
                   'CAGCAGCAGCAGCAGCAG',   'CACCACCACCACCACCAC',   'CTACTACTACTACTACTA',   'CTTCTTCTTCTTCTTCTT',   'CTGCTGCTGCTGCTGCTG',
                   'CTCCTCCTCCTCCTCCTC',   'CGACGACGACGACGACGA',   'CGTCGTCGTCGTCGTCGT',   'CGGCGGCGGCGGCGGCGG',   'CGCCGCCGCCGCCGCCGC'),
                  ('ATATATATATATATATATAT', 'AGAGAGAGAGAGAGAGAGAG', 'ACACACACACACACACACAC', 'TATATATATATATATATATA', 'TGTGTGTGTGTGTGTGTGTG', 'TCTCTCTCTCTCTCTCTCTC',
                   'GAGAGAGAGAGAGAGAGAGA', 'GTGTGTGTGTGTGTGTGTGT', 'GCGCGCGCGCGCGCGCGCGC', 'CACACACACACACACACACA', 'CTCTCTCTCTCTCTCTCTCT', 'CGCGCGCGCGCGCGCGCGCG'))

    def _is_GC_content_pass(self, seq, start, end, win, gc_low, gc_high):
        """
        >>> synth_obj._is_GC_content_pass(seq='AAAAAAAAAAAAAAAAAAGG', start=0, end=20, win=20, gc_low=0.10, gc_high=0.25)
        False
        >>> synth_obj._is_GC_content_pass(seq='AAAAAAAAAAAAAAAGGGGG', start=0, end=20, win=20, gc_low=0.10, gc_high=0.25)
        False
        >>> synth_obj._is_GC_content_pass(seq='AAAAAAAAAAAAAAAAGGGG', start=0, end=20, win=20, gc_low=0.10, gc_high=0.25)
        True
        """

        gc_count, i = 0.0, start
        while i < end:
            gc_count += seq[i] in ['G', 'C']
            if (gc_count / win) >= gc_high:
                return False
            if i >= (start + win - 1):
                if (gc_count / win) <= gc_low:
                    return False
                gc_count -= seq[i-win+1] in ['G', 'C']
            i += 1
        return True

    def _rule_GC_content(self, seq):
        """        
        >>> synth_obj._rule_GC_content(seq='GGGGGGGGGGGGGGAAAAAA')
        False
        >>> synth_obj._rule_GC_content(seq='GGGGGGAAAAAAAAAAAAAA')
        False
        >>> synth_obj._rule_GC_content(seq='GGGGGGGGGGAAAAAAAAAA')
        True
        >>> synth_obj._rule_GC_content(seq='GAAAAAAAAAAAAAAAAAGGGGGGGGGGGG')
        False
        >>> synth_obj._rule_GC_content(seq='GAAAAAAAAAAAAAAAAAGGGGGGGGGGGG')
        False
        >>> synth_obj._rule_GC_content(seq='GAGAGAGAGAGAGAGAGAGAGAGAGAGAGA')
        True
        """

        # Throw out insane sequences
        # if not self._is_GC_content_pass(seq=seq, start=0, end=len(seq), win=len(seq), gc_low=0.30, gc_high=0.68):
        #     return False
        # else:
        # Windowed GC content checkup
        if len(seq) >= 20  and not self._is_GC_content_pass(seq=seq, start=0, end=len(seq), win=20,  gc_low=0.15, gc_high=0.90):
            return False
        if len(seq) >= 100 and not self._is_GC_content_pass(seq=seq, start=0, end=len(seq), win=100, gc_low=0.28, gc_high=0.76):
            return False
        # Terminal GC content checkup
        if len(seq) > 60:
            # 5' end
            if not self._is_GC_content_pass(seq=seq, start=0, end=30, win=30, gc_low=0.24, gc_high=0.76):
                return False
            # 3' end
            if not self._is_GC_content_pass(seq=seq, start=len(seq)-30, end=len(seq), win=30, gc_low=0.24, gc_high=0.76):
                return False
        # All GC conditions satisfy
        return True

    def _is_melting_temp_pass(self, seq, tm_low, tm_high):
        """
        >>> synth_obj._is_melting_temp_pass(seq='AAAAAAAAAATTTTTTTTTT', tm_low=36.0, tm_high=73.0)
        False
        >>> synth_obj._is_melting_temp_pass(seq='GCGCGCGCGCGCATATATAT', tm_low=36.0, tm_high=73.0)
        True
        """
        
        if tm_low <= MeltingTemp.Tm_NN(seq, dnac1=250.0, dnac2=0.0, saltcorr=7) <= tm_high:
            return True
        return False

    def _rule_melting_temp(self, seq):
        """
        >>> synth_obj._rule_melting_temp(seq='AAAAAAAAAATTTTTTTTTT')
        False
        >>> synth_obj._rule_melting_temp(seq='GCGCGCGCGCGCATATATAT')
        True
        """

        i, win = 0, 20
        if len(seq) >= win:
            while i < len(seq)-win+1:
                if not self._is_melting_temp_pass(seq=seq[i:i+win], tm_low=36.0, tm_high=73.0):
                    return False
                i += 1
        return True

    def _is_slippage_pass(self, seq, start, end, rep_len, win):
        end = min([start+(rep_len*win), end])
        # Direct repeat slippage
        beg = start+1
        pattern = seq[start:start+rep_len]
        if seq.find(pattern, beg, end) > -1:
            return False
        # Inverted repeat slippage, causing a hairpin
        beg = start+rep_len
        pattern = pattern.translate(self.comp_table)[::-1]
        if seq.find(pattern, beg, end) > -1:
            return False
        return True

    def _rule_slippage(self, seq):
        start, rep_len, win, term_len = 0, 5, 3, 60
        # Check terminals of candidate sequence for tandem repeats
        if len(seq) > rep_len:
            # Too small for terminal checks, search fully
            if len(seq) <= term_len:
                start_max, end = len(seq)-rep_len, len(seq)
            # Initiate 5' terminal check
            else:                
                start_max, end = term_len-rep_len, term_len
            while start < start_max:
                if not self._is_slippage_pass(seq, start, end, rep_len, win):
                    return False
                start += 1
                # Switch to 3' terminal check
                if start == term_len-rep_len:
                    start, start_max, end = max([start, len(seq) - term_len]), len(seq)-rep_len, len(seq)
        return True

    def _is_runs_pass(self, seq, term_len):
        beg, end = term_len, len(seq)-term_len
        for run_tuple in self.runs_tuple:
            # Non-terminal region long enough to have runs
            if len(run_tuple[0]) <= len(seq)-(2*term_len):
                for run in run_tuple:
                    if seq.find(run, beg, end) > -1:
                        return False
            else:
                break
        return True

    def _rule_runs(self, seq):
        term_len = 60
        # Check for bad runs in non-terminal region only (compliments _rule_slippage)
        return self._is_runs_pass(seq, term_len)

    def _is_hairpin_pass(self, seq, stem, loop, max_mismatch, gc_high):
        i = 0
        while i < len(seq)-(stem+loop+stem)+1:
            j, mismatch, gc_count = 0, 0, 0.0
            while j < stem:
                if seq[i+j] in ['G', 'C']:
                    gc_count += 1.0
                if seq[i+(stem+loop+stem)-j-1] in ['G', 'C']:
                    gc_count += 1.0
                if seq[i+j] != seq[i+(stem+loop+stem)-j-1].translate(self.comp_table):
                    mismatch += 1
                if mismatch > max_mismatch:
                    break
                j += 1
            if j == stem:
                gc_content = (gc_count / 2.0) / stem
                if gc_content >= gc_high:
                    return False
            i += 1
        return True

    def _rule_palindrome(self, seq):
        stem, loop, max_mismatch, gc_high = 11, 0, 1, 0.0
        if len(seq) >= (stem+loop+stem):
            return self._is_hairpin_pass(seq, stem, loop, max_mismatch, gc_high)
        return True

    def _rule_hairpin(self, seq):
        # Type 1: Stem = 11bp, Loop = 3bp to 48bp, Mismatches = 2, Stem must have GC-content >= 0.80
        stem, max_mismatch, gc_high = 11, 2, 0.80
        for loop in xrange(3, 48+1):
            if len(seq) >= (stem+loop+stem):
                if not self._is_hairpin_pass(seq, stem, loop, max_mismatch, gc_high):
                    return False
            else:
                break

        # Type 2: Stem = 17bp, Loop = 3bp to 100bp, Mismatches = 3
        stem, max_mismatch, gc_high = 17, 3, 0.0
        for loop in xrange(3, 100+1):
            if len(seq) >= (stem+loop+stem):
                if not self._is_hairpin_pass(seq, stem, loop, max_mismatch, gc_high):
                    return False
            else:
                break

        if len(seq) > 500:
            # Type 3: Stem = 22bp, Loop = 100bp to 500bp, Mismatches = 5
            stem, max_mismatch, gc_high = 22, 5, 0.0
            for loop in xrange(100, 500+1):
                if len(seq) >= (stem+loop+stem):
                    if not self._is_hairpin_pass(seq, stem, loop, max_mismatch, gc_high):
                        return False
                else:
                    break

        return True

    def evaluate_short(self, seq):
        if self._rule_GC_content(seq):
            if self._rule_runs(seq):
                if self._rule_palindrome(seq):
                    if self._rule_hairpin(seq):
                        return True
        return False

    def evaluate_long(self, seq):
        if self._rule_GC_content(seq):
            if self._rule_slippage(seq):
                if self._rule_runs(seq):
                    if self._rule_palindrome(seq):
                        if self._rule_melting_temp(seq):
                            if self._rule_hairpin(seq):
                                return True
        return False

    def evaluate(self, seq):
        if len(seq) <= 500:
            return self.evaluate_short(seq)
        else:
            return self.evaluate_long(seq)



if __name__ == '__main__':
    import doctest
    doctest.testmod(extraglobs={'synth_obj': Synthesis()})

    seq = 'GGGCAGATGTCACTCAGTAGGATTTGATCAGTACGTGCATGATGCTGCGGGCAGATGTCACTCAGTAGGATTTGATCAGTACGTGCATGATGCTGC'*100#*1000
    # seq = 'AGTCGCAGACAAAAAAAAAAATTTTTTTTTTTCGATGCTAGTCG'

    # seq = 'GCGC GACG ATGA CCGG CAGA CCGG TCAT CGTC ATAT'
    # seq = 'GCGCGACGATGACCGGCAGACCGGTCATCGTCATAT'

    # synth_obj = synthesis()
    # t0 = time()
    # print synth_obj.evaluate(seq)
    # print 'Time elapsed: {}'.format(time()-t0)


    # t0 = time()
    # # print synth_obj._rule_GC_content(seq)
    # # print synth_obj._rule_melting_temp(seq)
    # # print synth_obj._rule_slippage(seq)
    # # print synth_obj._rule_runs(seq)
    # # print synth_obj._rule_palindrome(seq)
    # print synth_obj._rule_hairpin(seq)
    # print 'Time elapsed: {}'.format(time()-t0)