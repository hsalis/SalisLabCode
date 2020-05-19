#Commonly used Model Functions for use with the Non-Repetitive Parts Calculator (Maker Mode)
import re

def modelFunction_GC_Content(seq, GC_percentage_minimum, GC_percentage_maximum):
    GC_count = seq.count('G') + seq.count('C')
    return GC_percentage_minimum <= (GC_count * 100.0) / len(seq) <= GC_percentage_maximum

def modelFunction_Exclude_Sequences(seq, compiledRegExp):
    if compiledRegExp.search(seq) is None:
        return True
    else:
        return False

