
"""

Python wrapper for the Vienna RNA Package by Andreas R. Gruber, Ronny Lorenz,
Stephan H. Bernhart, Richard Neubock, and Ivo L. Hofacker (NAR, 2008).

Copyright 2015 Michael Klausen, Mads Bonde, Morten Sommer, all rights reserved.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

"""


FOLDER = None

#
# Folding functions
# -----------------

def mfe(seq):
    """Get the predicted minimum free energy."""
    if FOLDER is None:
        init()

    return FOLDER.mfe(seq)

#
# Utils
# -----

def init():
    """Initialize the module."""
    global FOLDER
    try:
        import RNA
        FOLDER = ViennaRNALib()
        return
    except ImportError:
        pass

    try:
        run_command([ViennaRNA2.RNAFOLD, '--version'])
        FOLDER = ViennaRNA2()
        return
    except Exception:
        pass

    raise Exception('No usable ViennaRNA installation found.')


def run_command(cmd, input_string=''):
    """Run the specified command and return output"""
    import subprocess
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stdin=subprocess.PIPE,
                         stderr=subprocess.PIPE, universal_newlines=True)

    out, stderr = p.communicate(input=input_string)
    if p.returncode:
        raise Exception('Cmd {} failed:\n{}'.format(cmd[0], stderr))
    return out


#
# ViennaRNA wrappers
# ------------------

class ViennaRNA2:

    """
    ViennaRNA version 2.X commandline wrapper.

    """

    RNAFOLD = 'RNAfold'

    def mfe(self, seq, d=2):
        """RNAfold."""
        output = run_command([self.RNAFOLD, '--noPS', '-d' + str(d)], seq)
        output = output.splitlines()
        if len(output) > 1:
            return self.output_to_brackets_and_dG(output[1])[1]
        else:
            return None, None

    def output_to_brackets_and_dG(self, outputline):
        """Single output line to brackets and dG"""
        outputline = outputline.split()
        brackets = outputline[0]
        dG = float(''.join(outputline[1:]).strip('()'))
        return brackets, dG


class ViennaRNALib:

    """
    ViennaRNA RNALib wrapper.

    Only works in Python version 2.X

    """

    def mfe(self, seq, d=2):
        """RNAfold MFE."""
        RNA.cvar.dangles = d
        RNA.cvar.cut_point = -1
        f = RNA.fold(seq)
        RNA.free_arrays()
        return f[1]
