
"""

EMOPEC python package.

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

from emopec._emopec import (get_expression,
                            predict_spacing,
                            make_library,
                            SEQS_MAX)
