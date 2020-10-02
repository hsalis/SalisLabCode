
"""

EMOPEC web interface

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

import re

import flask
import emopec


#Create a flask app
app = flask.Flask(__name__)


@app.route('/')
def index():
    return flask.render_template('emopec.html')


@app.route('/emopec', methods=['POST', 'GET'])
def make_library():
    if flask.request.method == 'GET':
        return flask.render_template('emopec.html')

    data          = flask.request.get_json()
    leader        = data.get('leader')
    cds           = (data.get('cds') or '').upper().replace('U', 'T')
    nlib          = data.get('n', 10)
    libdir        = data.get('libdir', 'both').lower()
    force_spacing = data.get('force_spacing', False)
    constraints   = data.get('constraints')

    if not leader or len(leader) < 15 or not re.match(r'^[ATGCUatgcu]+$', leader):
        err = 'Invalid leader: "{}" (must be DNA, at least 15 nts)'
        response = flask.jsonify(error=err.format(leader))
        response.status_code = 500
        return response
    else:
        leader = leader.upper().replace('U', 'T')

    if constraints and not re.match(r'^[ATWGRKDCMYHSVBNatwgrkdcmyhsvbn]*$', constraints):
        err = 'Invalid leader constraints: "{}" (must be IUPAC ambigous DNA)'
        response = flask.jsonify(error=err.format(constraints))
        response.status_code = 500
        return response

    if libdir not in ('both', 'up', 'down'):
        err = 'Invalid libdir value "{}" (expected "both", "up", or "down").'
        response = flask.jsonify(error=err.format(libdir))
        response.status_code = 500
        return response

    try:
        nlib = int(nlib)
        if not 1 <= nlib <= 50:
            raise ValueError('nlib must be between 1 and 50')
    except ValueError:
        err = 'Invalid nlib value "{}" (expected integer between 1 and 50).'
        response = flask.jsonify(error=err.format(nlib))
        response.status_code = 500
        return response

    if force_spacing:
        try:
            spc = int(force_spacing)
        except ValueError:
            err = 'Invalid value for forced spacing "{}".'
            response = flask.jsonify(error=err.format(force_spacing))
            response.status_code = 500
            return response

        upstream = leader[:-6-spc]
        sd       = leader[-6-spc:-spc]
        spacing  = leader[-spc:]
        expr = emopec.get_expression(sd, len(spacing))
    else:
        upstream, sd, spacing, expr = emopec.predict_spacing(leader)

    #Find max_expression for pct
    expr_max = emopec.get_expression(emopec.SEQS_MAX, len(spacing))
    expr_pct = expr / expr_max

    #Create the library
    current = 0. if libdir == 'both' else None
    target  = 0. if libdir == 'down' else 1.0
    lib = emopec.make_library(upstream, sd, spacing, cds, n=nlib, target=target,
        current=current, dg_seq=constraints)

    mkoligo = lambda _sd: _create_oligo(upstream, _sd, spacing, cds, 90)
    oligos = sorted([(x, mkoligo(_sd), pct, ddG) for _sd, x, ddG, pct in lib], reverse=True)
    #oligos.sort(reverse=True)

    warnings = []
    if oligos and '.' in oligos[0][1]:
        warnings.append('Sequence not long enough to create full oligo '
            '(Missing positions marked with ".").')

    if len(oligos) != nlib:
        warnings.append('Could not make the requested number of oligos due to constraints.')

    #Prepare return payload
    rtn = {
        'upstream': upstream,
        'sd':       sd,
        'spacing':  spacing,
        'cds':      cds,
        'library':  oligos,
        'expression': expr,
        'expression_percent': expr_pct,
        #'error': False,
        'warnings': warnings,
    }
    return flask.jsonify(**rtn)


def _create_oligo(upstream, new_sd, spacing, cds, oligo_len=90):
    """Create an oligo from a sequence.

        >>> _create_oligo('GCTGTATTTTTCCCTATACAAGTCGCTTAAGGCTTGCCAAC', 'AGGAGT',
        ...               'TTGCCGCC', 'ATGAAGTTTATCATTAAATTGTTCCCGGAAATCAC', 70)
        'TTCCCTATACAAGTCGCTTAAGGCTTGCCAACaggagtTTGCCGCCATGAAGTTTATCATTAAATTGTTC'

    Insufficient sequence::

        >>> _create_oligo('ACAAGTCGCTTAAGGCTTGCCAAC', 'AGGAGT', 'TTGCCGCC',
        ...               'ATGAAGTTTATCA', 70)
        '........ACAAGTCGCTTAAGGCTTGCCAACaggagtTTGCCGCCATGAAGTTTATCA...........'

    """
    l_hom = (oligo_len - len(new_sd)) // 2
    r_hom = l_hom + (oligo_len - len(new_sd)) % 2

    l_seq = upstream[-l_hom:]
    if len(l_seq) < l_hom:
        l_seq = '.' * (l_hom - len(l_seq)) + l_seq

    r_seq = (spacing + cds)[:r_hom]
    if len(r_seq) < r_hom:
        r_seq += '.' * (r_hom - len(r_seq))

    return l_seq + new_sd.lower() + r_seq


if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0')
