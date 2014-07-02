"""Module provides the functions to compute linear stability curves"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np


def write_header(f, pbm, res, params):
    """Write marginal curve file header"""

    # First header
    f.write('#Results for: ' + pbm + '\n')

    header = []

    # Resolution
    for i,r in enumerate(res):
        header.append('Res_' + str(i))

    # Parameters
    for k,v in sorted(params.items()):
        header.append(k)

    # Wave number
    header.append('k')

    # Critical values
    for i in range(1, params['mode']+1):
        header.append('Rac_' + str(i))

    # Critical frequencies
    for i in range(1, params['mode']+1):
        header.append('Omega_' + str(i))

    f.write('#'+'\t'.join(header) + '\n')


def write_results(f, res, kc, racs, omegas, params):
    """Write marginal curve point to file"""

    result = []

    # Resolution
    for r in res:
        result.append(str(r))

    # Parameters
    for k,v in sorted(params.items()):
        result.append(str(v))

    # Wave number
    result.append(str(kc))

    # Critical values
    for rac in racs:
        result.append(str(rac))

    # Critical frequencies
    for omega in omegas:
        result.append(str(omega))

    f.write('\t'.join(result) + '\n')
