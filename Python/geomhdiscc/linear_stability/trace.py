"""Module provides the functions to compute linear stability curves"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import geomhdiscc.linear_stability.io as io
import os.path

def marginal(pbm, res, ks, params, save_file = True, new_file = False, db_file = False):
    """Compute the marginal curve for a given problem"""

    # Check for file output
    if save_file:
        # Exact filename is provided
        if isinstance(save_file, str):
            filename = save_file
            save_file = True
        # Use default naming scheme
        else:
            filename = 'marginal_' + pbm + '.dat'

        # Open file
        if new_file:
            f = open(filename, 'w')
        else:
            f = open(filename, 'a')
    else:
        points = []

    # Write file header if required
    if save_file and (new_file or os.path.isfile(filename)):
        io.write_header(f, pbm, res, params)

    # Compute a unique index
    tot = 1
    key_size = {}
    for k,v in params.items():
        try:
            sze = len(v)
        except:
            sze = 1
        tot = sze*tot
        key_size[k] = sze

    # Loop over all parameters
    for i in range(tot):
        # Create parameters for current run
        idx = dict(zip(list(key_size.keys()),np.unravel_index(i, tuple(key_size.values()))))
        run = {}
        for k,v in params.items():
            try:
                val = v[idx[k]]
            except:
                val = v
            run[k] = val

        # Loop over the wave numbers
        for kc in ks:
            # Compute critical values
            racs, omegas = find_critical(pbm, res, kc, rac0, run)

            # Print run information
            if False:
                print('\t--> Converged to critical value: ')
                print('\t--> Converged to critical frequency: ')
            else:
                print('\t--> Failed to converge to critical value!')


            # Save results to file is required
            if save_file:
                io.write_results(f, res, kc, racs, omegas, run)
            else:
                if True:
                    points.append([kc]+ racs + omegas)


    # Close file if required
    if save_file:
        f.close()
    else:
        return points


def find_critical(pbm, res, kc, rac0, params):
    """Find the critical value for given parameters"""

    # Loop over the requested modes
    for n in range(1,params['mode']+1):
        # Compute approximation with Brent's root finding algorithm
        print('Start root finding algorithm')
