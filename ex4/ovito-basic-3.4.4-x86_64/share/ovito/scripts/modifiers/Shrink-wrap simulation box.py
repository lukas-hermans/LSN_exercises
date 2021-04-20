#
# Shrink-wrap simulation box:
#
# A user-defined modifier function which resets the simulation box geometry to 
# exactly match the current axis-aligned bounding box of the particle coordinates.
#

import numpy

def modify(frame, data):

    # There's nothing we can do if there are no particles. 
    if data.particles.count == 0: return

    # Compute bounding box of particle coordinates.
    coords_min = numpy.amin(data.particles.positions, axis=0)
    coords_max = numpy.amax(data.particles.positions, axis=0)

    # Adjust simulation cell matrix (three cell vectors and origin).
    data.cell_[:,:3] = numpy.diag(coords_max - coords_min)
    data.cell_[:, 3] = coords_min