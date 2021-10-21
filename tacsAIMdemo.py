# -*- coding: utf-8 -*-

from tacs import TACS

"""
TODO:
    - create tacs assembler: tacs
    - define the objective function: funcs[0]
    - evaluate the adjoint: adjoint
"""

#put this in pytacs or get comm and the tacs assembler object from them

# Compute the total derivative w.r.t. nodal locations, from tacs/tacs/crm/crm.py
fXptSens = assembler.createNodeVec()
product = assembler.createNodeVec()
assembler.evalXptSens(funcs[0], fXptSens)
assembler.evalAdjointResXptSensProduct(adjoint, product)
fXptSens.axpy(-1.0, product)

#total derivative of obj/cons func w.r.t. nodal positions 
fXptSens_nparray = fXptsSens.getArray()

# Gather from all mpi procs
# ..
global_array = comm.allgather(fXptsSens_nparray)

# Flatten the global vector
flattened_global_array = np.concatenate(global_array)

