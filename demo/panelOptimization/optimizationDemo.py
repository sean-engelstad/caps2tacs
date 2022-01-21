#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 10 09:44:02 2021

@author: sengelstad6
"""
from classDemo import Caps2Tacs

from paropt import ParOpt

tacsproj = Caps2Tacs("panel.csm")

tacsproj.setBounds(maxStress=100.0)

tacsproj.checkGradients()

options = {
    'algorithm': 'tr',
    'tr_init_size': 0.05,
    'tr_min_size': 1e-6,
    'tr_max_size': 10.0,
    'tr_eta': 0.25,
    'tr_infeas_tol': 1e-6,
    'tr_l1_tol': 1e-3,
    'tr_linfty_tol': 0.0,
    'tr_adaptive_gamma_update': True,
    'tr_max_iterations': 1000,
    'max_major_iters': 100,
    'penalty_gamma': 1e3,
    'qn_subspace_size': 10,
    'qn_type': 'bfgs',
    'abs_res_tol': 1e-8,
    'starting_point_strategy': 'affine_step',
    'barrier_strategy': 'mehrotra_predictor_corrector',
    'use_line_search': False}

options = {
    'algorithm': 'mma'}

# Set up the optimizer
opt = ParOpt.Optimizer(tacsproj, options)

#Set a new starting point
opt.optimize()
x, z, zw, zl, zu = opt.getOptimizedPoint()
fail, obj, con = tacsproj.evalObjCon(x[:])
print("Optimized x: ",x[:])
print("obj, con(x): ",obj,con)
