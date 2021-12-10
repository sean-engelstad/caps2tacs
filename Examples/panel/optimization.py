# -*- coding: utf-8 -*-
from caps2tacs import Caps2Tacs
from setupFunction import setupCAPS
from paropt import ParOpt
from mpi4py import MPI

#ParOpt Optimization Class
class Optimization(ParOpt.Problem):
    def __init__(self, capsProblem):
        self.problem = capsProblem
        
        nvar = 2 #number of design var
        ncon = 2
        nblock = 1
        super(Optimization, self).__init__(MPI.COMM_SELF, nvar, ncon, nblock)
        self.setBounds()
    def setBounds(self,maxStress=1.0,minStress=0.0):
        self.maxStress = maxStress
        self.minStress = minStress
    def getVarsAndBounds(self, x, lb, ub):
        """Get the variable values and bounds"""
        lb[:] = 1e-3
        ub[:] = 3.0
        x[:] = 0.95
        return
    def evalObjCon(self, x):
        """
        Return the objective, constraint and fail flag
        """
        #run the solver
        self.problem.fullSolve(x[:])

        fail = 0
        obj = self.problem.mass() #mass

        maxConstr = self.maxStress - self.problem.vm_stress()
        minConstr = self.problem.vm_stress() - self.minStress
        con = [maxConstr,minConstr] #vmstress

        return fail, obj, con

    def evalObjConGradient(self, x, g, A):
        """
        Return the objective, constraint and fail flag
        """
        #run the solver
        self.problem.fullSolve(x[:])
        
        fail = 0
        g[:] = self.problem.mass_grad()
        
        stress_grad = self.problem.vm_stress_grad()
        A[0][:] = -stress_grad
        A[1][:] = stress_grad
        
        return fail
        

runOpt = False

mySolver = Caps2Tacs("panel.csm", setupCAPS)

if (not(runOpt)):

    mySolver.printValues()
    #mySolver.checkGradients()
    
else:

    myOpt = Optimization(mySolver)

    myOpt.setBounds(maxStress=100.0)
    
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
    opt = ParOpt.Optimizer(myOpt, options)
    
    #Set a new starting point
    opt.optimize()
    x, z, zw, zl, zu = opt.getOptimizedPoint()
    fail, obj, con = myOpt.evalObjCon(x[:])
    print("Optimized x: ", x[:])
    print("obj, con(x): ",obj,con)