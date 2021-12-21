# -*- coding: utf-8 -*-
from caps2tacs import Caps2Tacs
#can't import setup Function from somewhere else, causes capsLock issue
#from setupFunction import setupCAPS
from paropt import ParOpt
from mpi4py import MPI

def setupCAPS(egadsAim,tacsAim):
    #setup function for panel.csm
    
	#Egads Aim section, for mesh
    egadsAim.input.Edge_Point_Min = 15
    egadsAim.input.Edge_Point_Max = 40
    
    egadsAim.input.Mesh_Elements = "Quad"
    
    egadsAim.input.Tess_Params = [.25,.01,15]
    
    #increase the precision in the BDF file
    tacsAim.input.File_Format = "Large"
    tacsAim.input.Mesh_File_Format = "Large"
    
    # Link the mesh
    tacsAim.input["Mesh"].link(egadsAim.output["Surface_Mesh"])
    
    # Set analysis type
    tacsAim.input.Analysis_Type = "Static"
    
    #materials section    
    madeupium    = {"materialType" : "isotropic",
                    "youngModulus" : 72.0E9 ,
                    "poissonRatio": 0.33,
                    "density" : 2.8E3,
                    "tensionAllow" :  20.0e7}
    
    tacsAim.input.Material = {"madeupium": madeupium}
    
    # Material properties section
    ribshell  = {"propertyType" : "Shell",
              "membraneThickness" : 0.02,
              "material"        : "madeupium",
              "bendingInertiaRatio" : 1.0, # Default
              "shearMembraneRatio"  : 5.0/6.0} # Default
    
    sparshell = {"propertyType" : "Shell",
              "membraneThickness" : 0.05,
              "material"        : "madeupium",
              "bendingInertiaRatio" : 1.0, # Default
              "shearMembraneRatio"  : 5.0/6.0} # Default
    
    tacsAim.input.Property = {"wing": ribshell}
    
    # constraint section
    constraints = {"groupName" : "wingRoot",
                  "dofConstraint" : 123456}
    
    tacsAim.input.Constraint = {"fixRoot": constraints}
    
    # Set load
    liftload = {"groupName"         : "bottomWing",
        "loadType"          : "GridForce",
        "forceScaleFactor"  : 1.0e5,
        "directionVector"   : [0.0, 1.0, 0.0]}
    
    # Set loads
    tacsAim.input.Load = {"lift": liftload }

#ParOpt Optimization Class
class Optimization(ParOpt.Problem):
    def __init__(self):
        self.problem = Caps2Tacs("afrl_wing4.csm", setupCAPS)
        
        self.nvar = 6 #number of design var
        ncon = 2
        nblock = 1
        super(Optimization, self).__init__(MPI.COMM_SELF, self.nvar, ncon, nblock)
        self.setBounds()
    def setBounds(self,maxStress=1.0,minStress=0.0):
        self.maxStress = maxStress
        self.minStress = minStress
    def getVarsAndBounds(self, x, lb, ub):
        """Get the variable values and bounds"""
        #area, aspectRatio, taper, twistAngle, leadingEdgeSweep, dihedral
        initialValues = [40.0,6.0,0.5,5.0,30.0,5.0]
        lowerBounds = [20.0, ]
        upperBounds = [100,10.0,10.0,10.0,50.0,20.0]
        for i in range(self.nvar):
            lb[i] = lowerBounds[i]
            ub[i] = upperBounds[i]
            x[i] = initialValues[i]
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
    def printObjCon(self, x):
        fail, obj, con = self.evalObjCon(x)
        print("\n")
        print("Objective Value: ",obj)
        print("Constraint Value(s): ",con)

## Optimization problem defined here ##

#run options: check, run
option = "check"
#option = "run"

myOpt = Optimization()


if (option == "check"):
    #myOpt.problem.checkGradients()
    x = [1.75002251e+01, 8.227+00, 1.00606643e-01, 2.77022526e+00, 2.87394287e+01, 9.47776381e+00]
    x2 = [ 0.10000001,  9.99998846,  1.60898103,  0.10006751, 10.9068694,   0.10015052]
    myOpt.problem.fullSolve(x2)
    myOpt.problem.viewGeometry(x2)
elif (option == "run"):
    myOpt.setBounds(maxStress=2.0,minStress=0.5)
    
    options = {
        'algorithm': 'mma'}
    
    # Set up the optimizer
    opt = ParOpt.Optimizer(myOpt, options)
    
    #Set a new starting point
    opt.optimize()
    x, z, zw, zl, zu = opt.getOptimizedPoint()
    
    #print optimized information
    myOpt.printObjCon(x)
    myOpt.problem.printDesignVariables(x[:])