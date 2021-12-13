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
    egadsAim.input.Edge_Point_Max = 20
    
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
    
    tacsAim.input.Material = {"Madeupium": madeupium}
    
    # Material properties section
    shell  = {"propertyType" : "Shell",
              "membraneThickness" : 0.006,
              "material"        : "madeupium",
              "bendingInertiaRatio" : 1.0, # Default
              "shearMembraneRatio"  : 5.0/6.0} # Default
    
    tacsAim.input.Property = {"plate": shell}
    
    # constraint section
    constraint = {"groupName" : "edge",
                  "dofConstraint" : 123456}
    
    tacsAim.input.Constraint = {"edgeConstraint": constraint}
    
    #loads section
    pressload = 1.0e5
    
    # Set load
    load = {"groupName" : "plate",
            "loadType" : "Pressure",
            "pressureForce" : pressload}
    
    # Set loads
    tacsAim.input.Load = {"appliedPressure": load }

#ParOpt Optimization Class
class Optimization(ParOpt.Problem):
    def __init__(self):
        self.problem = Caps2Tacs("panel.csm", setupCAPS)
        
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
    def printObjCon(self, x):
        fail, obj, con = self.evalObjCon(x)
        print("\n")
        print("Objective Value: ",obj)
        print("Constraint Value(s): ",con)

## Optimization problem defined here ##
        
myOpt = Optimization()

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