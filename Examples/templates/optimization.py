# -*- coding: utf-8 -*-
from caps2tacs import Caps2Tacs
from paropt import ParOpt

    def initParOpt(self):
        # Set up the topology optimization problem
        self.nvar = len(self.deskeys)
        ncon = 2
        nblock = 1
        #maybe also need numelems
        super(Caps2Tacs, self).__init__(MPI.COMM_SELF, self.nvar, ncon, nblock)

def setScales(self,mass=1.0,stress=1.0,x=1.0):
    #set the scale factors for each function and design
    self.massScale = 1.0
    self.stressScale = 1.0
    self.xScale = 1.0
    
def mass(self):
    return self.functions[1] /self.massScale
def vm_stress(self):
    return self.functions[0] /self.stressScale
def mass_grad(self,scale=1):
    return self.gradients[1,:] /self.massScale * self.xScale
def vm_stress_grad(self, x=None):
    return self.gradients[0,:] /self.stressScale * self.xScale

def checkGradients(self, x=None, functions=None, gradients=None, names = None,h=1e-4):
    if (functions is None):
        functions = [self.mass,self.vm_stress]
    if (gradients is None):
        gradients = [self.mass_grad, self.vm_stress_grad]
    if (names is None):
        names = ["mass","ks_vm_stress"]
    if (x is None):
        x = np.ones(2)
        #x = np.ones((self.nvar))
    
    #run the solution
    self.fullSolve(x)

    #compute the FD errors
    FDerrors = []
    for i in range(2):
        myfunc = functions[i]
        mygrad = gradients[i]
        name = names[i]
        
        x = np.array(x)        
        p = np.random.uniform(size=x.shape)
        p = p / np.linalg.norm(p)
        
        #print(p)
        
        func1 = myfunc(x-p*h)
        func2 = myfunc(x+p*h)
        
        fdGrad = (func2-func1)/2/h
        cgrad = mygrad(x)
        #print(mygrad,p)
        directDeriv = np.matmul(cgrad, p)
        
        FDerror = (directDeriv - fdGrad)/fdGrad
        FDerrors.append(FDerror)
    
    #print the FD errors
    for i in range(2):
        names = name[i]; FDerror = FDerrors[i]
        print(name + ' FD gradient error',FDerror)
        
def setupOptimization(self,maxStress):
    #can later be changed to allow more problems
    
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
    #run full solve for functions and gradients
    self.fullSolve(x[:])
    
    fail = 0
    obj = self.mass() #mass
    con = [1.0 - self.vm_stress()] #vmstress

    return fail, obj, con

def evalObjConGradient(self, x, g, A):
    """
    Return the objective, constraint and fail flag
    """
    #compute gradients
    self.fullSolve(x)
    
    fail = 0
    g[:] = self.mass_grad(x[:])
    A[0][:] = -self.vm_stress_grad(x[:])/self.maxStress

    return fail
