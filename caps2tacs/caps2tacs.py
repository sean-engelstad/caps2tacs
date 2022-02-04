#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 13:54:34 2022

@author: sengelstad6
"""

import os
import numpy as np
import pyCAPS

class Caps2Tacs:
    def __init__(self, csmFile, capsFunction, pytacsFunction, desvars):

        #order of design variables used for optimization
        self.desvars = desvars
            
        #initialize the tacs and egads aims
        self.initAims(csmFile, capsFunction)
        
        #store the run pytacs method
        self.pytacsFunction = pytacsFunction
        
        #specify the run settings
        self.reorder = True #reorder nodes
        self.printNodes = True #print nodes in the sens file
        
    def initAims(self, csmFile, capsFunction):
        #initialize egads and tacs AIMs and store them in self.tacs and self.egads
        
        #initialize the caps analysis object 
        self.caps = pyCAPS.Problem('myCAPS', capsFile=csmFile, outLevel=0)
        
        #store caps geometry and desvars
        self.geom = self.caps.geometry
        self.nvar = len(self.desvars)
        
        #generate egads aim
        self.egads = self.caps.analysis.create(aim="egadsTessAIM")
        
        #generate tacs aim
        self.tacs = self.caps.analysis.create(aim = "tacsAIM", name = "tacs")
        
        #put design variables into tacsAim
        # I now require you to put in the design variable dict
        # make sure the same order as your x vector
#        desDict = {}
#        for key in self.desvars:
#            desDict[key] = {}
#        self.tacs.input.Design_Variable = desDict
        
        #run the setupFunction to set mesh, geom, load, mat prop settings 
        #in egads and tacs aims
        capsFunction(self.egads, self.tacs)
        
    def solveStructuralProblem(self, D=None):
        #to solve each structural problem, we update our design variables
        #then we run TACS, compute its func and sens
        #then print & combine sensitivities with caps to get df/dD and store them
        

        #update design with new desvar
        if (not(D is None)):
            design = self.makeDesignDict(D)
            self.updateDesign(design)
        
        #run TACS with new design
        self.runTACS()
        
        #print sensitivity files to caps and combine them
        self.combineSensitivities()
        
        #store results - function values and full sensitivity
        self.storeResults()
    
    def buildMesh(self, D=None):
    	#used if you just want to build the .bdf and .dat files
    	#for the mesh for this design
    
    	#update design with new desvar
        if (not(D is None)):
            design = self.makeDesignDict(D)
            self.updateDesign(design)
            
        #run the tacs preanalysis
        self.tacs.preAnalysis()
        
    def makeDesignDict(self, D):
        D = np.array(D)
        designDict = {}
        for i in range(len(D)):
            designDict[self.desvars[i]] = D[i]
        return designDict
    def updateDesign(self, designDict, output=False):
        if (output): print("Design variables: ", designDict)
        for key in designDict:
            self.geom.despmtr[key].value = designDict[key]
    
    def runTACS(self):
        #build the BDF and data file with CAPS preanalysis
        #then run pytacs on the data file, computing func, sens
        
        #run tacs aim preanalysis to generate BDF and DAT file
        self.tacs.preAnalysis()
        
        #read the data file
        datFile = os.path.join(self.tacs.analysisDir, self.tacs.input.Proj_Name + '.dat')
        
        #compute function values and sensitivities for each SP in your Pytacs method
        self.pytacsFunction(self, datFile)
        
    def combineSensitivities(self):
        
        #run bookkeeping data to get nfunc, nnodes, etc.
        self.bookKeeping()
        
        #reorder the nodes from pytacs, thus rearranging Xpts sensitivity
        if (self.reorder): self.reorderNodes()
        
        #print the sensitivities from tacs to CAPS
        self.printSensitivity()
    
    def reorderNodes(self):
        #as the nodes get reordered via MPI when running pytacs
        #we have to reorder them back using the tacsNodeMap
        #this changes the final sensitivities
        
        #loop over each function that was evaluated to change sensitivities for that
        for key in self.funcKeys:
            
            #get sensitivity w.r.t. mesh aka Xpts sens
            dfdX = self.sens[key]['Xpts']
            dfdX = dfdX.reshape((self.nnodes,3))
            
            #thickness desvars, no reordering for this
            #dfdH = self.pytacsSens[key]['struct']
            
            #get the nodemap from pytacs
            bdfNodes = np.linspace(1.0,self.nnodes,self.nnodes)
            tacsNodeMap = self.FEASolver.meshLoader.getLocalNodeIDsFromGlobal(bdfNodes,nastranOrdering=True)
            
            #reorder nodes as dfdX_bdf
            dfdX_bdf = np.zeros((self.nnodes,3))
            for bdfind in range(self.nnodes):
                    #node 1 in bdf, corresponds to 0 bdfind, and plug in bdfind-1 to get tacsInd cause cyclically off by 1
                    tacsInd = tacsNodeMap[bdfind] #tacsNodeMap is cyclically off by 1
                    dfdX_bdf[bdfind,:] = dfdX[tacsInd,:]
            
            #update function sensitivity
            self.sens[key]['Xpts'] = dfdX_bdf
        
    def bookKeeping(self):
        #get bookkeeping info on func such as funcKeys, nfunctions, nnodes
        self.funcKeys = self.func.keys()
        self.nfunc = len(self.funcKeys)
        funcKeysList = list(self.funcKeys)
        
        #get number of nodes
        subSens = self.sens[funcKeysList[0]]['Xpts']
        self.nnodes = int(len(subSens)/3)
        print(self.nnodes)
    
    def printSensitivity(self):
        #print our dfdX sensitivity to CAPS
        #then this will be combined to produce df/dD sensitivity
        #available in dynout variables
        
        #where to print .sens file
        sensFilename = os.path.join(self.tacs.analysisDir, self.tacs.input.Proj_Name+".sens")
        
        #open the file
        with open(sensFilename, "w") as f:
            
            #write (nfunctions) in first line
            f.write("{}\n".format(self.nfunc))
            
            #for each function mass, stress, etc.
            for key in self.funcKeys:
                
                #get the pytacs/tacs sensitivity w.r.t. mesh for that function
                cSens = self.sens[key]['Xpts']
                
                #write the key,value,nnodes of the function
                f.write(key + "\n")
                f.write("{}\n".format(self.func[key]))
                f.write("{}\n".format(self.nnodes))

                #for each node, print nodeind, dfdx, dfdy, dfdz for that mesh element
                for nodeind in range(self.nnodes): # d(Func1)/d(xyz)
                    
                    #if we print nodes first, then print them
                    if (self.printNodes):
                        bdfind = nodeind + 1
                        f.write("{} {} {} {}\n".format(bdfind, cSens[nodeind,0], cSens[nodeind,1], cSens[nodeind,2]))
                    #otherwise just print df/dx, df/dy, df/dz for each element
                    else:
                        f.write("{} {} {}\n".format(cSens[nodeind,0], cSens[nodeind,1], cSens[nodeind,2]))
        
        #run aim postanalysis
        self.tacs.postAnalysis()
        print("finished postanalysis\n")
    
    def storeResults(self):
        #store the function and full df/dD sensitivities from CAPS aim dynout attributes
        
        #initialize gradient variable
        self.grad = {}
        
        #loop over each pytacs function
        for key in self.funcKeys:
            self.func[key] = self.tacs.dynout[key].value
            self.grad[key] = np.zeros((self.nvar))
            
            #loop over each design variable to get the full df/dD gradient
            ind = 0
            for desvar in self.desvars:
                self.grad[key][ind] = self.tacs.dynout[key].deriv(desvar)
                ind += 1
        #print("finished storing results\n")
        #print(self.grad)
    def printDesignVariables(self, desvar):
        ind = 0
        for desvarName in self.desvars:
            print("{}: {}".format(desvarName, desvar[ind]))
            ind += 1
    def checkGradients(self, x, functions, gradients, names,h=1e-4):
    	#computes finite difference check for functions and gradients
    	#that you provide to it
    
        #Example:
        #functions = [self.mass,self.vm_stress]
        #gradients = [self.mass_grad, self.vm_stress_grad]
        #names = ["mass","ks_vm_stress"]
        #x = np.ones(2)
            #x = np.ones((self.nvar))
        method = 1
        
        nfuncs = len(functions)
        
        if (method == 1):
            errors = []
            for i in range(nfuncs):
                myfunc = functions[i]
                mygrad = gradients[i]
                name = names[i]
                
                x = np.array(x)        
                p = np.random.uniform(size=x.shape)
                p = p / np.linalg.norm(p)
                
                func1 = myfunc(x-p*h)
                func2 = myfunc(x+p*h)
                
                fdGrad = (func2-func1)/2/h
                cgrad = mygrad(x)
                print(cgrad)
                #print(mygrad,p)
                directDeriv = np.dot(cgrad, p)
                
                error = (directDeriv - fdGrad)/fdGrad
                errors.append(error)
            for i in range(nfuncs):
                name = names[i]; error = errors[i]
                print(name + ' FD gradient error',error)
        elif (method == 2):
            fdGrads = {}
            cgrads = {}
            avgErrors = {}
            for i in range(nfuncs):
                myfunc = functions[i]
                mygrad = gradients[i]
                name = names[i]
                
                cgrad = mygrad(x)
                nx = len(cgrad)
                fdGrad = np.zeros((nx))
                
                
                for direc in range(nx):
                    ei = np.zeros((nx))
                    ei[direc] = 1
                    
                    fdGrad[direc] = (myfunc(x+ei*h)-myfunc(x-ei*h))/2/h
                
                avgError = np.linalg.norm(fdGrad - cgrad)/np.linalg.norm(fdGrad)
                
                fdGrads[str(i)] = fdGrad
                cgrads[str(i)] = cgrad
                avgErrors[str(i)] = avgError
                
            for i in range(nfuncs):
                print("grad compare for " + names[i])
                print("finite diff grad ",fdGrads[str(i)])
                print("chain rule grad ",cgrads[str(i)])
                print("avg comp rel error ",avgErrors[str(i)])
        elif (method == 3):
            #examine how well gradient predicts change in function
            errors = []
            nh = 1
            hvec = [1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e0]
            hvec = [1e-4]
            hvec = np.array(hvec)
            dfExact = np.zeros((nh))
            dfPred = np.zeros((nh))
            D = np.array(x) #desvars  

            for hi in range(nh):
                for i in range(nfuncs):
                    h = hvec[hi]
                    myfunc = functions[i]
                    mygrad = gradients[i]
                    name = names[i]
                    
                    p = np.random.uniform(size=D.shape)
                    p = p / np.linalg.norm(p)
                    
                    func1 = myfunc(D)
                    cgrad = mygrad(D)
                    func2 = myfunc(D+p*h)
                    
                    dfExact[hi] = func2-func1
                    dfPred[hi] = np.dot(cgrad, p) * h
                

            for fi in range(nfuncs):
                name = names[i]; error = errors[i]
                string1 = 'd' + name + '_exact'
                string2 = 'd' + name + 'pred'
                print(name + ' FD gradient error',error)

                fig, ax = plt.subplot(nfuncs,1,fi)
                ax.loglog(hvec,dfExact,'b-',label=string1)
                ax.loglog(hvec, dfPred,'g-',label=string2)
                ax.xlabel('step size h in deltaD')
                ax.ylabel('change in f(D)')
                ax.legend()
            plt.show()
