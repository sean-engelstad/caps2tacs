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
    def __init__(self, csmFile, capsFunction, pytacsFunction, printNodes=True):
            
        #initialize the tacs and egads aims
        self.initAims(csmFile, capsFunction)
        
        #store the run pytacs method
        self.pytacsFunction = pytacsFunction
        
        #specify the run settings
        self.reorder = True #reorder nodes
        self.printNodes = printNodes #print nodes in the sens file
        
    def initAims(self, csmFile, capsFunction):
        #initialize egads and tacs AIMs and store them in self.tacs and self.egads
        
        #initialize the caps analysis object 
        self.caps = pyCAPS.Problem('myCAPS', capsFile=csmFile, outLevel=0)
        
        #store caps geometry and desvars
        self.geom = self.caps.geometry
        self.desvars = self.geom.despmtr.keys()
        self.nvar = len(self.desvars)
        
        #generate egads aim
        self.egads = self.caps.analysis.create(aim="egadsTessAIM")
        
        #generate tacs aim
        self.tacs = self.caps.analysis.create(aim = "tacsAIM", name = "tacs")
        
        #put design variables into tacsAim
        desDict = {}
        for key in self.desvars:
            desDict[key] = {}
        self.tacs.input.Design_Variable = desDict
        
        #run the setupFunction to set mesh, geom, load, mat prop settings 
        #in egads and tacs aims
        capsFunction(self.egads, self.tacs)
        
    def solveStructuralProblem(self, desvar=None):
        #to solve each structural problem, we update our design variables
        #then we run TACS, compute its func and sens
        #then print & combine sensitivities with caps to get df/dD and store them
        
        #update design with new desvar
        if (not(desvar is None)):
            self.updateDesign(desvar)
        
        #run TACS with new design
        self.runTACS()
        
        #print sensitivity files to caps and combine them
        self.combineSensitivities()
        
        #store results - function values and full sensitivity
        self.storeResults()
        
    def updateDesign(self, desvar, output=False):
        desvar = np.array(desvar)
        if (output): print("Design variables: ", desvar)
        ind = 0
        for key in self.geom.despmtr.keys():
            self.geom.despmtr[key].value = desvar[ind]
            ind += 1    
    
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
            
            #write (nfunctions, nnodes) in first lien
            f.write("{} {}\n".format(self.nfunc, self.nnodes))
            
            #for each function mass, stress, etc.
            for key in self.funcKeys:
                
                #get the pytacs/tacs sensitivity w.r.t. mesh for that function
                cSens = self.sens[key]['Xpts']
                
                #write the value of the function
                f.write(key + "\n")
                f.write("{}\n".format(self.func[key]))
                
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
            for deskey in self.tacs.input.Design_Variable.keys():
                self.grad[key][ind] = self.tacs.dynout[key].deriv(deskey)
                ind += 1
        #print("finished storing results\n")
        #print(self.grad)
    def printDesignVariables(self, desvar):
        ind = 0
        for desvarName in self.desvars:
            print("{}: {}".format(desvarName, desvar[ind]))
            ind += 1