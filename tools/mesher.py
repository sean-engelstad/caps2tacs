# -*- coding: utf-8 -*-
import os
from caps2tacs import Caps2Tacs
from paropt import ParOpt
from mpi4py import MPI
from tacs.pytacs import pyTACS
from tacs import functions

##put your caps2tacs function here
##or modify this one
def capsFunction(egadsAim,tacsAim):
    #setup function for panel.csm
    
	#Egads Aim section, for mesh
    egadsAim.input.Edge_Point_Min = 3
    egadsAim.input.Edge_Point_Max = 5
    
    egadsAim.input.Mesh_Elements = "Quad"
    
    egadsAim.input.Tess_Params = [.25,.01,15]
    
    #increase the precision in the BDF file
    tacsAim.input.File_Format = "Small"
    tacsAim.input.Mesh_File_Format = "Small"
    
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
    
    propDict = {}
    for i in range(8):
    	propDict["plate" + str(i+1)] = shell
    
    tacsAim.input.Property = propDict
    
    # constraint section
    con1 = {"groupName" : "edge1",
                  "dofConstraint" : 123}
    con2 = {"groupName" : "edge2",
              "dofConstraint" : 12}
    
    tacsAim.input.Constraint = {"con1": con1,
    "con2" : con2}
    
    #loads section
    pressload = 1.0e5
    
    # Set load
    load = {"groupName" : "plate",
            "loadType" : "Pressure",
            "pressureForce" : pressload}
    
    # Set loads
    #tacsAim.input.Load = {"appliedPressure": load }
    
    #design variables
    tacsAim.input.Design_Variable = {"plateLength" : {},"plateWidth" : {},"stiffHeight" : {}}

#write the design variable names to go with your design var inputs
desvarList = ["plateLength", "plateWidth","stiffHeight"]
#write your design variables
D = [2.0,3.0,0.1]

interface = Caps2Tacs("stiffPanel4.csm", capsFunction, pytacsFunction, desvarList)
interface.buildMesh(D)
interface.tacs.geom.view()
