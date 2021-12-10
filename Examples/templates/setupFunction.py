from caps2tacs import Caps2Tacs

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

#start caps to tacs problem
myprob = Caps2Tacs("panel.csm",setupCAPS)

x = [1.0,0.5]
myprob.fullSolve(x)

myprob.printValues()