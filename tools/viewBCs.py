    #initialize pytacs with that data file
    obj.FEASolver = pyTACS(datFile)
        
    # Set up TACS Assembler
    obj.FEASolver.initialize()
    
    #choose the functions to evaluate
    evalFuncs = ['wing_mass', 'ks_vmfailure','compliance']

    vec = obj.FEASolver.assembler.createVec()
    vec.getArray()[:] = 1.0
    obj.FEASolver.assembler.applyBCs(vec)
    obj.FEASolver.assembler.setVariables(vec)
    obj.FEASolver.outputViewer.writeToFile('bcs.f5')