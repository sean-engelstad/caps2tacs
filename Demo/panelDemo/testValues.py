#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 10 09:47:08 2021

@author: sengelstad6
"""

from classDemo import Caps2Tacs

tacsproj = Caps2Tacs("panel.csm")

tacsproj.setBounds(maxStress=1.0)

f = []; m = []; s = []; g = {}
cvec = [1.0,0.1,0.001]

i = 0
for c in cvec:
    x = [c,c]

    f.append([tacsproj.function(x,getVec=False)])
    m.append(tacsproj.mass(x))
    s.append(tacsproj.vm_stress(x))
    g[str(i)] = tacsproj.mass_grad(x)
    i += 1
for i in range(len(f)):
    c = cvec[i]
    print("x: ",[c,c])
    print("Func dict mass,stress: ",f[i])
    print("mass(x): ",m[i])
    print("stress(x): ",s[i])
    print("mass grad(x): ",g[str(i)])