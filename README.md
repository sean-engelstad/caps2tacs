# CAPS2TACS 
is an interface code between ESP/CAPS, a software to build parametric geometries, and TACS, the Toolkit for the Analysis of Composite Structures. The software CAPS2TACS allows optimization to be performed of objective functions in TACS a structure solver with respect to design parameters in ESP/CAPS

## Installation/Setup
The software uses a modified release of ESP/CAPS called ESPbeta which hasn't been updated yet. Modifications include changes in tensAllow, shearAllow, compressAllow in nastranUtils.c, an issue with the design sensitivity product from both caps and tacs, and a modification which fixes capsLock issues. The latest version of ESPbeta has not been updated yet, so the example code will probably not work yet for other users. The latest version of ESPbeta should be released to the following site: https://acdl.mit.edu/ESP/archive/.  The user must also install PyNastran to use the Nastran-based interface from ESP/CAPS to TACS with pip install pyNastran.

The latest version of TACS can be found at https://github.com/smdogroup/tacs. Build instructions for TACS are located there as well. We use a new part of tacs    called pytacs which reads .BDF and .DAT files (Nastran-format) of finite element meshes into tacs. The structural problem is then setup with pytacs.

## How to Use
How to use CAPS2TACS: the main class can be found in $ROOT/caps2tacs/caps2tacs.py. The main class Caps2Tacs computes function values in TACS such as wing mass, ks_vm_failure, etc. in terms of the parametric geometries from ESP. Gradients are also provided for optimization. To setup a parametric geometry in ESP/CAPS, build a .csm file such as panel.csm as found in the panel example $ROOT/Examples/panel/panel.csm. The panel example also includes a simple optimization case demonstrating mass minimization with stress-constraints for a simple panel. 
