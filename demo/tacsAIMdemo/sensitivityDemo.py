import os
NumberOfNode = 20
filename = os.path.join('WorkDir', "tacs.sens")
with open(filename, "w") as f:
    f.write("2 {}\n".format(NumberOfNode)) # Two functionals Number Of Nodes
    f.write("Func1\n") # 1st functional
    f.write("42\n")    # Value of Func1
    for i in range(NumberOfNode): # d(Func1)/d(xyz)
        f.write("{} {} {}\n".format(i, 2*i, 3*i))

    f.write("Func2\n") # 2nd functiona;
    f.write("21\n")    # Value of Func2
    for i in range(NumberOfNode): # d(Func2)/d(xyz)
        f.write("{} {} {}\n".format(3*i, 2*i, i))