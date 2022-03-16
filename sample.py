# build crystal
import math
import matplotlib.pyplot as plt
import numpy as np

natom = 100
alatt = 10

class Atom():
    def __init__(self):
        self.x_cor = None
        self.y_cor = None
    
    def pbc_boundary(self):
        if(self.x_cor<0): self.x_cor += alatt
        if(self.x_cor >= alatt): self.x_cor -= alatt
        if(self.y_cor<0): self.y_cor += alatt
        if(self.y_cor >= alatt): self.y_cor -= alatt
    

def build_crystral():
    atoms = []
    for ii in range(0, natom): # id
        this_atom = Atom()
        this_atom.x_cor = ii % alatt
        this_atom.y_cor = int(ii / alatt)
        this_atom.pbc_boundary()
        atoms.append(this_atom)
    return atoms

def calculate_distance(atomii, atomjj):
    dr = math.sqrt((atomii.x_cor - atomjj.x_cor)**2 + (atomii.y_cor - atomjj.y_cor)**2) 
    return dr

atoms = build_crystral()
number_density = natom / (alatt**2)
dr = 0.1 
rmax = 10
layer_number = int(rmax/dr) # total number of layers
atom_flag = atoms[20]

gr = np.zeros(layer_number, dtype='int') # number of atoms per layer

for atom in atoms:
    drij = calculate_distance(atom_flag, atom)
    if (drij < rmax):
        layer = int(drij/dr)
        r = layer * dr # inner radius
        gr[layer] += 1

count = len(gr)
x = []
y = []
for ii in range(0, count):
    r = ii * dr
    x.append(r)
    s = 2 * math.pi * r * dr 
    if(s > 0):
        rho = gr[ii] / s / natom
    else: rho = 0
    y.append(rho / number_density)

fig = plt.figure(figsize=(20,10))
ax0 = fig.add_subplot(1,2,1)
for atom in atoms:
    plt.scatter(atom.x_cor,atom.y_cor, c='blue', marker='o')
plt.scatter(atom_flag.x_cor, atom_flag.y_cor, c='red', marker='o')
for r in range(1,5):
    circle = plt.Circle((atom_flag.x_cor, atom_flag.y_cor), r, color='red', fill=False)
    plt.gcf().gca().add_artist(circle)
plt.grid()

ax0 = fig.add_subplot(1,2,2)
plt.plot(x,y)
plt.show()
