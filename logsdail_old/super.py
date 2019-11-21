import ase.io.vasp
x = 2 
y = 2
z = 2
cell = ase.io.vasp.read_vasp("POSCAR_CUBE")
ase.io.vasp.write_vasp("POSCAR.NEW",cell*(x,y,z),direct=True,sort=False)
