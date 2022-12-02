# Implements a (rather crude) GUI for adsorbate placer.
# To do - comment and perhaps re-code adsorbate placer to be OOP
# for easier re-use of variables.

from tkinter import *

def rotationbox_gui(RotationObj, output="output.xyz"):

    import copy

    global x, y, z
    global _RotationObj, _output_name

    _RotationObj = copy.deepcopy(RotationObj)

    master = Tk()

    x = Scale(master, from_=0, to=360, label='x')
    y = Scale(master, from_=0, to=360, label='y')
    z = Scale(master, from_=0, to=360, label='z')

    x.grid(row=0,column=0)
    y.grid(row=0,column=1)
    z.grid(row=0,column=2)

    _output_name=output

    button = Button(master, text="Print Position", command=print_position)
    button.grid(row=1, column=1)
    button2 = Button(master, text="Start ASE", command=start_ase_gui)
    button2.grid(row=1, column=2)

    mainloop()

def start_ase_gui():

    from ase.gui.gui import GUI

    gui = GUI()

    x_rot, y_rot, z_rot = get_rotations()

    _RotationObj.rotate([x_rot, y_rot, z_rot])

    gui.new_atoms(_RotationObj.ads_and_site)
    gui.repeat_poll(main_ase_gui, 1000)

    gui.run()

def main_ase_gui(gui):

    x_rot, y_rot, z_rot = get_rotations()

    _RotationObj.rotate([x_rot, y_rot, z_rot])

    gui.atoms.positions = _RotationObj.ads_and_site.positions
    gui.draw()

def print_position():

    from ase.io import write

    print(f"Writing Atoms Object to... {_output_name}")
    write(_output_name,_RotationObj.ads_and_site)

def get_rotations():

    x_rot = x.get()
    y_rot = y.get()
    z_rot = z.get()

    return x_rot, y_rot, z_rot