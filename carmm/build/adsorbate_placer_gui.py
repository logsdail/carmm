# Implements a (rather crude) GUI for adsorbate placer.
# To do - comment and perhaps re-code adsorbate placer to be OOP
# for easier re-use of variables.

from tkinter import *

def rotationbox_gui(RotationObj, output="output.xyz"):
    '''

    A graphical view of the adsorbate-site system,
    allowing you to visually rotate the adsorbate in place
    and produce a structure file with that configuration.
    Begin the visualisation by pressing "Start ASE",
    rotate the adsorbate using the x, y and z sliders,
    and write the structure to the desired filename
    using the "Print Position" button.

    Args:
        RotationObj: RotationBox object
            The desired adsorbate-site system as a RotationBox() object from build/adsorbate_placer.py
        output: string
            The desired filename to print the structure to (include desired file extension)

    Returns:
        On press of "Print Position" button, file containing Atoms object

    '''

    import copy

    global _x, _y, _z
    global _RotationObj, _output_name

    _RotationObj = copy.deepcopy(RotationObj)

    master = Tk()

    _x = Scale(master, from_=0, to=360, label='x')
    _y = Scale(master, from_=0, to=360, label='y')
    _z = Scale(master, from_=0, to=360, label='z')

    _x.grid(row=0,column=0)
    _y.grid(row=0,column=1)
    _z.grid(row=0,column=2)

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

    x_rot = _x.get()
    y_rot = _y.get()
    z_rot = _z.get()

    return x_rot, y_rot, z_rot