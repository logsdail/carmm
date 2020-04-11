"""
Stream input commands to lammps to perform desired simulations
"""
from ase.parallel import paropen
from ase.utils import basestring as asestring
from ase.calculators.lammps.unitconvert import convert

# "End mark" used to indicate that the calculation is done
CALCULATION_END_MARK = "__end_of_ase_invoked_calculation__"


def lammps_create_atoms(fileobj, parameters, atoms, prismobj):
    """Create atoms in lammps with 'create_box' and 'create_atoms'

    :param fileobj: open stream for lammps input
    :param parameters: dict of all lammps parameters
    :type parameters: dict
    :param atoms: Atoms object
    :type atoms: Atoms
    :param prismobj: coordinate transformation between ase and lammps
    :type prismobj: Prism

    """
    if parameters["verbose"]:
        fileobj.write("## Original ase cell\n".encode("utf-8"))
        fileobj.write(
            "".join(
                [
                    "# {0:.16} {1:.16} {2:.16}\n".format(*x)
                    for x in atoms.get_cell()
                ]
            ).encode("utf-8")
        )

    fileobj.write("lattice sc 1.0\n".encode("utf-8"))

    # Get cell parameters and convert from ASE units to LAMMPS units
    xhi, yhi, zhi, xy, xz, yz = convert(prismobj.get_lammps_prism(),
            "distance", "ASE", parameters.units)

    if parameters["always_triclinic"] or prismobj.is_skewed():
        fileobj.write(
            "region asecell prism 0.0 {0} 0.0 {1} 0.0 {2} ".format(
                xhi, yhi, zhi
            ).encode("utf-8")
        )
        fileobj.write(
            "{0} {1} {2} side in units box\n".format(xy, xz, yz).encode(
                "utf-8"
            )
        )
    else:
        fileobj.write(
            "region asecell block 0.0 {0} 0.0 {1} 0.0 {2} "
            "side in units box\n".format(xhi, yhi, zhi).encode("utf-8")
        )

    symbols = atoms.get_chemical_symbols()
    try:
        # By request, specific atom type ordering
        species = parameters["specorder"]
    except AttributeError:
        # By default, atom types in alphabetic order
        species = sorted(set(symbols))

    species_i = {s: i + 1 for i, s in enumerate(species)}

    fileobj.write(
        "create_box {0} asecell\n" "".format(len(species)).encode("utf-8")
    )
    for sym, pos in zip(symbols, atoms.get_positions()):
        # Convert position from ASE units to LAMMPS units
        pos = convert(pos, "distance", "ASE", parameters.units)
        if parameters["verbose"]:
            fileobj.write(
                "# atom pos in ase cell: {0:.16} {1:.16} {2:.16}\n"
                "".format(*tuple(pos)).encode("utf-8")
            )
        fileobj.write(
            "create_atoms {0} single {1} {2} {3} units box\n".format(
                *((species_i[sym],) + tuple(prismobj.vector_to_lammps(pos)))
            ).encode("utf-8")
        )


def write_lammps_in(lammps_in, parameters, atoms, prismobj,
                    lammps_trj=None, lammps_data=None):
    """Write a LAMMPS in_ file with run parameters and settings."""

    def write_model_post_and_masses(fileobj, parameters):
        # write additional lines needed for some LAMMPS potentials
        if 'model_post' in parameters:
            mlines = parameters['model_post']
            for ii in range(0,len(mlines)):
                fileobj.write(mlines[ii].encode('utf-8'))

        if "masses" in parameters:
            for mass in parameters["masses"]:
                # Note that the variable mass is a string containing
                # the type number and value of mass separated by a space
                fileobj.write("mass {0} \n".format(mass).encode("utf-8"))

    if isinstance(lammps_in, asestring):
        fileobj = paropen(lammps_in, "wb")
        close_in_file = True
    else:
        # Expect lammps_in to be a file-like object
        fileobj = lammps_in
        close_in_file = False

    if parameters["verbose"]:
        fileobj.write("# (written by ASE)\n".encode("utf-8"))

    # Write variables
    fileobj.write(
        (
            "clear\n"
            'variable dump_file string "{0}"\n'
            'variable data_file string "{1}"\n'
        )
        .format(lammps_trj, lammps_data)
        .encode("utf-8")
    )

    if "package" in parameters:
        fileobj.write(
            (
                "\n".join(
                    ["package {0}".format(p) for p in parameters["package"]]
                )
                + "\n"
            ).encode("utf-8")
        )

    # setup styles except 'pair_style'
    for style_type in ("atom", "bond", "angle",
                       "dihedral", "improper", "kspace"):
        style = style_type + "_style"
        if style in parameters:
            fileobj.write('{} {} \n'.format(style, parameters[style]).encode("utf-8"))

    # write initialization lines needed for some LAMMPS potentials
    if 'model_init' in parameters:
        mlines = parameters['model_init']
        for ii in range(0,len(mlines)):
            fileobj.write(mlines[ii].encode('utf-8'))

    # write units
    if 'units' in parameters:
       units_line = 'units ' + parameters['units'] + '\n'
       fileobj.write(units_line.encode('utf-8'))
    else:
       fileobj.write('units metal\n'.encode('utf-8'))

    pbc = atoms.get_pbc()
    if "boundary" in parameters:
        fileobj.write(
            "boundary {0} \n".format(parameters["boundary"]).encode("utf-8")
        )
    else:
        fileobj.write(
            "boundary {0} {1} {2} \n".format(
                *tuple("sp"[int(x)] for x in pbc)
            ).encode("utf-8")
        )
    fileobj.write("atom_modify sort 0 0.0 \n".encode("utf-8"))
    for key in ("neighbor", "newton"):
        if key in parameters:
            fileobj.write(
                "{0} {1} \n".format(key, parameters[key]).encode("utf-8")
            )
    fileobj.write("\n".encode("utf-8"))

    # write the simulation box and the atoms
    if not lammps_data:
        lammps_create_atoms(fileobj, parameters, atoms, prismobj)
    # or simply refer to the data-file
    else:
        fileobj.write("read_data {0}\n".format(lammps_data).encode("utf-8"))

    # Write interaction stuff
    fileobj.write("\n### interactions\n".encode("utf-8"))
    if "kim_interactions" in parameters:
        fileobj.write("{}\n".format(parameters["kim_interactions"]).encode("utf-8"))
        write_model_post_and_masses(fileobj, parameters)

    elif ("pair_style" in parameters) and ("pair_coeff" in parameters):
        pair_style = parameters["pair_style"]
        fileobj.write("pair_style {0} \n".format(pair_style).encode("utf-8"))
        for pair_coeff in parameters["pair_coeff"]:
            fileobj.write(
                "pair_coeff {0} \n" "".format(pair_coeff).encode("utf-8")
            )
        write_model_post_and_masses(fileobj, parameters)

    else:
        # simple default parameters
        # that should always make the LAMMPS calculation run
        fileobj.write(
            "pair_style lj/cut 2.5 \n"
            "pair_coeff * * 1 1 \n"
            "mass * 1.0 \n".encode("utf-8")
        )

    if "group" in parameters:
        fileobj.write(
            (
                "\n".join(["group {0}".format(p) for p in parameters["group"]])
                + "\n"
            ).encode("utf-8")
        )

    fileobj.write("\n### run\n" "fix fix_nve all nve\n".encode("utf-8"))

    if "fix" in parameters:
        fileobj.write(
            (
                "\n".join(["fix {0}".format(p) for p in parameters["fix"]])
                + "\n"
            ).encode("utf-8")
        )

    fileobj.write(
        "dump dump_all all custom {1} {0} id type x y z vx vy vz "
        "fx fy fz\n"
        "".format(lammps_trj, parameters["dump_period"]).encode("utf-8")
    )
    fileobj.write(
        "thermo_style custom {0}\n"
        "thermo_modify flush yes format float %23.16g\n"
        "thermo 1\n".format(" ".join(parameters["thermo_args"])).encode(
            "utf-8"
        )
    )

    if "timestep" in parameters:
        fileobj.write(
            "timestep {0}\n".format(parameters["timestep"]).encode("utf-8")
        )

    if "minimize" in parameters:
        fileobj.write(
            "minimize {0}\n".format(parameters["minimize"]).encode("utf-8")
        )
    if "run" in parameters:
        fileobj.write("run {0}\n".format(parameters["run"]).encode("utf-8"))
    if not (("minimize" in parameters) or ("run" in parameters)):
        fileobj.write("run 0\n".encode("utf-8"))

    fileobj.write(
        'print "{0}" \n'.format(CALCULATION_END_MARK).encode("utf-8")
    )
    # Force LAMMPS to flush log
    fileobj.write("log /dev/stdout\n".encode("utf-8"))

    fileobj.flush()
    if close_in_file:
        fileobj.close()
