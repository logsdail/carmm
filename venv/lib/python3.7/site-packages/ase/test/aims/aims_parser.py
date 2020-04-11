from numpy.linalg import norm
from ase.io import read


def run():
    atoms = read("aims.out", format="aims-output")

    # find total energy in aims.out
    key = "| Total energy corrected        :"
    with open("aims.out") as f:
        line = next(l for l in f if key in l)
        ref_energy = float(line.split()[5])

    assert norm(atoms.get_total_energy() - ref_energy) < 1e-12

    # find force in aims.out
    key = "Total atomic forces (unitary forces cleaned) [eV/Ang]:"
    with open("aims.out") as f:
        next(l for l in f if key in l)
        line = next(f)
        ref_force = [float(l) for l in line.split()[2:5]]

    assert norm(atoms.get_forces()[0] - ref_force) < 1e-12

    # find stress in aims.out
    key = "Analytical stress tensor - Symmetrized"
    with open("aims.out") as f:
        next(l for l in f if key in l)
        # scroll to significant lines
        for _ in range(4):
            next(f)
        line = next(f)
        ref_stress = [float(l) for l in line.split()[2:5]]

    assert norm(atoms.get_stress(voigt=False)[0] - ref_stress) < 1e-12

    # find atomic stress in aims.out
    key = "Per atom stress (eV) used for heat flux calculation"
    with open("aims.out") as f:
        next(l for l in f if key in l)
        # scroll to boundary
        next(l for l in f if "-------------" in l)

        line = next(f)
        xx, yy, zz, xy, xz, yz = [float(l) for l in line.split()[2:8]]
        ref_stresses = [[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]]

    assert norm(atoms.get_stresses()[0] - ref_stresses) < 1e-12


def write_output():
    output = "  Basic array size parameters:\n  | Number of species                 :        1\n  | Number of atoms                   :        8\n  | Number of lattice vectors         :        3\n  | Max. basis fn. angular momentum   :        2\n  | Max. atomic/ionic basis occupied n:        3\n  | Max. number of basis fn. types    :        3\n  | Max. radial fns per species/type  :        5\n  | Max. logarithmic grid size        :     1346\n  | Max. radial integration grid size :       42\n  | Max. angular integration grid size:      302\n  | Max. angular grid division number :        8\n  | Radial grid for Hartree potential :     1346\n  | Number of spin channels           :        1\n\n\n  Input geometry:\n  | Unit cell:\n  |        5.42606753        0.00000000        0.00000000\n  |        0.00000000        5.42606753        0.00000000\n  |        0.00000000        0.00000000        5.42606753\n  | Atomic structure:\n  |       Atom                x [A]            y [A]            z [A]\n  |    1: Species Si            0.03431851       -0.09796859        0.09930953\n  |    2: Species Si            5.44231311        2.73920529        2.78205416\n  |    3: Species Si            2.75321969        0.10000784        2.72715717\n  |    4: Species Si            2.73199531        2.68826367       -0.08575931\n  |    5: Species Si            1.34757448        1.42946424        1.25761431\n  |    6: Species Si            1.35486030        4.13154987        4.06589071\n  |    7: Species Si            4.04177845        1.27675199        4.00805480\n  |    8: Species Si            3.99821025        4.01092826        1.42388121\n\n  +-------------------------------------------------------------------+\n  |              Analytical stress tensor - Symmetrized               |\n  |                  Cartesian components [eV/A**3]                   |\n  +-------------------------------------------------------------------+\n  |                x                y                z                |\n  |                                                                   |\n  |  x        -0.01478211      -0.01327277      -0.00355870           |\n  |  y        -0.01327277      -0.01512112      -0.01367280           |\n  |  z        -0.00355870      -0.01367280      -0.01534158           |\n  |                                                                   |\n  |  Pressure:       0.01508160   [eV/A**3]                           |\n  |                                                                   |\n  +-------------------------------------------------------------------+\n\n  ESTIMATED overall HOMO-LUMO gap:      0.21466369 eV between HOMO at k-point 1 and LUMO at k-point 1\n\n  Energy and forces in a compact form:\n  | Total energy uncorrected      :         -0.630943948216411E+05 eV\n  | Total energy corrected        :         -0.630943948568205E+05 eV  <-- do not rely on this value for anything but (periodic) metals\n  | Electronic free energy        :         -0.630943948919999E+05 eV\n  Total atomic forces (unitary forces cleaned) [eV/Ang]:\n  |   1         -0.104637839735875E+01          0.500412824184706E+00         -0.439789552504239E+00\n  |   2         -0.155820611394662E+00         -0.476557335046913E+00         -0.655396151432312E+00\n  |   3         -0.193381405004926E+01         -0.122454085397628E+01         -0.169259060410046E+01\n  |   4          0.404969041951871E-01          0.457139849737633E+00         -0.128445757910440E+00\n  |   5          0.109984435024380E-01         -0.165609149153507E+00          0.114351292468512E+01\n  |   6          0.663029766776301E+00         -0.814079627100908E-01          0.384378715376525E-04\n  |   7          0.213211510059627E+01          0.918575437083381E+00          0.189666102862743E+01\n  |   8          0.289372843732474E+00          0.719871898810707E-01         -0.123990325236629E+00\n\n\n    - Per atom stress (eV) used for heat flux calculation:\n        Atom   | Stress components (1,1), (2,2), (3,3), (1,2), (1,3), (2,3)\n      -------------------------------------------------------------------\n             1 |     0.9843662637E-01   -0.1027274769E+00    0.7237959330E-01   -0.3532042840E+00    0.2563317062E+00   -0.3642257991E+00\n             2 |     0.1244911861E+00   -0.4107147872E-01   -0.1084329966E+00    0.1201650287E+00   -0.1716383020E+00   -0.4669712541E-01\n             3 |    -0.1019986539E+01   -0.7054557814E+00   -0.8410240482E+00   -0.3714228752E+00   -0.4921256188E+00   -0.7970402772E+00\n             4 |    -0.5372048581E+00   -0.2498902919E+00   -0.2260340202E+00   -0.4368600591E+00    0.8622059429E-01    0.9182206824E-01\n             5 |    -0.3268304136E-01   -0.1853638313E+00    0.8046857169E-01   -0.3825550863E+00    0.3088175411E+00   -0.2399437437E+00\n             6 |    -0.2682129292E+00   -0.3832959470E+00   -0.5895171406E+00   -0.8151368635E-02    0.5046578049E-01   -0.6756388823E+00\n             7 |    -0.6970248515E+00   -0.6819450154E+00   -0.9123466446E+00   -0.5254451278E+00   -0.5070403877E+00   -0.6281674944E+00\n             8 |    -0.2933806554E-01   -0.6593089867E-01    0.7360641037E-01   -0.1629233327E+00   -0.9955320981E-01    0.4755870988E+00\n      -------------------------------------------------------------------\n\n\n          Have a nice day.\n------------------------------------------------------------\n"

    with open("aims.out", "w") as f:
        f.write(output)


write_output()
run()
