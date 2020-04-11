from ase.io.jsonio import read_json


class CLICommand:
    """Describe and/or visualize ASE data.

    Examples:

        $ # Show band path stored in bandpath.json
        $ ase show bandpath.json
        $ # Show band structure stored in bs.json
        $ ase show bs.json
      """

    @staticmethod
    def add_arguments(parser):
        parser.add_argument('file', nargs='+')

    @staticmethod
    def run(args, parser):

        import matplotlib.pyplot as plt
        colors = 'bgrcmyk'

        for i, fname in enumerate(args.file):
            obj = read_json(fname)
            objtype = obj.ase_objtype
            print('{}: {}'.format(fname, objtype))
            print(obj)

            kw = {}
            if obj.ase_objtype == 'bandstructure':
                ax = plt.gca()
                kw.update(show=False,
                          ax=ax,
                          colors=colors[i])
            # plot() should be uniform among plottable objects
            obj.plot(**kw)

        plt.show()
