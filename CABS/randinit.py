from CABS.atom import Atoms


class RandomInitialStructure(Atoms):

    # this class tries to build initial random structure from the source variable
    # source can be SEQUENCE, SEQUENCE:SECONDARY or filename

    def __init__(self, source):

        try:
            # try reading from file
            with open(source, 'rb') as f:
                source = f.read()

        except IOError:
            # source is considered string from now on
            pass

        # remove white spaces and make capital
        source = source.replace(' ', '').replace('\n', '').upper()

        super(RandomInitialStructure, self).__init__(source)
        self.change_chid('X', 'A')

    @property
    def pdb(self):
        return self.random_conformation()


if __name__ == '__main__':
    import sys
    ris = RandomInitialStructure(sys.argv[1])
    print(ris.pdb)
