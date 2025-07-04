"""
Random initial structure generator with type annotations.
"""

from typing import Union, BinaryIO
from CABS.structures.atom import Atoms


class RandomInitialStructure(Atoms):
    """
    Random initial structure generator.
    
    This class tries to build initial random structure from the source variable.
    Source can be SEQUENCE, SEQUENCE:SECONDARY or filename.
    """

    def __init__(self, source: Union[str, bytes]) -> None:
        """
        Initialize random initial structure.
        
        Arguments:
            source: sequence string, sequence with secondary structure, or filename
        """
        try:
            # try reading from file
            with open(source, 'rb') as f:
                source = f.read().decode('utf-8')
        except (IOError, TypeError, UnicodeDecodeError):
            # source is considered string from now on
            if isinstance(source, bytes):
                source = source.decode('utf-8')

        # remove white spaces and make capital
        source = source.replace(' ', '').replace('\n', '').upper()

        super().__init__(source)
        self.change_chid('X', 'A')

    @property
    def pdb(self) -> str:
        """Generate PDB string for random conformation."""
        return self.random_conformation()


if __name__ == '__main__':
    import sys
    ris = RandomInitialStructure(sys.argv[1])
    print(ris.pdb)
