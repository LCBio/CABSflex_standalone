"""Module for handling distance restraints"""

from CABS.utils import pep2pep1
import random


class Restraint:
    """Class represents single distance restraint"""
    def __init__(self, line, is_side_chain=False):
        words = line.split()
        if len(words) == 4:
            i1, i2, dist, weight = words
            wmin = wmax = float(weight)
            width = 0.0 if is_side_chain else 1.0
        elif len(words) == 5:
            i1, i2, dist, width, weight = words
            wmin = wmax = float(weight)
            width = float(width)
        elif len(words) == 6:
            i1, i2, dist, width, wmin, wmax = words
            wmin = float(wmin)
            wmax = float(wmax)
            width = float(width)
        else:
            raise Exception('Invalid restraint format: "%s"' % line)

        self.id1 = pep2pep1(i1)
        self.id2 = pep2pep1(i2)
        self.distance = float(dist)
        self.width = width
        self.weight_min = wmin
        self.weight_max = wmax
        self.sg = is_side_chain

    def __repr__(self):
        s = '%s %s %.4f %.4f' % (self.id1, self.id2, self.distance, self.width)
        if self.weight_min == self.weight_max:
            s += ' %.2f' % self.weight_max
        else:
            s += ' %.2f %.2f' % (self.weight_min, self.weight_max)
        if self.sg:
            s += ' SG'
        return s

    def update_id(self, ids):
        self.id1 = ids[self.id1]
        self.id2 = ids[self.id2]
        return self


class Restraints:
    """Container for Restraint(s)"""
    def __init__(self, restraints, sg=False):
        self.data = []
        if restraints:
            self.data.extend(Restraint(line, sg) for line in restraints)

    @classmethod
    def from_parser(cls, restraints, sg=False):
        as_strings = [' '.join(str(w) for w in r) for r in restraints]
        return cls(as_strings, sg)

    @classmethod
    def from_file(cls, filename, sg=False):
        with open(filename) as f:
            as_strings = [line for line in f]
        return cls(as_strings, sg)

    def __repr__(self):
        return '\n'.join(str(r) for r in self.data)

    def __iadd__(self, other):
        self.data.extend(other.data)
        return self

    def update_id(self, ids):
        for restr in self.data:
            restr.update_id(ids)
        return self

    def reduce_by(self, factor):
        if 0. < factor < 1.:
            _count = int(len(self.data) * (1. - factor))
            self.data = random.sample(self.data, _count)

    def retain_percentage(self, percentage):
        if 0 < percentage < 100:
            _count = int(len(self.data) * (percentage / 100))
            self.data = random.sample(self.data, _count)


if __name__ == '__main__':
    pass
