# coding: utf-8
# Distributed under the terms of the MIT License.

import random
import collections
import inspect

class RandomPopSet(collections.abc.MutableSet):

    def __init__(self, iterable=None):
        self.map = {}
        self.list = []
        if iterable is not None:
            self |= iterable

    def __len__(self):
        return len(self.map)

    def __iter__(self):
        return self.list.__iter__()

    def __contains__(self, item):
        return item in self.map

    def add(self, item):
        if item not in self.map:
            self.map[item] = len(self.list)
            self.list.append(item)

    def push(self, item):
        return self.add(item)

    def discard(self, item):
        if item in self.map:
            index = self.map[item]
            self.list[index], self.list[-1] = self.list[-1], self.list[index]
            self.map[self.list[index]] = index
            del self.list[-1]
            del self.map[item]

    def get_random(self):
        if self.list:
            random_ind = random.randint(0, len(self.list)-1)
            return self.list[random_ind]

    def random_pop(self):
        item = self.get_random()
        self.discard(item)

        return item

def raise_not_defined():
    fileName = inspect.stack()[1][1]
    line = inspect.stack()[1][2]
    method = inspect.stack()[1][3]

    print("*** Method not implemented: %s at line %s of %s" % (method, line, fileName))
    sys.exit(1)
