import numpy as np
from copy import copy

class Reactions(object):
    def __init__(self, reactions = None, info = None):
        self._reactions = list(reactions)
        self.info = info

    def append(self, reaction):
        self._reactions.append(reaction)

    def extend(self, reactions):
        self._reactions.extend(reactions)

    def index(self, name):
        for i, reaction in enumerate(self._reactions):
            if reaction.name == name:
                return i
        else:
            return None

    def remove(self, name = None, index = None):
        if name is not None and index is not None:
            raise Exception('Both name and index cannot be specified.')
        elif name is None and index is None:
            raise Exception('Either name or index must be specified.')
        elif name is not None:
            index = self.index(name)
        del self._reactions[index]