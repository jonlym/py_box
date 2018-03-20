try:
    import numpy as np
except:
    pass
    
class Site(object):
    def __init__(self, x = None, y = None, z = None, position = None, neighbors = None, info = None):
        if position is not None:
            self.set_position(position)
        else:
            self.x = x
            self.y = y
            self.z = z
        self.neighbors = neighbors
        self.info = info

    def set_position(self, position):
        self.x = position[0]
        self.y = position[1]
        self.z = position[2]

    def get_position(self):
        return np.array([self.x, self.y, self.z])
