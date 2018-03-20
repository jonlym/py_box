import numpy as np
from py_box3 import any_alpha

class dos_local_charge:
    """Stores the localized charge found in the OUTCAR file after a density of states run."""
    def __init__(self, ss = [], ps = [], ds = [], totals = [], atoms = None):
        self.ss = ss
        self.ps = ps
        self.ds = ds
        self.totals = totals
        self.atoms = atoms

    @classmethod
    def from_OUTCAR(cls, OUTCAR_path = 'OUTCAR', atoms = None):
        with open(OUTCAR_path, 'r') as outcar:
            record_charge = False
            ss = []
            ps = []
            ds = []
            totals = []
            for i, line in enumerate(outcar):
                if ' total charge     ' in line:
                    record_charge = True
                    continue
                if record_charge:
                    if 'tot  ' in line:
                        break
                    elif (not any_alpha(line) and ('--' not in line) and len(line) > 3):
                        data = [np.float(x) for x in line.split(' ') if x != '']
                        ss.append(data[1])
                        ps.append(data[2])
                        ds.append(data[3])
                        totals.append(data[1] + data[2] + data[3])
        return cls(ss = ss, ps = ps, ds = ds, totals = totals, atoms = atoms)

    def __str__(self):
        lines = []
        lines.append("Element[index]        s        p        d        Total\n")
        lines.append("------------------------------------------------------\n")
        for atom, s, p, d, total in zip(self.atoms, self.ss, self.ps, self.ds, self.totals):
            atom_info = '{}[{}]'.format(atom.symbol, atom.index)
            pad_length = len('Element[index]') - len(atom_info)
            lines.append('{}{}{:9.2f}{:9.2f}{:9.2f}{:13.2f}\n'.format(atom_info, ' '*pad_length, s, p, d, total))
        return ''.join(lines)
