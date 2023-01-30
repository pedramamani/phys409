import numpy as np
from toolbox.constants import *
from rb_atom_config import *


class RbAtom:
    def __init__(self, Rbx):  # Rbx = 0, 1 for Rb85, Rb87 respectively
        self.hf_count = HF_COUNT_VALS[Rbx]
        self.I = I_VALS[Rbx]
        self.gI = gI_VALS[Rbx]
        self.Ahf = Ahf_VALS[Rbx]

        mS_vals = np.arange(-S, S + 1)
        mI_vals = np.arange(-self.I, self.I + 1)
        self.basis = [(mS, mI) for mS in mS_vals for mI in mI_vals]

    def energy_levels(self, B0):
        HB = muB * B0 * np.array([[self._HB_coeff(mS, mS_, mI, mI_)
                                   for mS, mI in self.basis] for mS_, mI_ in self.basis])
        Hhf = self.Ahf * np.array([[self._Hhf_coeff(mS, mS_, mI, mI_)
                                    for mS, mI in self.basis] for mS_, mI_ in self.basis])
        H = HB + Hhf
        E_levels = np.linalg.eigvals(H)
        E_levels.sort()
        return E_levels

    def rf_frequencies(self, B0, both_levels=False):
        Es = self.energy_levels(B0)
        fs = []
        for i in range(self.hf_count[0], self.hf_count[0] + self.hf_count[1] - 1):
            fs.append((Es[i + 1] - Es[i]) / h)
        if both_levels:
            for i in range(0, self.hf_count[0] - 1):
                fs.append((Es[i + 1] - Es[i]) / h)
        fs.sort()
        return np.array(fs)

    def _HB_coeff(self, mS, mS_, mI, mI_):
        return (gS * mS + self.gI * mI) * kdelta(mS, mS_) * kdelta(mI, mI_)

    def _Hhf_coeff(self, mS, mS_, mI, mI_):
        SpIm = np.sqrt((S - mS) * (S + mS + 1) * (self.I + mI) * (self.I - mI + 1)) \
               * kdelta(mS + 1, mS_) * kdelta(mI - 1, mI_)
        SmIp = np.sqrt((S + mS) * (S - mS + 1) * (self.I - mI) * (self.I + mI + 1)) \
               * kdelta(mS - 1, mS_) * kdelta(mI + 1, mI_)
        SzIz = mS * mI * kdelta(mS, mS_) * kdelta(mI, mI_)
        return (SpIm + SmIp) / 2 + SzIz


def kdelta(a, b):
    if a == b:
        return 1
    return 0
