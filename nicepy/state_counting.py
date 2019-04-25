import numpy as _np
from scipy.integrate import quad as _quad


class StateCounting:
    def __init__(self, mass, ve, je, zpe, energy, deg=1):
        self.mass = mass
        self.ve = ve
        self.je = je
        self.zpe = zpe
        self.emax = energy
        self.pmax = None
        self.energies = None
        self.momentum = None
        self.states = None
        self.counts = None
        self.max_count = None
        self.deg = deg
        self.int_value = None
        self.int_error = None

        self._internal_count()

    @staticmethod
    def _internal_energy(ve, je, v, j):
        energy = ve * (v + 0.5) + je * j * (j + 1)
        #         self.wexe*(v+0.5)**2
        return energy

    def _internal_count(self):
        vmax = int(self.emax / self.ve - 0.5)
        states = []
        energies = []
        self.counts = []

        for v in range(vmax + 1):
            j = 0
            while self.zpe < self._internal_energy(self.ve, self.je, v, j) < self.emax:
                states.append((v, j))
                energies.append(self._internal_energy(self.ve, self.je, v, j))
                j += 1

        self.energies, self.states = zip(*sorted(zip(energies, states)))

        count = 0
        for v, j in self.states:
            count += self.deg * 2 * j + 1
            self.counts.append(count)

        self.max_count = count

    def _function(self, p):
        idx = _np.argmin(_np.abs([self.pmax.magnitude - p - i.to(self.pmax.units).magnitude for i in self.momentum]))

        a = 4 * _np.pi * p ** 2
        b = self.counts[idx]

        return a * b

    def external_count(self, atomic_mass):
        mu = self.mass * atomic_mass / (self.mass + atomic_mass)
        self.pmax = _np.sqrt(2 * mu * self.emax)
        self.momentum = [_np.sqrt(2 * mu * e) for e in self.energies]

        self.int_value, self.int_error = _quad(self._function, 0, self.pmax.magnitude)
