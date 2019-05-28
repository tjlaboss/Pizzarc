from pressure_vessel import PressureVessel
from constants import *
from math import sqrt
from scipy.optimize import fsolve

ahex = sqrt(3)*3/2*RAD_MAJ**2
r_equiv = sqrt(ahex/PI)


def stress(t):
	rpv = PressureVessel(r_equiv, t, PRESSURE_MOD, 0.1)
	#print(rpv.get_max_stress())
	return rpv.get_max_stress()

boiler_t = fsolve(lambda t: stress(t) - CLAD_YIELDS["SS316"]/2, 20)[0]
print(boiler_t)

