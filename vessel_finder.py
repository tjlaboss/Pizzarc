# Vessel finder

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from pressure_vessel import PressureVessel

Z = 1
P_OUT = 10  # MPa
P_IN0_MIN = 0.69  # Standard BWR fill gas
P_IN0_MAX = 10
P_IN1 = P_OUT

FUEL_EXPANSION = 1.3  # 30% growth by volume
DR = FUEL_EXPANSION**(1/3)


def vgap(Rgap, Rfuel):
	return Z*np.pi*(Rgap**2 - Rfuel**2)


def find_rgap(rf0, p_in0):
	def solve_it(rg):
		vg0 = vgap(rg, rf0)
		vg1 = vgap(rg, DR*rf0)
		return vg0/vg1 - P_IN1/p_in0
	
	rgap = fsolve(solve_it, 1.1*rf0)
	return rgap[0]


if __name__ == "__main__":
	
	pressures = np.linspace(P_IN0_MIN, 0.8*P_IN0_MAX, num=6)
	rfuels = np.arange(1.0, 2.1, step=0.25)
	thickness_ratios = np.linspace(0.02, 0.10, num=9)
	nt = len(thickness_ratios)
	#thickness_ratios = [0.1]
	result_array = np.zeros((len(pressures), len(rfuels), nt))
	for i, p_in in enumerate(pressures):
		print("Initial pressure: {:.1f}".format(p_in))
		for j, rf in enumerate(rfuels):
			rg = find_rgap(rf, p_in)
			print("\tRfuel: {:.2f},\tRgap: {:.2f}".format(rf, rg))
			thicknesses = [foo*rg for foo in thickness_ratios]
			for k, thicc in enumerate(thicknesses):
				pv0 = PressureVessel(rg, thicc, p_in, P_OUT)
				pv0max = pv0.get_max_stress()
				pv1 = PressureVessel(rg, thicc, P_OUT, P_OUT)
				pv1max = pv1.get_max_stress()
				print("\t\tt_clad: {:.2f};\tstress: {:.1f} -> {:.1f}".format(
					thicc, pv0max, pv1max))
				result_array[i, j, k] = max(pv0max, pv1max)
	
	fig = plt.figure()
	
	cmin, cmax = result_array.min(), result_array.max()
	for k, thicc in enumerate(thickness_ratios):
		ax = plt.subplot(3, 3, k+1)
		these = result_array[...,k]
		c = plt.imshow(these, cmap='nipy_spectral')
		plt.ylabel("Fill gass pressure (MPa)")
		plt.yticks(range(len(pressures)), pressures.round(2))
		plt.xlabel("$R_{fuel}$ (cm)")
		plt.xticks(range(len(rfuels)), rfuels)
		plt.title("t = {:.0%}".format(thicc))
		plt.clim(cmin, cmax)
		plt.colorbar(c)
	
	plt.show()

