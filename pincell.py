import openmc
from scipy.optimize import fsolve
from constants import *



def vgap(Rgap, Rfuel):
	return PI*(Rgap**2 - Rfuel**2)

class Pincell:
	def __init__(self, rfuel, clad_type, matdict):
		assert clad_type in CLADS
		self.clad_type = clad_type
		self.fuel_mat = matdict["Fuel"]
		self.gap_mat = None
		self.clad_mat = matdict[clad_type]
		self.rfuel = rfuel
		self.rgap = self._get_rgap()
		self.rclad = CLAD_RATIOS[clad_type]*self.rgap
		self.mod_mat = matdict["Mod"]
	
	def _get_rgap(self):
		def solve_it(rg):
			vg0 = vgap(rg, self.rfuel)
			vg1 = vgap(rg, DR*self.rfuel)
			return vg0/vg1 - PRESSURE_MOD/PRESSURE_GAP
		rgap = fsolve(solve_it, 1.1*self.rfuel)
		return rgap[0]
	
	def build(self):
		# All cylinders
		fcyl = openmc.ZCylinder(R=self.rfuel, name="Fuel Cylinder")
		gcyl = openmc.ZCylinder(R=self.rgap,  name="Gap Cylinder")
		ccyl = openmc.ZCylinder(R=self.rclad, name="Clad Cylinder")
		# Innermost ring: plain fuel
		fring = openmc.Cell()
		fring.region = -fcyl
		fring.fill = self.fuel_mat
		# Next ring: gap
		gring = openmc.Cell()
		gring.region = +fcyl & -gcyl
		gring.fill = self.gap_mat  # probably void
		# Next ring: clad
		cring = openmc.Cell()
		cring.region = +gcyl & -ccyl
		cring.fill = self.clad_mat
		# Outside: moderator
		moderator = openmc.Cell()
		moderator.region = +ccyl
		moderator.fill = self.mod_mat
		# Finalize
		pincell_universe = openmc.Universe(name="Pincell")
		pincell_universe.add_cells((fring, gring, cring, moderator))
		return pincell_universe
	

class GuideTube:
	def __init__(self, pincell):
		self.pincell = pincell
		self.r_inner = pincell.rclad
		self.r_outer = CLAD_RATIOS[pincell.clad_type]*self.r_inner
		
	
	def build(self):
		# All cylinders
		icyl = openmc.ZCylinder(R=self.r_inner, name="Inside Tube")
		ocyl = openmc.ZCylinder(R=self.r_outer, name="Tube Cylinder")
		# Innermost ring: water area
		iring = openmc.Cell()
		iring.region = -icyl
		iring.fill = self.pincell.mod_mat
		# Next ring: Tube
		tring = openmc.Cell()
		tring.region = +icyl & -ocyl
		tring.fill = self.pincell.clad_mat
		# Outside: moderator
		moderator = openmc.Cell()
		moderator.region = +ocyl
		moderator.fill = self.pincell.mod_mat
		# Finalize
		tubular_universe = openmc.Universe(name="Guide Tube")
		tubular_universe.add_cells((iring, tring, moderator))
		return tubular_universe
		
