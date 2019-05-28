import os
import numpy as np
import openmc
from openmc import mgxs
from openmc.stats import Box
from materials import all_materials, colormap
from constants import *
from pincell import Pincell, GuideTube
from math import *


class Slice3D:
	def __init__(self, radius, pitch, clad_type,
	             rodded=False, matdict=all_materials):
		assert clad_type in CLADS
		self.radius = radius
		self.pitch = pitch
		self.clad_type = clad_type
		self.matdict = matdict
		self._geometry = None
		self._lattice = None
		self.pincell = Pincell(self.radius, self.clad_type, self.matdict)
		self.gtube = GuideTube(self.pincell, rodded)
	
	def get_lattice(self):
		upin = self.pincell.build()
		ugtb = self.gtube.build()
		hlat = openmc.HexLattice()
		hlat.pitch = [self.pitch]*2
		
		unis = [
			[upin]*78,
			[upin]*72,
			[upin]*66,
			[upin]*60,
			[upin]*54,
			[upin]*48,
			[upin]*42,
			[upin]*36,
			[upin]*30,
			[upin]*24,
			[upin]*18,
			[upin]*12,
			[upin]*6,
			[ugtb]
		]
		hlat.center = (0, -self.pitch*(len(unis) - 1), 0)
		# Center row
		unis[0][39] = ugtb
		unis[3][30] = ugtb
		unis[6][21] = ugtb
		unis[9][12] = ugtb
		unis[12][3] = ugtb
		unis[13][0] = upin
		unis[12][0] = ugtb
		unis[9][0] = ugtb
		unis[6][0] = ugtb
		unis[3][0] = ugtb
		unis[0][0] = ugtb
		# Next row
		unis[0][36] = ugtb
		unis[3][27] = ugtb
		unis[6][18] = ugtb
		unis[9][9] = ugtb
		unis[10][4] = ugtb
		unis[9][3] = ugtb
		unis[6][3] = ugtb
		unis[3][3] = ugtb
		# Another row
		unis[0][33] = ugtb
		unis[3][24] = ugtb
		unis[6][15] = ugtb
		unis[7][10] = ugtb
		unis[6][6] = ugtb
		unis[3][6] = ugtb
		unis[0][6] = ugtb
		# Another row
		unis[0][30] = ugtb
		unis[3][21] = ugtb
		unis[4][16] = ugtb
		unis[4][14] = ugtb
		unis[3][9] = ugtb
		unis[0][9] = ugtb
		# Corner
		unis[0][26] = ugtb
		
		hlat.universes = [unis]
		hlat.outer = upin
		print(hlat.show_indices(hlat.num_rings))
		return hlat
	
	def build(self):
		x0 = openmc.XPlane(x0=0, boundary_type="reflective")
		
		right_inner_wall = openmc.Plane(A=cos(pi/3), B=-cos(pi/6),
		                                D=(RAD_MAJ - STEEL_THICK)*cos(pi/6))
		right_outer_wall = openmc.Plane(A=cos(pi/3), B=-cos(pi/6),
		                                D=RAD_MAJ*cos(pi/6))
		right_refl_edge = openmc.Plane(boundary_type="reflective",
		                               A=cos(pi/6), B=cos(pi/3), D=-GAP/2*cos(pi/6)*0)
		ymin = openmc.YPlane(y0=-RAD_MAJ, boundary_type="vacuum", name="YMIN")
		
		zfuel = openmc.ZPlane(z0=WIDTH-STEEL_THICK)
		zmin = openmc.ZPlane(z0=0, boundary_type="reflective", name="ZMIN")
		zmax = openmc.ZPlane(z0=WIDTH, boundary_type="vacuum", name="ZMAX")
		
		ru = openmc.Universe(name="root universe")
		radialu = openmc.Universe(name="radial universe")
		self._lattice = self.get_lattice()
		dist = RAD_MAJ - STEEL_THICK - (sqrt(3)/2)*self.pitch/1.1
		right_inner_water = openmc.Plane(A=cos(pi/3), B=-cos(pi/6),
		                                 D=dist*cos(pi/6))
		
		# Fueled area
		inner = openmc.Cell()
		inner.region = -right_inner_water
		inner.fill = self._lattice
		radialu.add_cell(inner)
		# Water buffer
		buffer = openmc.Cell()
		buffer.region = +right_inner_water & -right_inner_wall
		buffer.fill = all_materials["Mod"]
		radialu.add_cell(buffer)
		# Reactor pressure vessel
		rpv = openmc.Cell()
		rpv.region = +right_inner_wall & -right_outer_wall
		rpv.fill = all_materials["SS316"]
		radialu.add_cell(rpv)
		outside = openmc.Cell()
		outside.region = +right_outer_wall
		outside.fill = all_materials["Air"]
		radialu.add_cell(outside)
		
		# Make it axially finite
		axialu = openmc.Universe()
		cmain = openmc.Cell()
		cmain.region = +zmin & -zfuel
		cmain.fill = radialu
		axialu.add_cell(cmain)
		cwall = openmc.Cell()
		cwall.region = +zfuel & -zmax
		cwall.fill = self.matdict["SS316"]
		axialu.add_cell(cwall)
		
		# Root Universe
		root_cell = openmc.Cell(name="root cell")
		root_cell.fill = axialu
		root_cell.region = +x0 & +ymin & +zmin & -zmax & -right_refl_edge
		ru.add_cell(root_cell)
		self._geometry = openmc.Geometry()
		self._geometry.root_universe = ru
	
	def make_tallies(self, folder_name):
		tals = openmc.Tallies()
		# 8-group equal-lethargy bin energy tally
		etal = openmc.Tally()
		efilter = openmc.EnergyFilter([0] + list(np.logspace(-3, 7, 9)))
		etal.filters = [efilter]
		etal.scores = ["flux"]
		tals.append(etal)
		# MGXS tallies
		lib = mgxs.Library(self._geometry)
		lib.energy_groups = mgxs.EnergyGroups([0, 20E6])
		lib.mgxs_types = ["fission", "capture"]
		lib.domain_type = "universe"
		lib.domains = [self._geometry.root_universe]
		lib.by_nuclide = True
		lib.build_library()
		lib.dump_to_file("mgxs", folder_name)
		lib.add_to_tallies_file(tals)
		return tals
	
	def make_plots(self):
		p1 = openmc.Plot(name="Plot-XY")
		p1.basis = "xy"
		p1.color_by = "material"
		p1.width = [RAD_MAJ*1.02]*2
		p1.pixels = [1200, 1200]
		p1.origin = [WIDTH/4, -WIDTH/4, 0]
		
		p2 = openmc.Plot(name="Plot-XZ")
		p2.basis = "xz"
		p2.color_by = "material"
		p2.width = [RAD_MAJ*1.02, WIDTH*1.02]
		p2.pixels = [600, 1200]
		p2.origin = [WIDTH/2, -WIDTH/4, WIDTH/2]
		return [p1]
	
	def export_to_xml(self):
		folder_name = "{clad_type}/radius{radius:.2f}_pitch{pitch:.2f}/".format(**vars(self))
		if not os.path.isdir(folder_name):
			os.makedirs(folder_name)
		self._geometry.export_to_xml(folder_name + "geometry.xml")
		mfile = openmc.Materials()
		for mat in all_materials.values():
			mfile.append(mat)
		mfile.export_to_xml(folder_name + "materials.xml")
		pfile = openmc.Plots()
		pfile += self.make_plots()
		# noinspection PyTypeChecker
		for p in pfile:
			if p.color_by == "material":
				p.colors = colormap
		pfile.export_to_xml(folder_name + "plots.xml")
		tfile = self.make_tallies(folder_name)
		tfile.export_to_xml(folder_name + "tallies.xml")
		sfile = openmc.Settings()
		sfile.particles = 10000
		sfile.batches = 50
		sfile.inactive = 20
		sfile.source = openmc.Source(space=Box([-RAD_MIN/2, -RAD_MAJ, -10], [RAD_MIN/2, 0, 10]))
		sfile.export_to_xml(folder_name + "settings.xml")
		print("Exported to:", folder_name)


if __name__ == '__main__':
	bar = Slice3D(2.25, 9.8, "SS316", rodded=False)
	bar.build()
	bar.export_to_xml()
