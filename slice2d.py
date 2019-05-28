import os
import openmc
from openmc.stats import Box
from materials import all_materials, colormap
from constants import *
from pincell import Pincell, GuideTube
from math import *


class Slice2D:
	def __init__(self, radius, pitch, clad_type,
	             rodded=False, matdict=all_materials):
		assert clad_type in CLADS
		self.radius = radius
		self.pitch = pitch
		self.clad_type = clad_type
		self.matdict = matdict
		self._geometry = None
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
		unis[4][27] = ugtb
		unis[8][15] = ugtb
		unis[9][0] = ugtb
		unis[4][0] = ugtb
		unis[0][0] = ugtb
		# Next row
		unis[0][35] = ugtb
		unis[4][23] = ugtb
		unis[9][4] = ugtb
		unis[8][11] = ugtb
		unis[4][4] = ugtb
		# Another row
		unis[0][31] = ugtb
		unis[4][19] = ugtb
		unis[5][13] = ugtb
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
		zmin = openmc.ZPlane(z0=-10, boundary_type="periodic", name="ZMIN")
		zmax = openmc.ZPlane(z0=+10, boundary_type="periodic", name="ZMAX")
		
		ru = openmc.Universe(name="root universe")
		radialu = openmc.Universe(name="radial universe")
		lattice = self.get_lattice()
		dist = RAD_MAJ - STEEL_THICK - sqrt(3)/2*self.pitch
		right_inner_water = openmc.Plane(A=cos(pi/3), B=-cos(pi/6),
		                                D=dist*cos(pi/6))
		# Fueled area
		inner = openmc.Cell()
		inner.region = -right_inner_water
		inner.fill = lattice
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
		# Root Universe
		root_cell = openmc.Cell(name="root cell")
		root_cell.fill = radialu
		root_cell.region = +x0 & +ymin & +zmin & -zmax & -right_refl_edge
		ru.add_cell(root_cell)
		self._geometry = openmc.Geometry()
		self._geometry.root_universe = ru
	
	def make_plots(self):
		p1 = openmc.Plot()
		p1.color_by = "material"
		p1.width = [RAD_MAJ*1.02]*2
		p1.pixels = [1200, 1200]
		p1.origin = [WIDTH/4, -WIDTH/4, 0]
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
		sfile = openmc.Settings()
		sfile.particles = 10000
		sfile.batches = 50
		sfile.inactive = 20
		sfile.source = openmc.Source(space=Box([-RAD_MIN/2, -RAD_MAJ, -10], [RAD_MIN/2, 0, 10]))
		sfile.export_to_xml(folder_name + "settings.xml")
		print("Exported to:", folder_name)


if __name__ == '__main__':
	bar = Slice2D(2, 8, "Zr4", rodded=True)
	bar.build()
	bar.export_to_xml()
