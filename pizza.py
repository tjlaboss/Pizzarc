from openmc import *
from openmc.stats import Box
from materials import all_materials, colormap
from pincell import Pincell
from constants import *


# PINCELL STUFF
RADIUS = 1.0
PITCH = 3.0
HPITCH = PITCH/2.0


x0 = XPlane(x0=0, boundary_type="reflective")
y0 = YPlane(y0=0)
z0 = ZPlane(z0=0)

right_inner_wall = Plane(A=cos(pi/3), B=-cos(pi/6), D=(RAD_MAJ - STEEL_THICK)*cos(pi/6))
right_outer_wall = Plane(A=cos(pi/3), B=-cos(pi/6), D=RAD_MAJ*cos(pi/6))
right_refl_edge = Plane(boundary_type="reflective",
                        A=cos(pi/6), B=cos(pi/3), D=-GAP/2*cos(pi/6) * 0)
ymin = YPlane(y0=-RAD_MAJ, boundary_type="vacuum", name="YMIN")
zmin = ZPlane(z0=-10, boundary_type="periodic",    name="ZMIN")
zmax = ZPlane(z0=+10, boundary_type="periodic",    name="ZMAX")

ru = Universe(0, name="root universe")
radialu = Universe(1, name="radial universe")

upin = Pincell(RADIUS, "SS316", all_materials).build()

hlat = HexLattice()
hlat.center = [0]*3
hlat.pitch = (PITCH, PITCH)
hlat.universes = [[[upin]]]
hlat.outer = upin

# Right slice
inner = Cell()
inner.region = -right_inner_wall
inner.fill = hlat
radialu.add_cell(inner)
rpv = Cell()
rpv.region = +right_inner_wall & -right_outer_wall
rpv.fill = all_materials["SS316"]
radialu.add_cell(rpv)
outside = Cell()
outside.region = +right_outer_wall
outside.fill = all_materials["Air"]
radialu.add_cell(outside)
# Root Universe
root_cell = Cell(name="root cell")
root_cell.fill = radialu
root_cell.region = +x0 & +ymin & +zmin & -zmax & -right_refl_edge
ru.add_cell(root_cell)


def make_plots():
	p1 = Plot()
	p1.color_by = "material"
	p1.width = [RAD_MAJ*1.02]*2
	p1.pixels = [1200, 1200]
	p1.origin = [WIDTH/4, -WIDTH/4, 0]
	
	p2 = Plot()
	p2.color_by = "material"
	p2.width = p1.width
	p2.pixels = p1.pixels
	p2.origin = [0, 0, -5]
	
	p3 = Plot()
	p3.color_by = "cell"
	p3.width = p1.width
	p3.pixels = p1.pixels
	return [p1, p3]



def export_to_xml():
	gfile = Geometry()
	gfile.root_universe = ru
	gfile.export_to_xml()
	mfile = Materials()
	for mat in all_materials.values():
		mfile.append(mat)
	mfile.export_to_xml()
	pfile = Plots()
	pfile += make_plots()
	# noinspection PyTypeChecker
	for p in pfile:
		if p.color_by == "material":
			p.colors = colormap
	pfile.export_to_xml()
	sfile = Settings()
	sfile.particles = 10000
	sfile.batches = 10
	sfile.inactive = 5
	sfile.source = Source(space=Box([-RAD_MIN/2, -RAD_MAJ, -10], [RAD_MIN/2, 0, 10]))
	sfile.export_to_xml()


if __name__ == "__main__":
	print("Exporting to xml...")
	export_to_xml()
	print("...done.")

