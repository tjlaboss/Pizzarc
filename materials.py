import openmc


all_materials = {}
colormap = {}
# Moderator at 10 MPa, 600 K
water = openmc.Material(name="Water @ 10 MPa")
water.add_element("H", 2/3)
water.add_nuclide("O16", 1/3)
water.set_density(density=0.49, units="g/cc")
all_materials["Mod"] = water
# Uranium metal (10% Moly by weight)
m2 = openmc.Material(name="Uranium Moly10")
m2.add_element("U", 0.9, enrichment=19.99, percent_type='wo')
m2.add_element("Mo", 0.1, percent_type='wo')
m2.set_density("g/cc", 19.0)
all_materials["Fuel"] = m2
# Zirc 4 copied from BEAVRS
zr4 = openmc.Material(name="Zircaloy-4")
zr4.set_density('g/cc', 6.55)
zr4.add_nuclide('O16', 0.00125, 'wo')
zr4.add_element('Cr', 0.0010, 'wo')
zr4.add_element('Fe', 0.0021, 'wo')
zr4.add_element('Zr', 0.98115, 'wo')
zr4.add_element('Sn', 0.0145, 'wo')
colormap[zr4] = "gray" #[0.8]*3
all_materials["Zr4"] = zr4
# SS316 - Atlas Steels (worldstainless.com)
ss = openmc.Material(name="Stainless Steel 316")
ss.set_density("g/cc", 8.0)
ss.add_element("Cr", 17, 'wo')
ss.add_element("Mo", 2, 'wo')
ss.add_element("Ni", 12, 'wo')
ss.add_element("Fe", 69, 'wo')
colormap[ss] = "darkgray" #[0.5]*3
all_materials["SS316"] = ss
# Plain old air
air = openmc.Material(name="air")
air.set_density("g/cc", 0.001)
air.add_nuclide("O16", 0.232, 'wo')
air.add_nuclide("N14", 0.755, 'wo')
air.add_nuclide("Ar40",0.013, 'wo')
colormap[air] = "tan"
all_materials["Air"] = air
# Copper
copper = openmc.Material(name="copper")
copper.add_element("Cu", 1)
copper.set_density("g/cc", 9)
colormap[copper] = "orange"
all_materials["Copper"] = copper

