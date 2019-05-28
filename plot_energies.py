import os
import openmc
from openmc import mgxs
import matplotlib.pyplot as plt
import numpy as np

os.environ["OPENMC_CROSS_SECTIONS"] = \
	'/home/share/nukedata/mcnp_endfb71_hdf5/cross_sections.xml'
FOLDER = "SS316/radius2.25_pitch9.80/"
STATEPOINT = "statepoint.50.h5"
sp = openmc.StatePoint(FOLDER + STATEPOINT)
lib = mgxs.Library.load_from_file(directory=FOLDER)
lib.load_from_statepoint(sp)

summ = sp.summary
ru = summ.geometry.root_universe
mgdata = lib.get_xsdata(ru, "U235-fission", nuclide="U235")
micro_xs = mgdata.fission[0]/1712E-6  # TODO: Correctly calculate number density


efilt = sp.filters[1]
etal = sp.get_tally(scores=["flux"], filters=[efilt])

talvals = np.array(etal.mean).squeeze()
stdvals = np.array(etal.std_dev).squeeze()

print(talvals)
print(stdvals)
print(efilt.bins)
bins = efilt.bins[:-1]

fig = openmc.plot_xs(openmc.Nuclide("U235"), ['fission'], temperature=600)
#fig = openmc.plot_xs(openmc.Nuclide("U238"), ['capture'], temperature=600)
ax = fig.gca()
ax.step(bins, talvals, label="$\phi(E)$", color="blue")
#ax.step(bins, talvals+stdvals, color="salmon", label="$\pm\sigma$")
#ax.step(bins, talvals-stdvals, color="salmon")
#ax.step([0, 2E7], [micro_xs]*2, color="red")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("$E$")
ax.set_ylabel("$\phi(E)$")

plt.legend(["U235 fission", "$\phi(E)$", "$\sigma_g(E)$"])
plt.title("Light Water")
plt.show()
