import math

PI = math.pi

PRESSURE_MOD = 10
PRESSURE_GAP =  5

CLADS = ("Zr4", "SS316")
CLAD_RATIOS = {"Zr4"  : 1.07,
               "SS316": 1.05}
CLAD_YIELDS = {"Zr4"  : 160,
               "SS316": 205}

DR = 1.3**(1/3)

# dimensions in cm
WIDTH = 240
GAP = 5
STEEL_THICK = 11
# hexagon stuff
RAD_MAJ = WIDTH/2
RAD_MIN = RAD_MAJ/math.sqrt(3)
WALL2 = WIDTH/8
# Thermal expansion joint
EXPJOINT_THICK = GAP
EXPJOINT_HEIGHT = 5
EXPJOINT_MIN = -RAD_MAJ/2 - EXPJOINT_HEIGHT
EXPJOINT_MAX = EXPJOINT_MIN + EXPJOINT_HEIGHT
