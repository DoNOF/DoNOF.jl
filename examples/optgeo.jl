using DoNOF

mol = """
0 1
 H   0.0 0.0  -0.404
 H   0.0 0.0   0.404
"""

bset, p = DoNOF.molecule(mol, "cc-pvdz", spherical = true)

p.RI = true

p.ipnof = 5
p.ista = 0

p.threshgorb = 10^-5

DoNOF.optgeo(mol, bset, p)
