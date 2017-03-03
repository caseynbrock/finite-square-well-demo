import matplotlib.pyplot as plt
import square_well
import sys

E = float(sys.argv[1])
a = float(sys.argv[2])
V0 = float(sys.argv[3])

pseudo_pot = square_well.Potential(V0, a)
pseudo_wf = square_well.WaveFunction(E, pseudo_pot)
square_well.animate_scattering(pseudo_wf)
