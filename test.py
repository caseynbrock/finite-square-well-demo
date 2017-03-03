import matplotlib.pyplot as plt
import square_well
import subprocess

V0=-2.0
a=2. 
E=0.5
pot = square_well.Potential(V0, a)
wf = square_well.WaveFunction(E, pot)
#square_well.animate_scattering(wf)
print  wf.find_a_pseudopotential(-2)
pseudo_V_list = wf.find_some_pseudopotentials(-100, 10, 10)

for V in pseudo_V_list:
    print "\nV0:  ", V
    subprocess.Popen(["python", "quick_plot.py", str(E), str(a), str(V)])
