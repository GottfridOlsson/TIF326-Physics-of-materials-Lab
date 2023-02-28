import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def gaussian(x, A, mu, sigma, C):
    return C + A * np.exp(-0.5*((x - mu) / sigma)**2)

for sample in ["130C", "160C", "190C", "220C", "250C"]:
    
    file = f"SAXS Al7075/data/{sample}_Al7075_WAXS.grad"
    data = np.genfromtxt(file, skip_header=20, skip_footer=20, delimiter=", ")
    q = data[:,0]
    I = data[:,1]

    # Calculate expected q for different peaks
    a = 4.0479
    # FCC structure factor: h, k, l all even or odd


    fig = plt.figure(figsize=(8,4))
    ax = fig.add_subplot(1,1,1)
    ax.set_title(f"{sample}")
    ax.plot(q, I, label="Al7075, Aged 160C 2h")

    # Aluminium planes
    Al_a = 4.0479
    Al_planes = [
        [1, 0, 0],
    ]
    for i, plane in enumerate(Al_planes):
        h, k, l = (1, 0, 0)
        q_plane = 2 * np.pi * np.sqrt(h**2 + k**2 + l**2) / Al_a
        ax.vlines(q_plane, 0, 0.001, alpha=0.3)
        ax.annotate(f"Al ({plane[0]}{plane[1]}{plane[2]})",
            (q_plane+0.01, 0.0001), alpha=0.5)
        
    # MgZn2 planes
    MgZn2_a = 5.16
    MgZn2_c = 5.16

    ax.set_xlabel("Momentum transfer, $q$ "  + "(Å$^{-1}$)")
    ax.set_ylabel("Intensity, $I$ (A.U.)")
    ax.legend()
    ax.set_ylim(0, 0.0015)
    #ax.set_xlim(0.5, 2.0)
    #fig.tight_layout()
    fig.savefig(f"SAXS Al7075/output/{sample}_Al7075_WAXS.pdf")