import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import plot_functions as p_f

p_f.set_LaTeX_and_CMU(True)
p_f.set_font_size(11, 9, 9)

def gaussian(x, A, mu, sigma, C):
    return C + A * np.exp(-0.5*((x - mu) / sigma)**2)


for sample in ["130C", "160C", "190C", "220C", "250C"]:
    
    file = f"SAXS Al7075/data/{sample}_Al7075_WAXS.grad"
    data = np.genfromtxt(file, skip_header=20, skip_footer=20, delimiter=", ")
    q = data[:,0]
    I = data[:,1]

    fig = plt.figure(figsize=(p_f.cm_2_inch(16),p_f.cm_2_inch(8)))
    ax = fig.add_subplot(1,1,1)
    #ax.set_title(f"{sample}")
    ax.plot(q, I, label="Al7075, aged 160C 2h")

    # Aluminium planes
    #Al_a = 4.0479
    #Al_planes = [
    #    [1, 0, 0],
    #]
    #for i, plane in enumerate(Al_planes):
    #    h, k, l = plane
    #    q_plane = 2 * np.pi * np.sqrt(h**2 + k**2 + l**2) / Al_a
    #    ax.vlines(q_plane, 0, 0.001, alpha=0.3)
    #    ax.annotate(f"Al ({plane[0]}{plane[1]}{plane[2]})",
    #        (q_plane+0.01, 0.0001), alpha=0.5)
        
    # MgZn2 planes
    MgZn2_a = 5.220
    MgZn2_c = 8.566
    MgZn2_planes = [
        [0, 0, 0, 1],
        [1, 0,"$\\bar{1}$", 1],
        [0, 0, 0, 2],
        [1, 0,"$\\bar{1}$", 0],
    ]
    for plane_no, plane in enumerate(MgZn2_planes):
        h, k, i, l = plane
        d_plane = 1 / np.sqrt(
            (4/3) * (h**2 + h * k + k**2) / MgZn2_a**2 + l**2/MgZn2_c**2
        )
        q_plane = 2 * np.pi / d_plane
        ax.vlines(q_plane, 0, 0.0007+0.00018*plane_no, "k", alpha=0.3)
        ax.annotate(f"MgZn$_2$({plane[0]}{plane[1]}{plane[2]}{plane[3]})",
            (q_plane+0.01, 0.0006+0.00018*plane_no), alpha=0.5)

    ax.set_xlabel("Momentum transfer, $q$ "  + "(Å$^{-1}$)")
    ax.set_ylabel("Intensity, $I$ (arb. unit)")
    ax.legend()
    ax.set_ylim(0, 0.0015)
    ax.set_xlim(0.0, 2.0)
    fig.tight_layout()
    fig.savefig(f"SAXS Al7075/output/{sample}_Al7075_WAXS.pdf")
