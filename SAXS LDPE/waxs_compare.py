import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import plot_functions as p_f

p_f.set_LaTeX_and_CMU(True)
p_f.set_font_size(11, 9, 9)

def gaussian(x, A, mu, sigma, C):
    return C + A * np.exp(-0.5*((x - mu) / sigma)**2)

for direction in ["CD", "MD"]:

    fig = plt.figure(figsize=(8,4))
        
    # Read SAXS data
    file_deformed = f"SAXS LDPE/data/{direction}_deformed_WAXS.grad"
    data_deformed = np.genfromtxt(file_deformed, skip_header=20, skip_footer=20, delimiter=", ")
    q_deformed = data_deformed[:,0]
    I_deformed = data_deformed[:,1]

    file_undeformed = f"SAXS LDPE/data/{direction}_undeformed_WAXS.grad"
    data_undeformed = np.genfromtxt(file_undeformed, skip_header=20, skip_footer=20, delimiter=", ")
    q_undeformed = data_undeformed[:,0]
    I_undeformed = data_undeformed[:,1]
        
    # Shift to same height
    I_deformed *= np.min(I_undeformed) / np.min(I_deformed)

    ax = fig.add_subplot(1,1,1)
    ax.plot(q_undeformed, I_undeformed, label=f"LDPE, {direction}, Undeformed")
    ax.plot(q_deformed, I_deformed, label=f"LDPE, {direction}, Deformed")

    # Peak guesses
    mu_guess    = np.array([1.53, 1.69, 1.96])
    fit_width   = np.array([0.05, 0.03, 0.04])

    for sample in ["undeformed", "deformed"]:

        if sample == "undeformed": 
            q = q_undeformed
            I = I_undeformed
        else:
            q = q_deformed
            I = I_deformed

        # Iterate over peaks
        for i, mu in enumerate(mu_guess):

            # Refine guess (center on max value)
            width = fit_width[i]
            i_max = (q < mu + width) & (q > mu - width) # where to look for peak
            q_max = q[np.argmax(I * i_max)]

            # Fit gaussian
            i_fit = (q < mu + width) & (q > mu - width)
            param, cov = curve_fit(gaussian, q[i_fit], I[i_fit],
                bounds=([0  , q_max-width, 0,        0], 
                        [0.1, q_max+width, 10*width, 0.1]))
            A, mu, sigma, C = param
            max = C + A
            q_mid = mu
            FWHM = - 2.355 * sigma
            d = 2 * np.pi / q_mid

            i_plot = np.abs(q - mu) < 3 * sigma
            label = None
            #if sample == "undeformed" and i == 0: label="Gaussian fit"
            #ax.plot(q[i_plot], gaussian(q[i_plot], A, mu, sigma, C), "k--", label=label)
            #ax.annotate(
            #    f"q = {mu:.3f} $\pm$ {sigma:.3f}",
            #    (mu-0.1, A + C + 0.001))

    # Calculate expected q for different peaks
    a = 7.414
    b = 4.942
    c = 2.5473
    # BCC structure factor: h+k+l=even
    planes = [
        [2, 0, 0],
        [1, 1, 0],
        [2, 1, 0],
    ]
    q_planes = np.zeros(len(planes))

    for i, plane in enumerate(planes):
        h, k, l = plane
        q_plane = 2 * np.pi * np.sqrt((h/a)**2 + (k/b)**2 + (l/c)**2)

        ax.vlines(q_plane, -0.005, 0.06, alpha=0.3)
        ax.annotate(
            f"({planes[i][0]}{planes[i][1]}{planes[i][2]})",
            (q_plane+0.005, -0.003),
            alpha=0.5)

    plot_functions.set_LaTeX_and_CMU(True)
    ax.set_xlabel("Momentum transfer, $q$ "  + "(Å$^{-1}$)")
    ax.set_ylabel("Intensity, (arb. unit)")
    ax.legend()
    ax.set_ylim(-0.005, 0.06)
    ax.set_xlim(1.4, 2.1)
    plt.show()

    fig.tight_layout()
    fig.savefig(f"SAXS LDPE/output/LDPE_{direction}_WAXS.pdf")