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

    # Fit to slope
    
    # Guesses
    mu_guess    = np.array([0.75, 1.40, 1.48, 1.57])
    fit_width = np.array([0.05, 0.03, 0.03, 0.03])

    # Refine guess by centering on max value
    for i, mu in enumerate(mu_guess):
        width = fit_width[i]
        i_max = (q < mu + width) & (q > mu - width) # where to look for peak
        mu_guess[i] = q[np.sum(q <= mu - width) + np.argmax(I[i_max])]

    mu_fit    = np.zeros(4)
    sigma_fit = np.zeros(4)
    A_fit     = np.zeros(4)
    C_fit     = np.zeros(4)

    i_gauss_tot = np.zeros(len(q), dtype=np.bool)

    # Fit smaller gaussians
    for i, mu in enumerate(mu_guess):
        width = fit_width[i]
        
        # Select data to fit gaussian
        i_gauss = (q < mu + width) & (q > mu - width)
        i_gauss_tot += i_gauss

        # Preform gaussian fit
        A0      = 0.05
        mu0     = mu_guess[i]
        sigma0  = 0.01
        C0      = 0.01
        param, cov = curve_fit(gaussian, q[i_gauss], I[i_gauss], [A0, mu0, sigma0, C0])
        A, mu, sigma, C = param
        mu_fit[i]       = mu
        sigma_fit[i]    = sigma
        A_fit[i]        = A
        C_fit[i]        = C

    # Calculate expected q for different peaks
    a = 4.0479
    # FCC structure factor: h, k, l all even or odd
    planes = [
        #[1, 0, 0],
        #[0, 1, 0],
        #[0, 0, 1],

        #[2, 0, 0],
        #[0, 2, 0],
        #[0, 0, 2],

        #[3, 0, 0],
        #[0, 3, 0],
        #[0, 0, 3],

        #[0, 1, 1],
        #[1, 0, 1],
        #[1, 1, 0],

        #[1, 1, 1],

        #[2, 1, 0],
        #[2, 0, 1],
        #[0, 2, 1],
        #[1, 2, 0],
        #[1, 0, 2],
        #[0, 1, 2],

    ]
    q_planes = np.zeros(len(planes))

    for i, plane in enumerate(planes):
        h, k, l = plane
        q_planes[i] = 2 * np.pi * np.sqrt(h**2 + k**2 + l**2) / a

    # Plot
    fig = plt.figure(figsize=(16/2.54,9/2.54))
    ax = fig.add_subplot(1,1,1)

    ax.set_title(f"{sample}")
    ax.plot(q, I, label="WAXS data")
    ax.plot(q[i_gauss_tot], I[i_gauss_tot], "r.", label="WAXS data used for fit")
    
    for i, q_plane in enumerate(q_planes):
        ax.vlines(q_plane, 0, 0.001, alpha=0.3)
        ax.annotate(
            f"({planes[i][0]}{planes[i][1]}{planes[i][2]})",
            (q_plane+0.01, 0.0001),
            alpha=0.5)

    for i, mu in enumerate(mu_fit):
        sigma = sigma_fit[i]
        i_plot = np.abs(q - mu) < 5 * sigma
        label = None
        if i == 0: label="Gaussian fit"
        ax.plot(
            q[i_plot], 
            gaussian(q[i_plot], A_fit[i], mu, sigma, C_fit[i]), 
            "k--",
            label=label)
        ax.annotate(
            f"{mu:.3f}$\\,\\pm\\,${sigma:.3f}",
            (mu-0.01, A_fit[i] + C_fit[i] + 0.00002))

    ax.set_xlabel("Momentum transfer, $q$ "  + "(Å$^{-1}$)")
    ax.set_ylabel("Intensity, $I$ (A.U.)")
    #ax.grid()
    ax.legend()
    ax.set_ylim(0, 0.0015)
    ax.set_xlim(0.5, 2.0)

    fig.tight_layout()
    fig.savefig(f"SAXS Al7075/output/{sample}_Al7075_WAXS.pdf")
    #plt.show()