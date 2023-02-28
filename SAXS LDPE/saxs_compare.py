import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import plot_functions as p_f

p_f.set_LaTeX_and_CMU(True)
p_f.set_font_size(13, 11.5, 10)

def gaussian(x, A, mu, sigma, C):
    return C + A * np.exp(-0.5*((x - mu) / sigma)**2)

def exponential(x, A, alpha, C):
    return C + A * np.exp(-x/alpha)

for direction in ["CD", "MD"]:

    fig = plt.figure(figsize=(8,6))
        
    # Read SAXS data
    file_deformed = f"SAXS LDPE/data/{direction}_deformed_SAXS.grad"
    data_deformed = np.genfromtxt(file_deformed, skip_header=20, skip_footer=20, delimiter=", ")
    q_deformed = data_deformed[:,0]
    I_deformed = data_deformed[:,1]

    file_undeformed = f"SAXS LDPE/data/{direction}_undeformed_SAXS.grad"
    data_undeformed = np.genfromtxt(file_undeformed, skip_header=20, skip_footer=20, delimiter=", ")
    q_undeformed = data_undeformed[:,0]
    I_undeformed = data_undeformed[:,1]

    # Scale to same intensity
    I_deformed *= I_undeformed[0] / I_deformed[0]

    # Fit exponential
    i_fit = (q_undeformed < 0.03) | (q_undeformed > 0.15)
    param, cov = curve_fit(exponential, np.log10(q_undeformed[i_fit]), np.log10(I_undeformed[i_fit]))
    A, alpha, C = param

    # For plotting fit
    q_fit = np.linspace(0.01, 0.3)
    I_fit = 10**(exponential(np.log10(q_fit), A, alpha, C))

    # Plot inset with raw data
    ax1 = fig.add_subplot(3,3,1)
    ax1.loglog(q_undeformed, I_undeformed, label="Undeformed")
    ax1.loglog(q_deformed, I_deformed, label="Deformed")
    ax1.loglog(q_fit, I_fit, "k--", label="Exponential fit")
    ax1.set_title("Raw SAXS data")

    # Subtract intensity in log space
    I_deformed = 10**(np.log10(I_deformed) - exponential(np.log10(q_deformed), A, alpha, C))
    I_undeformed = 10**(np.log10(I_undeformed) - exponential(np.log10(q_undeformed), A, alpha, C))

    # Fit gaussians
    i_fit = (q_undeformed > 0.03) & (q_undeformed < 0.15)
    param, cov = curve_fit(gaussian, q_undeformed[i_fit], I_undeformed[i_fit])
    A, mu, sigma, C = param
    I_gauss_undeformed = gaussian(q_fit, A, mu, sigma, C)
    max_undeformed = C + A
    q_mid_undeformed = mu
    FWHM_undeformed = 2.355 * sigma
    d_undeformed = 2 * np.pi / q_mid_undeformed

    i_fit = (q_deformed > 0.03) & (q_deformed < 0.15)
    param, cov = curve_fit(gaussian, q_deformed[i_fit], I_deformed[i_fit])
    A, mu, sigma, C = param
    I_gauss_deformed = gaussian(q_fit, A, mu, sigma, C)
    max_deformed = C + A
    q_mid_deformed = mu
    FWHM_deformed = 2.355 * sigma
    d_deformed = 2 * np.pi / q_mid_deformed

    # Plot main plot with fits
    ax2 = fig.add_subplot(3,3,(4,9))
    ax2.plot(q_undeformed, I_undeformed, label=f"LDPE, {direction}, Undeformed")
    ax2.annotate(
        f"$q$ = {q_mid_undeformed:.3f}" + "$\,$Å$^{-1}$\n" +
        f"$d$ = {0.1*d_undeformed:.1f}" + "$\,$nm", 
        (q_mid_undeformed + 0.3*FWHM_undeformed, max_undeformed-0.5))
    ax2.plot(q_deformed, I_deformed, label=f"LDPE, {direction}, Deformed")
    ax2.annotate(
        f"$q$ = {q_mid_deformed:.3f}" + "$\,$Å$^{-1}$\n" +
        f"$d$ = {0.1*d_deformed:.1f}" + "$\,$nm", 
        (q_mid_deformed + 0.5*FWHM_deformed, max_deformed-0.5))
    ax2.plot(q_fit, I_gauss_undeformed, "k--", label="Fitted curve")
    ax2.plot(q_fit, I_gauss_deformed, "k--")
    ax2.set_xlabel("Momentum transfer, $q$ "  + "(Å$^{-1}$)")
    ax2.set_ylabel("Intensity (arb. unit)")
    ax2.legend(loc="upper left")
    ax2.grid()
    ax2.set_ylim(0.5, 6)

    fig.tight_layout()
    fig.savefig(f"SAXS LDPE/output/LDPE_{direction}_SAXS.pdf")