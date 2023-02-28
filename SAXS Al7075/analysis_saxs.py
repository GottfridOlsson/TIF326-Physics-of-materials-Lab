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

file_ref = f"SAXS Al7075/data/Reference_Al7075_SAXS.grad"
data_ref = np.genfromtxt(file_ref, skip_header=20, skip_footer=20, delimiter=", ")
q_ref = data_ref[:,0]
I_ref = data_ref[:,1]

# Read SAXS data
for sample in ["130C", "160C", "190C", "220C", "250C"]:
    file = f"SAXS Al7075/data/{sample}_Al7075_SAXS.grad"
    data = np.genfromtxt(file, skip_header=20, skip_footer=20, delimiter=", ")
    q = data[:,0]
    I = data[:,1]

    # Fit exponential
    i_fit = (q < 0.015) | (q > 0.2)
    param, cov = curve_fit(exponential, np.log10(q[i_fit]), np.log10(I[i_fit]))
    A, alpha, C = param

    # For plotting fit
    q_fit = np.linspace(0.01, 0.3)
    I_fit = 10**(exponential(np.log10(q_fit), A, alpha, C))

    # Plot fit
    fig = plt.figure(figsize=(8,6))
    ax1 = fig.add_subplot(3,3,1)
    ax1.set_title(f"Raw SAXS data")
    ax1.loglog(q_ref, I_ref, label="Reference")
    ax1.loglog(q, I, label=f"Aged, {sample} 2h")
    ax1.loglog(q_fit, I_fit, "k--", label="Fitted curve")


    I_sub = 10**(np.log10(I) - exponential(np.log10(q), A, alpha, C))

    # Fit gaussians
    i_fit = (q > 0.00) & (q < 0.5)
    param, cov = curve_fit(gaussian, q[i_fit], I_sub[i_fit],
        bounds=([0, 0, -100, 0], [10, 0.3, 0, 2]))
    A, mu, sigma, C = param
    I_gauss = gaussian(q_fit, A, mu, sigma, C)
    max = C + A
    q_mid = mu
    FWHM = - 2.355 * sigma
    d = 2 * np.pi / q_mid

    print(f"{sample}: d = {0.1*d:.1f} nm")

    ax2 = fig.add_subplot(3,3,(4,9))
    ax2.plot(0, 0, label="Al7075, no aging")
    ax2.plot(q, I_sub, label="Al7075, aged 160C 2h")
    ax2.plot(q_fit, I_gauss, "k--", label="Fitted curve")
    ax2.annotate(
        f"$q$ = {q_mid:.3f}" + "$\\,$Å$^{-1}$\n" +
        f"$d$ = {0.1*d:.1f}" + "$\\,$nm", 
        (q_mid + 0.4*FWHM, max-0.5))
    ax2.set_xlabel("Momentum transfer, $q$ (Å$^{-1}$)")
    ax2.set_ylabel("Intensity (arb. unit)")
    ax2.grid()
    ax2.legend(loc="upper left")
    ax2.set_ylim(0.5, 5)



    fig.tight_layout()
    fig.savefig(f"SAXS Al7075/output/{sample}_Al7075_SAXS.pdf")