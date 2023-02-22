import numpy as np
import matplotlib.pyplot as plt
        
# Read SAXS data
for sample in ["Reference", "130C", "160C", "190C", "220C", "250C"]:
    file = f"SAXS Al7075/data/{sample}_Al7075_SAXS.grad"
    data = np.genfromtxt(file, skip_header=20, skip_footer=20, delimiter=", ")
    q = data[:,0]
    I = data[:,1]

    # Fit to slope
    i_fit = q < 0.02
    param_line = np.polyfit(np.log(q[i_fit]), np.log(I[i_fit]), 1, w=I[i_fit])
    I_fit = np.exp(np.polyval(param_line, np.log(q)))


    # Plot fit
    fig = plt.figure()
    ax = fig.add_subplot(2,1,1)
    ax.set_title(f"Al7075, {sample}")
    ax.loglog(q, I, label="SAXS data")
    ax.loglog(q, I_fit, "k--", label="Linear fit to log-log data")
    ax.set_ylabel("Intensity, $I$ (A.U.)")
    ax.grid()
    ax.legend()

    # Subtract slope from data to get peak
    I_sub = I - I_fit

    # Select data to fit gaussian
    i_gauss = (I_sub > 0) & (q > 0.02) & (q < 0.15)

    # Preform gaussian fit (fit parabola to log of data)
    param = np.polyfit(q[i_gauss], np.log(I_sub[i_gauss]), 2, w=I_sub[i_gauss]**2)
    a, b, c = param
    mu = - b / (2 * a)
    sigma = 1 / np.sqrt(- 2 * a)

    print(f"{sample}: ({mu:.4f} pm {sigma:.4f}) Å^-1")

    I_gauss_fit = np.exp(np.polyval(param, q))

    ax = fig.add_subplot(2,1,2)
    ax.plot(q, I_sub, label="SAXS data, linear fit subtracted")
    ax.plot(q[i_gauss], I_sub[i_gauss], "r", label="Data used for Gaussian fit")
    ax.plot(q, I_gauss_fit, "k--", 
            label=f"Gaussian fit, ({mu:.4f} $\pm$ {sigma:.4f}) " + "Å$^{-1}$")
    ax.set_xlabel("Momentum transfer, $q$ (Å$^{-1}$)")
    ax.set_ylabel("Intensity, $I$ (A.U.)")
    ax.grid()
    ax.legend()
    ax.set_ylim(-0.1, 0.4)

    fig.tight_layout()
    fig.savefig(f"SAXS Al7075/output/{sample}_Al7075_SAXS.pdf")