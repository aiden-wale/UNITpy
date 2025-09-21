
import uonidtoolbox as unit
import numpy as np
import matplotlib.pyplot as plt

# plt.rcParams.update({
#     "text.usetex": True,
#     "font.family": "serif"
# })

# TODO: deal with MIMO, MISO, SIMO
def showbode(G: list):

    nModels = len(G)
    if nModels < 1: return

    for i in range(0, nModels):
        G[i] = unit._utils.m2f(G[i])
    #endfor

    use_hertz   = False
    rad_2_hz    = 0.5/np.pi if use_hertz else 1.0

    fig, axs = plt.subplots(2, 1)
    for i in range(0, nModels):
        # Magnitude [dB]
        axs[0].semilogx(G[i].w*rad_2_hz, 20*np.log10(np.abs(G[i].G)))

        # Phase [deg]
        axs[1].semilogx(G[i].w*rad_2_hz, np.rad2deg(np.unwrap(np.angle(G[i].G))))
    #endfor

    axs[0].set_title("Estimated system: u_1 to y_1")
    axs[0].set_ylabel("Magnitude [dB]")
    axs[1].set_ylabel("Phase [deg]")
    x_axis_label = "Frequency [Hz]" if use_hertz else "Frequency [rad/s]"
    axs[1].set_xlabel(x_axis_label)

    plt.show()
#endfunction




