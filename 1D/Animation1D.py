from operator import pos
import numpy as np 
from matplotlib import pyplot as plt
from matplotlib import animation
from matplotlib.ticker import FormatStrFormatter

def main():
    # ============================================================================= 
    # Needed parameters 
    # =============================================================================
    
    h_bar = 1.055e-34
    m = 9.11e-31 # let's consider an electron
    # cell lenght and x-domain width (in meters)
    L = 500 * 1e-10
    # Initial wave-vector in reciprocal metres
    k0 = 3*1e10
    # Initial momentum & velocity 
    p0 = k0 * h_bar
    # Initial Energy  in J 
    E0 = (k0 * h_bar)**2 / (2 * m)

    # Potential barrier params
    # Specify potential barrier range in x coordinate 
    a = 253 * 1e-10
    b = 255 * 1e-10

    # Change 'f' to vary the potential barrier height V0 (in this way it's easy to check if E0 > V0 or E0 < V0)

    f  = 0.5
    V0 = f * E0


    # Simulations params 

    dt = 1e-16
    N_steps = 2000

    # =============================================================================
    # Animation
    # =============================================================================

    # Load saved data
    data = np.loadtxt("psi_1D.txt")
    positions = np.loadtxt("x_1D.txt")

    # Calculate energy and potential in eV
    V = V0 / (1.6*1e-19)
    E = E0 / (1.6*1e-19) 

    # Set the figure
    fig = plt.figure(figsize = (12 , 8))
    
    # Set the main plot 
    limy = np.amax(data)*0.4 # Adjust y-axis limit in order to have better visualization
    ax = plt.axes(xlim=(0, L), ylim=(0, limy))
    #ax.xaxis.set_major_formatter(FormatStrFormatter('%.0e'))
    ax.set_ylabel('$|\psi(x)|^2$',fontsize = 16)
    ax.set_xlabel('x [m]',fontsize = 16)
    line, = ax.plot([], [], lw=2, color = '#A600FF', label = '$|\psi(x)|^2$') # WF line
    title = ax.set_title("")
    
    # Create a second y-axis to visualize the potential 
    ax2 = ax.twinx()
    lim = E + 5 # Set y-axis limit to be the energy of the particle + 5eV
    ax2.set_ylim(0,lim)
    ax2.set_ylabel(r'$V\ (eV)$', fontsize = 16)

    # Draw the potential barrier
    ax2.vlines(a, 0, V, color='black', zorder=2, label = 'Potential barrier')
    ax2.vlines(b, 0, V, color='black', zorder=2)
    ax2.hlines(V, a, b, color='black', zorder=2)
    
    # Draw horizontal line to represent the energy of the particle 
    ax2.hlines(E,0,a,color='#07E230',linestyle = "dashed",lw= 1, zorder=2, label = 'Energy of the particle')
    plt.legend(loc = 'best', fontsize = 16)

    # Initialization function: plot the background of each frame
    def init():
        line.set_data([], [])
        title.set_text("")
        return line, title

    # Animation function
    def animate(i):
        x = positions 
        y = data[:,i] # Ad ogni iterazione prendo una colonna cioè un istante temporale 
        line.set_data(x, y)
        time = dt * i
        title.set_text('Elapsed time: {:6.2f} fs'.format(time * 1e+15))
        return line, title

    # Call the animator.  Blit=True means only re-draw the parts that have changed 
    anim = animation.FuncAnimation(fig, animate, init_func=init,frames=N_steps, interval=20, blit=True)

    anim.save('PotBarr_1D_{:.1f}.mp4'.format(f), fps=40, extra_args=['-vcodec', 'libx264'],dpi=200)
    plt.show()


    return 



if __name__ == "__main__":
    main()
