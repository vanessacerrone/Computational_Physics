import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import gridspec
from matplotlib.ticker import FormatStrFormatter

def main():
    # =============================================================================
    # Needed parameters 
    # =============================================================================
    L = 8 # Well of width L
    dy = 0.1 # Spatial step size 
    dt = dy**2/4 # Temporal step size, for stability 
    Nx = int(L/dy) + 1 # Number of points on the x axis
    Ny = int(L/dy) + 1 # Number of points on the y axis
    Nt = 800 # Number of time steps.

    # Parameters of the double slit.
    w = 0.6 # Width of the walls of the double slit/amplitude of the barrier in x
    s = 0.2 # Separation between the edges of the slits/amplitude of the barrier in y
    a = 0.2 # Aperture of the slits.

    # Indices that parameterize the double slit in the space of points.
    # Horizontal axis
    x1 = int(1/(2*dy)*(L-w)) # Left edge
    x2 = int(1/(2*dy)*(L+w)) # Right edge
    # Vertical axis
    y1 = int(1/(2*dy)*(L+s) + a/dy) # Lower edge of the lower slit
    y2 = int(1/(2*dy)*(L+s))        # Upper edge of the lower slit
    y3 = int(1/(2*dy)*(L-s))        # Lower edge of the upper slit
    y4 = int(1/(2*dy)*(L-s) - a/dy) # Upper edge of the upper slit



    # =============================================================================
    # Animation
    # =============================================================================
    # Load saved data

    loaded_mod_psis = np.loadtxt("2Slit2D_mod_psis_data.txt")
    
    # The loaded_mod_psis array is an auxiliary 2D array, we need to return it to its original form

    mod_psisshape2 = Ny-2

    # Obtain mod_psis array

    mod_psis = loaded_mod_psis.reshape( 
        loaded_mod_psis.shape[0], loaded_mod_psis.shape[1] // mod_psisshape2, mod_psisshape2) 
    
    ## For deleting the auxiliary 2D array.
    # del loaded_mod_psis

    # Create a screen to show Y_projection
    # Posizione dello schermo 
    c = 6
    index = int(c/dy)
    yy = np.zeros((Ny-2,Nt))

    for n in np.arange(0,Nt):
        for i in np.arange(0, Ny-2): 
            yy[i][n] = (mod_psis[n][i][index])**2 
    
    # in order n = temporal step, i = y-value index on the grid, j = x-value index on the grid 



    # We create the figure
    
    fig = plt.figure(figsize=(16, 8)) 
    gs = fig.add_gridspec(10, 35) 
    ax = fig.add_subplot(gs[:,:21])
    ax2 = fig.add_subplot(gs[:,25:])

    ax.set_xlim(0,L)
    ax.set_ylim(0,L)
    ax.set_ylabel('y[a.u]',fontsize = 16)
    ax.set_xlabel('x[a.u]',fontsize = 16)
    ax.vlines(c, 0, L, color='white',ls = 'dashed',lw =2, zorder=2, label = 'Screen')
    img = ax.imshow(mod_psis[0], extent=[0,L,0,L], cmap=plt.get_cmap("jet"), vmin=0, vmax=np.max(mod_psis)*0.6, zorder=1, interpolation="sinc") # Here the modulus of the 2D wave function shall be represented
    
    # create an axes on the right side of ax. The width of cax will be 5%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch
    
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)

 
    # Plot settings for the screen 

    limy = np.amax(yy) # Adjust y-axis limit in order to have better visualization
    ax2.set_xlim(0, limy),
    ax2.set_ylim(0,L)  
    ax2.set_xlabel('$|\psi(y,t)|^2$',fontsize = 16)
    ax2.set_ylabel('y [a.u]',fontsize = 16)
    ax2.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    line, = ax2.plot([], [], lw=2, color = '#3339FF', label = '$|\psi(y,t)|^2$') # WF line
    title2 = ax2.set_title("")
    

    

    print('\nInsert 0 for potential barrier')
    print('\nInsert 1 for single slit')
    print('\nInsert 2 for double slit')
    print('\nInsert 3 for double slit & hard wall')


    # Read input 
    try:
        n=int(input('\n'))
    except ValueError:
        print('\nMust be a number')

    # Choose the right potential configuration
    if n == 0:
        # Paint the walls of the double slit with rectangles
        wall_middle = Rectangle((x1*dy,y3*dy), w, (y2-y3)*dy, color='white', zorder=50, alpha=0.5)
        ax.add_patch(wall_middle)
        name = 'PotBar'
        title = 'Potential barrier' 

    elif n == 1:
        wall_bottom = Rectangle((x1*dy,0),     w, y4*dy,      color='white', zorder=50, alpha=0.5) # (x0, y0), width, height
        wall_top    = Rectangle((x1*dy,y1*dy), w, y4*dy,      color='white', zorder=50, alpha=0.5)
        ax.add_patch(wall_bottom)  # Add the rectangular patches to the plot
        ax.add_patch(wall_top)
        name = '1Slit'
        title = 'Single slit with potential $V_0$' 

    elif n == 2:
        wall_bottom = Rectangle((x1*dy,0),     w, y4*dy,      color='white', zorder=50, alpha=0.5) # (x0, y0), width, height
        wall_middle = Rectangle((x1*dy,y3*dy), w, (y2-y3)*dy, color='white', zorder=50, alpha=0.5)
        wall_top    = Rectangle((x1*dy,y1*dy), w, y4*dy,      color='white', zorder=50, alpha=0.5)

        ax.add_patch(wall_bottom)
        ax.add_patch(wall_middle)
        ax.add_patch(wall_top)
        name = '2Slit'
        title = 'Double slit with potential $V_0$' 

    elif n == 3:
        wall_bottom = Rectangle((x1*dy,0),     w, y4*dy,      color='white', zorder=50, alpha=1, fill = True ) # (x0, y0), width, height
        wall_middle = Rectangle((x1*dy,y3*dy), w, (y2-y3)*dy, color='white', zorder=50, alpha=1, fill = True)
        wall_top    = Rectangle((x1*dy,y1*dy), w, y4*dy,      color='white', zorder=50, alpha=1, fill = True)
        
        ax.add_patch(wall_bottom)
        ax.add_patch(wall_middle)
        ax.add_patch(wall_top)
        name = '2Slit_HW'
        title = 'Double slit with infinite potential barrier' 

    else:
        print('Try again')

    ax.set_title("{0}".format(title), fontsize = 16)
    cbar = fig.colorbar(img,cax=cax)

    # Define the animation function for FuncAnimation.
    
    def init():
        img.set_data(mod_psis[0])
        line.set_data([], [])
        title2.set_text("")
        return img, line, title2,

    def animate(i):

        img.set_data(mod_psis[i]) # Fill img with the modulus data of the wave function
        img.set_zorder(1)
        #wall_bottom = ax.fill_between( x =(x1*dy,x1*dy+w),x1 = y1*dy, x2 = y1*dy+y4*dy, color = "white")

        
        y = np.linspace(0, L, Ny-2)
        ppsi = yy[:,i] # Ad ogni iterazione prendo una colonna cioè un istante temporale 
        line.set_data(ppsi, y)
        time = dt * i
        title2.set_text('Elapsed time: {:6.2f} a.u.'.format(time))
   

        return img, line, title2,


       
    
    
    anim = FuncAnimation(fig, animate, interval=5,init_func=init, frames=Nt, repeat=False, blit=True) # Generate the animation
    

    anim.save('{0}_2D_wProjY.mp4'.format(name), extra_args=['-vcodec', 'libx264'],dpi=100, fps=40)
    
    plt.show() # Show the animation.
    
    
    return

if __name__ == "__main__":
    main()







