import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import make_axes_locatable

def main():

    # =============================================================================
    # Needed parameters 
    # =============================================================================

    L = 8 # Well width 
    dy = 0.1 # Spatial step size 
    dt = dy**2/4 # Temporal step size, condition for stability 
    Nx = int(L/dy) + 1 # Number of points on the x axis
    Ny = int(L/dy) + 1 # Number of points on the y axis
    Nt = 800 # Number of time steps

    ## Parameters of the double slit
    w = 0.6 # Width of the walls of the double slit/amplitude of the barrier in x
    s = 2 # Separation between the edges of the slits/amplitude of the barrier in y
    a = 0.2 # Aperture of the slits

    ## Indexes that parameterize the double slit in the space of points
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
    loaded_mod_psis = np.loadtxt("PotBar2D_mod_psis_data.txt")
    
    # The loaded_mod_psis array is an auxiliary 2D array, we need to return it to its original form
    mod_psisshape2 = Ny-2

    # Obtain mod_psis array
    mod_psis = loaded_mod_psis.reshape( 
        loaded_mod_psis.shape[0], loaded_mod_psis.shape[1] // mod_psisshape2, mod_psisshape2) 
    
    ## For deleting the auxiliary 2D array
    # del loaded_mod_psis


    # Create the figure
    fig = plt.figure(figsize = (8 , 8)) 
    ax = fig.add_subplot(xlim=(0,L), ylim=(0,L)) # Add the subplot to the figure
    ax.set_ylabel('y[a.u]',fontsize = 16)
    ax.set_xlabel('x[a.u]',fontsize = 16)
    img = ax.imshow(mod_psis[0], extent=[0,L,0,L], cmap=plt.get_cmap("plasma"), vmin=0, vmax=np.max(mod_psis)*0.5, zorder=1, interpolation="sinc") # Represent the modulus of the 2D wave function 
    
    # Create an axes on the right side of ax. The width of cax will be 5% of ax and the padding between cax and ax will be fixed at 0.05 inch
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)


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
        wall_middle = Rectangle((x1*dy,y3*dy), w, (y2-y3)*dy, color="white", zorder=50, alpha=1) 
        ax.add_patch(wall_middle)
        name = 'PotBar'
        title = 'Potential barrier' 

    elif n == 1:
        wall_bottom = Rectangle((x1*dy,0),     w, y4*dy,      color="white", zorder=50, alpha=1) # (x0, y0), width, height
        wall_top    = Rectangle((x1*dy,y1*dy), w, y4*dy,      color="white", zorder=50, alpha=1)
        ax.add_patch(wall_bottom)  # Add the rectangular patches to the plot
        ax.add_patch(wall_top)
        name = '1Slit'
        title = 'Single slit with infinite potential barrier' 

    elif n == 2:
        wall_bottom = Rectangle((x1*dy,0),     w, y4*dy,      color="white", zorder=50, alpha=0.5) # (x0, y0), width, height
        wall_middle = Rectangle((x1*dy,y3*dy), w, (y2-y3)*dy, color="white", zorder=50, alpha=0.5)
        wall_top    = Rectangle((x1*dy,y1*dy), w, y4*dy,      color="white", zorder=50, alpha=0.5)

        ax.add_patch(wall_bottom)
        ax.add_patch(wall_middle)
        ax.add_patch(wall_top)
        name = '2Slit'
        title = 'Double slit with potential $V_0$' 

    elif n == 3:
        wall_bottom = Rectangle((x1*dy,0),     w, y4*dy,      color="white", zorder=50, alpha=1) # (x0, y0), width, height
        wall_middle = Rectangle((x1*dy,y3*dy), w, (y2-y3)*dy, color="white", zorder=50, alpha=1)
        wall_top    = Rectangle((x1*dy,y1*dy), w, y4*dy,      color="white", zorder=50, alpha=1)
        
        ax.add_patch(wall_bottom)
        ax.add_patch(wall_middle)
        ax.add_patch(wall_top)
        name = '2Slit_HW'
        title = 'Double slit with infinite potential barrier' 

    else:
        print('Try again')

    ax.set_title("{0}".format(title), fontsize = 16)
    cbar = fig.colorbar(img,cax=cax)

    # Define the animation function for FuncAnimation
 
    def animate(i):

        img.set_data(mod_psis[i]) # Fill img with the modulus data of the wave function
        img.set_zorder(1)
        
        return img, # Return the result 

    anim = FuncAnimation(fig, animate, interval=5, frames=Nt, repeat=False, blit=True) # Generate the animation
    
    
    anim.save('{0}_2D.mp4'.format(name), extra_args=['-vcodec', 'libx264'],dpi=150, fps=60)

    plt.show() # Show the animation.
    
    
    return

if __name__ == "__main__":
    main()



