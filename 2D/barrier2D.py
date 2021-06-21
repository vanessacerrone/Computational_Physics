import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Rectangle
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve

# Conventions: hbar = m = 1

# Wave function for the initial time t=0

def psi0(x, y, x0, y0, sigma=0.5, kx=15*np.pi, ky = 0):

    '''

     function that return gaussian wave packet with predefined input parameters 
     we assume the electron not to have initial momentum in the y direction -> ky=0

    '''
    
    return np.exp(-1/2*((x-x0)**2 + (y-y0)**2)/sigma**2)*np.exp(-1j*(kx*(x-x0))+ky*(y))

def main():
    
    # =============================================================================
    # Parameters 
    # =============================================================================

    L = 8 # Well of width L
    dy = 0.1 # Spatial step size
    dx = dy
    dt = dy**2/4 # Temporal step size
    Nx = int(L/dx) + 1 # Number of points on the x-axis
    Ny = int(L/dy) + 1 # Number of points on the y-axis
    Nt = 800 # Number of time steps

    alpha_x = -dt/(2j*dy**2) # Constant to simplify expressions
    alpha_y = -dt/(2j*dy**2) # Constant to simplify expressions

    # Initial position of the Gaussian wave function
    x0 = L/5
    y0 = L/2

    # Parameters of the double slit
    w = 0.6 # Width of the walls of the double slit/amplitude of the barrier in x
    s = 2 # Separation between the edges of the slits/amplitude of the barrier in y  ( ≈ 2/3 for pot. barr., o(0.1) for slits)
    a = 0.2 # Aperture of the slits

    ## Indexes that parameterize the double slit 
    # Horizontal axis
    x1 = int(1/(2*dy)*(L-w)) # Left edge
    x2 = int(1/(2*dy)*(L+w)) # Right edge
    # Vertical axis
    y1 = int(1/(2*dy)*(L+s) + a/dy) # Lower edge of the lower slit
    y2 = int(1/(2*dy)*(L+s))        # Upper edge of the lower slit
    y3 = int(1/(2*dy)*(L-s))        # Lower edge of the upper slit
    y4 = int(1/(2*dy)*(L-s) - a/dy) # Upper edge of the upper slit

    # Value of the potential (if finite)
    v0 = 20000
    v = np.zeros((Ny,Ny), complex) 
        
    Ni = (Nx-2)*(Ny-2)  # Number of unknown factors v[i,j], i = 1,...,Nx-2, j = 1,...,Ny-2

    x = np.linspace(0, L, Ny-2) # Array of spatial points
    y = np.linspace(0, L, Ny-2) # Array of spatial points

    x, y = np.meshgrid(x, y)

    psis = [] # To store the wave function at each time step

    psi = psi0(x, y, x0, y0) # Initialise the wave function with the Gaussian
    psi[0,:] = psi[-1,:] = psi[:,0] = psi[:,-1] = 0 # Boundary conditions at the edges of the simulation box (infinite potential well)
    psis.append(np.copy(psi)) # Store the wave function of the initial time step


    print('\nInsert 0 for potential barrier')
    print('\nInsert 1 for single slit')
    print('\nInsert 2 for double slit')
    print('\nInsert 3 for double slit & hard wall')


    # Read input from keyboard
    try:
        n=int(input('\n'))
    except ValueError:
        print('\nMust be a number')

    # Enter the appropriate function
    if n == 0:
        v[y3:y2,x1:x2] = v0
        name = 'PotBar'
    
    elif n == 1:
        v = np.zeros((Ny,Ny), complex) 
        #v[0:y4, x1:x2] = v0
        #v[y3:y2,x1:x2] = 0
        #v[y1:,  x1:x2] = v0
        name = '1Slit_HW'

    elif n == 2:
        v[0:y4, x1:x2] = v0
        v[y3:y2,x1:x2] = v0
        v[y1:,  x1:x2] = v0
        name = '2Slit'

    elif n == 3:
        v = np.zeros((Ny,Ny), complex) 
        name = '2Slit_HW'
    else:
        print('Try again')




    # =============================================================================
    # Construction of the matrices of the system of equations
    # =============================================================================

    # Matrices for the Crank-Nicolson calculus. The problem A·x[n+1] = b = M·x[n] will be solved at each time step
    A = np.zeros((Ni,Ni), complex)
    M = np.zeros((Ni,Ni), complex)

    # Fill the A and M matrices
    for k in range(Ni):     
        
        # k = (i-1)*(Ny-2) + (j-1)
        i = 1 + k//(Ny-2)
        j = 1 + k%(Ny-2)
        
        # Main central diagonal
        A[k,k] = 1 + 2*alpha_x + 2*alpha_y + 1j*dt/2*v[i,j]
        M[k,k] = 1 - 2*alpha_x - 2*alpha_y - 1j*dt/2*v[i,j]
        
        if i != 1: # Lower lone diagonal
            A[k,(i-2)*(Ny-2)+j-1] = -alpha_y 
            M[k,(i-2)*(Ny-2)+j-1] = alpha_y
            
        if i != Nx-2: # Upper lone diagonal
            A[k,i*(Ny-2)+j-1] = -alpha_y
            M[k,i*(Ny-2)+j-1] = alpha_y
        
        if j != 1: # Lower main diagonal
            A[k,k-1] = -alpha_x 
            M[k,k-1] = alpha_x 

        if j != Ny-2: # Upper main diagonal
            A[k,k+1] = -alpha_x
            M[k,k+1] = alpha_x



    # =============================================================================
    # Solve Schrödinger's equation
    # =============================================================================

    A_s = csc_matrix(A) # Compressed Sparse Column format

    # Solve the matrix system at each time step in order to obtain the wave function

    for i in range(1,Nt):
        psi_vect = psi.reshape((Ni)) # We adjust the shape of the array to generate the matrix b of independent terms
        b = np.matmul(M,psi_vect) # Calculate the array of independent terms (matmul = matrix product of two arrays)
        psi_vect = spsolve(A_s,b) # Solve the system for this time step
        psi = psi_vect.reshape((Nx-2,Ny-2))  # Restore the shape of the wave function array 

        if n == 3 :    # If infinite potential wall we cancel the wave function inside the walls of the double slit
            psi[0:y4, x1:x2] = 0
            psi[y3:y2,x1:x2] = 0
            psi[y1:,  x1:x2] = 0
        
        if n == 1 :    # If infinite potential wall we cancel the wave function inside the walls of the single slit
            psi[0:y4, x1:x2] = 0
            psi[y1:,  x1:x2] = 0

        psis.append(np.copy(psi)) # Save the result


    # Calculate the modulus of the wave function at each time step
    mod_psis = [] # For storing the modulus of the wave function at each time step
    for wavefunc in psis:
        re = np.real(wavefunc) # Real part
        im = np.imag(wavefunc) # Imaginary part
        mod = np.sqrt(re**2 + im**2) # We calculate the modulus
        mod_psis.append(mod) 
        # We save the calculated modulus in a 3D array defined by three indexes [n][i][j]
        # n = temporal step, i = y-value index on the grid, j = x-value index on the grid 

    
    

    ## If you need to save memory
    del psis
    del M
    del psi_vect
    del A
    del A_s
    del b
    del im 
    del re
    del psi
    
    # Change shape to save data in a txt file 
    # Asarray converts the input into an array 
    mod_psis_reshaped = np.asarray(mod_psis).reshape(np.asarray(mod_psis).shape[0], -1) 
    np.savetxt("{0}2D_mod_psis_data.txt".format(name), mod_psis_reshaped)
    
    return 
    

if __name__ == "__main__":
    main()
