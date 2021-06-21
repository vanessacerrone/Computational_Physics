import numpy as np 
from matplotlib import pyplot as plt
from matplotlib import animation
from solve import solve 

def main():


    # =============================================================================
    # Parameters 
    # =============================================================================

    
    # Fundamental constants
    h_bar = 1.055e-34
    m = 9.11e-31 # let's consider an electron
    # Cell lenght and x-domain width (in meters)
    L = 500 * 1e-10
    # Sigma in metres
    sigma = (20)*1e-10
    # Initial position of wave packet
    x0 = 150 * 1e-10
    # Initial wave-vector in reciprocal metres
    k0 = 3*1e10
    # Initial momentum & velocity 
    p0 = k0 * h_bar
    v0 = p0 / m
    # Initial Energy  in J 
    E0 = (k0 * h_bar)**2 / (2 * m)

    # Potential barrier parameters
    # Specify potential barrier range in x coordinate 
    a = 253 * 1e-10
    b = 255 * 1e-10

    # Change 'f' to vary the potential barrier height V0 (in this way it's easy to check if E0 > V0 or E0 < V0)
    f  = 1.1
    V0 = f * E0


    # X domain and simulations params 
    Nx = 1000
    h = L/(Nx-1) # There are Nx points and Nx -1 spaces between them 
    Nmat = Nx-2 #  Useful to set the matrix

    # Time domain
    dt = 1e-16
    N_steps = 2000


    
    psi0 = np.zeros(Nmat, dtype=complex) # Array to store psi at t = 0
    psi1 = np.zeros(Nmat, dtype=complex) # Array to store time-evolved psi

    # Arrays to store only the non-zeros elements of the matrix 
    diag = np.zeros(Nmat, dtype=complex)
    upper = np.zeros(Nmat, dtype=complex)
    lower = np.zeros(Nmat, dtype=complex)
    f = np.zeros(Nmat, dtype=complex) # solution 

    # Initialise the wave function with the Gaussian  

    norm=0.
    for i in np.arange(1,Nx-1):
        x = h*i
        psi0[i-1]= np.exp(k0*x*1.0j)*np.exp(-np.power(x-x0,2)/(2*np.power(sigma,2)))
        norm += np.power(np.abs(psi0[i-1]),2)
    norm = norm*L/Nx
    norm = np.sqrt(norm)

    # Normalize initial wave function 
    for i in np.arange(0,Nmat):
        psi0[i]=psi0[i]/norm

    
    # =============================================================================
    # Construction of the matrix M and solution of the system 
    # =============================================================================

    # Set entries above (upper) and below (lower) the main diagonal
    for i in np.arange(0,Nmat): 
        upper[i]=1.0
        lower[i]=1.0
        
    # Set entries in the main diagonal
    for i in np.arange(0,Nmat): 
        Vx = 0
        x = h*(i+1)
        if x >= a and x <= b:
            Vx = V0
        diag[i]=4.j*m*np.power(h,2)/(h_bar*dt) - 2. - 2*m*np.power(h,2)*Vx / (h_bar**2)


    positions = np.zeros(Nmat) # Array to store positions value on the x axis grid

    for i in np.arange(0,Nmat):
        x = h*(i+1)
        positions[i] = x

    # Array to store WF at each time step: each column corresponds to a time step (N_steps)
    # Each row corresponds to the WF at that time step (Nmat) 

    data = np.zeros((Nmat,N_steps)) 

    # Solve SC. equation for each time step
    for n in np.arange(0,N_steps):
        for i in np.arange(1,Nmat-1):
            Vx = 0
            x = h*(i+1)
            if x >= a and x <= b:
                Vx = V0
            f[i] = -psi0[i+1] + 2.*psi0[i] - psi0[i-1] + (4.j*m*np.power(h,2) / (h_bar*dt))*psi0[i] + 2.*m*np.power(h,2)*Vx*psi0[i] / (h_bar**2)

        f[0] = -psi0[1] + 2.*psi0[0] + (4.j*m*np.power(h,2) / (h_bar*dt))*psi0[0] 
        f[Nmat-1] = -psi0[Nmat-2] + 2.*psi0[Nmat-1] +(4.j*m*np.power(h,2) / (h_bar*dt))*psi0[Nmat-1]

        # Solve the equation with the tridiagonal matrix method  
        psi1 = solve(Nmat, diag, upper, lower,f, psi1)

        # Store the WF at this time step
        for i in np.arange(0,Nmat):
            psi0[i] = psi0[i] * np.power(L/Nx,0.5)  # We have to normalize!
            data[i][n]=np.power(np.abs(psi0[i]),2) 
            psi0[i] = psi1[i] # Update WF with the solution for the next step

        # Check if WF is normalized at each time step
        # norm = 0.
        # for i in np.arange(0,Nmat):
        #     norm+=pow(abs(psi0[i]),2)
        # norm=norm*L/Nx
        # print(n,norm)
    
    # Save data(WF) and x values in a txt file so it's easier to acces it 
    np.savetxt("psi_1D.txt", data)
    np.savetxt('x_1D.txt', positions)


    return 



if __name__ == "__main__":
    main()


