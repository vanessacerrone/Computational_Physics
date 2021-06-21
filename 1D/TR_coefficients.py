import numpy as np 
from matplotlib import pyplot as plt


# Function that computes Transmission and Reflection probabilities if function of time 

def main():
    
    dt = 1e-16
    N_steps = 2000
    a = 253 * 1e-10
    b = 255 * 1e-10
    Nx = 1000
    L = 500 * 1e-10

    # Load data
    data = np.loadtxt("psi_1D.txt")
    
    # Array to store transmission and reflection probability 
    P_t = np.zeros(N_steps)  # probability of WF in the x range [b,L]
    P_r = np.zeros(N_steps)  # probability of WF in the x range [0,a]

    # Indexes that parameterize the range limits in the space of points
    i_a = int((Nx-1)* a / (L))
    i_b = int((Nx-1)* b / (L))

    for n in np.arange(0,N_steps):
        P_r[n] = np.sum(data[0:i_a,n]) # Reflection probability
        P_t[n] = np.sum(data[i_b:,n])  # Transmission probability

    # Set the figure
    fig = plt.figure(figsize = (10 , 8))
    # Set x axis values 
    t = np.arange(0,N_steps) * dt * 1e15
    # Plot
    plt.plot(t , P_r , label = '$P_{reflect}$', color = '#FF4B00')
    plt.plot(t, P_t , label = '$P_{transmit}$', color = '#0451FF')
    plt.ylabel('Probability',fontsize = 16, loc = 'top')
    plt.xlabel('$t [fs]$',fontsize = 16, loc = 'right')
    plt.title('Transmission and reflection probabilities against time', fontsize = "20", pad=20 )
    #plt.tight_layout()
    plt.legend(fontsize = 16, loc = "best")
    plt.show()
    fig.savefig("TR_coefficients.png", dpi=400)
    return 


if __name__ == "__main__":
    main()

