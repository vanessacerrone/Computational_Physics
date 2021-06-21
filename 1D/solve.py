import numpy as np


# Function that solves an eigenvalue equation with tridiagonal matrix 

def solve(n,diag,upper,lower,b,sol): 
    '''
     function that solves Mx = b where M is a tridiagonal matrix 
     x - initially contains the input vector v, and returns the solution sol
     n - number of equations (length of vector x)
     lower - subdiagonal 
     diag -  main diagonal
     upper - superdiagonal  
     
    '''  
    # Array of zeros to store previously defined coefficients
    alfa = np.zeros(n, dtype=complex)
    beta = np.zeros(n, dtype=complex)
    
    # Initialize first term 
    alfa[0] = -diag[0]/upper[0]
    beta[0] = b[0]/upper[0]

    # Calculate the remaining values 
    for i in np.arange(1,n):
        alfa[i] = (-lower[i-1]/(upper[i]*alfa[i-1])-diag[i]/upper[i])
        beta[i] = b[i]/upper[i]+lower[i-1]*beta[i-1]/(upper[i]*alfa[i-1])
    
    # Calculate first term of the solution, that will be used to iterate the algorythm 
    sol[n-1]=(b[n-1]/lower[n-2]+beta[n-2]/alfa[n-2])/(1./alfa[n-2]+diag[n-1]/lower[n-2])

    # Perform the back substitution to calculate the solution
    for i in np.arange(n-1,0,-1): 
        sol[i-1]=sol[i]/alfa[i-1]-beta[i-1]/alfa[i-1]

    return(sol)