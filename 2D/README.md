# 2D time dependent Schrödinger equation
## Run the script that solves the eigenvalue equations and follow the istructions
Run the [Barrier2D.py](https://github.com/vanessacerrone/Computational_Physics/blob/main/1D/Barrier1D.py) module:

```
$ python3 Barrier2D.py
```
You can choose whether you want to simulate a finite potential barrier of predefined height, a single slit experiment with infinite potential barrier or a double slit experiment with either a finite or infinite potential. You will be asked to enter a number to choose the configuration: 

    - if n = 0: Potential barrier 
    - if n = 1: 1 Slit
    - if n = 2 : 2Slits with finite potential
    - if n = 3: 2Slits with hard walls 
  

The necessary data will be saved in txt files that will be loaded in the animation script. 
## Run the script to save the animation 
Run the [Animation2D.py](https://github.com/vanessacerrone/Computational_Physics/blob/main/1D/Animation1D.py) module:

```
$ python3 Animation1D.py
```

## Run the script to plot transmission and reflection coefficients 
Run the [TR_coefficients.py](https://github.com/vanessacerrone/Computational_Physics/blob/main/1D/TR_coefficients.py) module:

```
$ python3 TR_coefficients.py
```

**_Warning_: Make sure you have defined the same parameters in the various scripts, e.g. the height and width of the barrier!**


