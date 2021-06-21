# 2D time dependent Schrödinger equation
## Run the script that solves the eigenvalue equations and follow the istructions
Run the [Barrier2D.py](https://github.com/vanessacerrone/Computational_Physics/blob/main/2D/barrier2D.py) module:

```
$ python3 Barrier2D.py
```
You can choose whether you want to simulate a finite potential barrier of predefined height, a single slit experiment with infinite potential barrier or a double slit experiment with either a finite or infinite potential. You will be asked to enter a number to choose the configuration: 

    - if n = 0: potential barrier 
    - if n = 1: 1 slit
    - if n = 2 : 2 slits with finite potential
    - if n = 3: 2 slits with hard walls 
  

The necessary data will be saved in txt files that will be loaded in the animation script. 
## Run the script to save the animation 
Run the [Animation2D.py](https://github.com/vanessacerrone/Computational_Physics/blob/main/2D/Animation2D.py) module:

```
$ python3 Animation2D.py
```
You will be asked to enter the same number according to the chosen configuration. 
## Run the animation script in which we add a screen to visualize the y-projection of the wave function 
Run the [slits_projection.py](https://github.com/vanessacerrone/Computational_Physics/blob/main/2D/slits_projection.py) module:

```
$ python3 slits_projection.py
```

**_Warning_: Make sure you have defined the same parameters in the various scripts, e.g. the height and width of the barrier, and that you have chosen the same configurations!**


