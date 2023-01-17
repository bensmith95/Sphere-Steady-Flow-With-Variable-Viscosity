# Sphere-Steady-Flow-With-Variable-Viscosity

MATLAB code that solves and compiles the steady flow around a rotating sphere with variable viscosity of the form discussed in [Miller _et al_ (2020)](https://doi.org/10.1063/1.5129220). The package contains code that solves the heated, rotating sphere boundary layer, the new model of the impinging region and the subsequent eruption of the radial jet.

All numerical methods are the same as those in the [Steady Sphere Flow](https://github.com/bensmith95/Sphere-Steady-Flow) code. Full details can be seen in Smith _et al_ (2023) [in draft].

To run the code, download all the files making sure the directory structure remains intact. Then simply enter _BaseFlow(**Re**,**Pr**,**&lambda;**)_ in the command window, where _**Re**_ is the square root of the Reynolds number, **Pr** is the Prandtl number (~0.72 is the value for air) and **&lambda;** is a sensitivity constant that represents a change in fluid (if set to 0 then the temperture independent case is solved). The data will be saved in a directory called _Flows_ that will be created automatically. To see the results run the file _figures.m_.

<p align="center">
  <img width="650" src="https://user-images.githubusercontent.com/29705711/212893164-1682e079-f5c9-45a6-9103-eeadd5f8255c.png">
</p>
