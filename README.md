# Numerical implementation of an invariant torus for a model of Tokamak.

Description of repository scripts:

* Control-of-Magnetic-Field-Lines.jl

It generates the initial torus $K_{0}(\theta) = (\theta + x_{p}(\theta), y_{p}(\theta))$ by integrating the vector field over time multiples of 2pi.

* Torus.jl

It Build the two-dimensional torus, for this a grid is generated with the coordinates $\theta$ and $\varphi$. The periodic components at the point $\theta_{i}$ and $\varphi_{j}$ of the torus are obtained by integrating the vector field for time $t = \varphi_{j}$ and initial condition $(\theta_{i} - frac{\omega \varphi_{j}}{\alpha} + xp{\theta_{i} - frac{\omega \varphi_{j}}{\alpha}}, xp{\theta_{i} - frac{\omega \varphi_{j}}{\alpha}})$. Where $\omega$ is the internal frequencies and $\alpha$ is the external frequencies.

* Graph-Torus.jl

Graph the one-dimensional torus.

* Auxiliary-Functions.jl

In this script are the necessary functions for the scripts described above.
