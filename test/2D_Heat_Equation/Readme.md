# Benchmark: 2D Diffusion


## Problem description

This code simulates the anisotropic 2D heat equation,

$$\frac{\partial u}{\partial t} = \nabla \cdot (D \nabla u) + b(t, \mathbf{x})$$

where $D$ is a diagonal matrix with entries $k_x$ and $k_y$. The system is
evolved for $t \in [0, t_f]$ on the rectangular domain
$(x,y) \equiv \mathbf{x} \in [\mathbf{0}, \mathbf{x}_{\text{max}}]^2$, with the
initial condition

$$u(0,\mathbf{x}) = \sin^2(\pi x) \sin^2(\pi y),$$

and stationary boundary conditions

$$\frac{\partial u}{\partial t}(t,0,y) = \frac{\partial u}{\partial t}(t,x_{\text{max}},y) = \frac{\partial u}{\partial t}(t,x,0) = \frac{\partial u}{\partial t}(t,x,y_{\text{max}}) = 0.$$

The source term is given by

$$b(t,\mathbf{x}) = -2 \pi \sin^2(\pi x) \sin^2(\pi y) \sin(\pi t) \cos(\pi t) - k_x 2 \pi^2 (\cos^2(\pi x) - \sin^2(\pi x)) \sin^2(\pi y) \cos^2(\pi t) - k_y 2 \pi^2 (\cos^2(\pi y) - \sin^2(\pi y)) \sin^2(\pi x) \cos^2(\pi t).$$

Under this setup, the problem has the analytical solution

$$u(t,\mathbf{x}) = \sin^2(\pi x) \sin^2(\pi y) \cos^2(\pi t).$$