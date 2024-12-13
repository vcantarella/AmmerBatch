
"""
    ammer_ode(du, u, p, t)

RHS ODE function, describing the reactions in the batch reactors.
Done according to requirements to solvers in DifferentialEquations.jl.

# Arguments
- `p`: Store model parameters. These are:
    `r_max = p[1]`: Maximum reaction rate.
    `K_doc = p[2]`: Half-saturation constant for DOC.
    `K_nit = p[3]`: Half-saturation constant for nitrate.
    `K_nitam = p[4]`: Half-saturation constant for nitrite.
    `K_n2o = p[5]`: Half-saturation constant for N2O.
- `kwargs...`: Description of keyword arguments.

# Returns
we update the differential equation (du) in place.

# Examples
```julia
du = zeros(3)
u = [1.0, 2.0, 3.0]
p = [1.0, 2.0, 3.0, 4.0, 5.0]
t = 0.0
AmmerModel(du, u, p, t)
```
"""
function dummy_project_function(x, y)
    return x + y
end
