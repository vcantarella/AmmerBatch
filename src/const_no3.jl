using IfElse
"""
    const_denitrification(du, u, p, t, H)

This function represents the right-hand side (RHS) of an ordinary differential equation (ODE) in the context of the SciML ecosystem. It models the rates of change of concentrations of various substances over time.

# Arguments
- `du::Vector{Float64}`: A vector that will be modified in place to contain the rates of change of the concentrations of the substances at the current time step.
- `u::Vector{Float64}`: A vector containing the concentrations of the substances at the current time step. Specifically, it includes:
  - `u[1]`: Concentration of NO3- (nitrate)
  - `u[2]`: Concentration of DOC (dissolved organic carbon)
  - `u[3]`: Concentration of N2O (nitrous oxide) in water (mol/L)
  - `u[4]`: Gas concentration of N2O in ppmv (nitrous oxide)
  - `u[5]`: Volume of water (L)
  - `u[6]`: Volume of gas (L)
- `p::Vector{Float64}`: A vector of constant rates applied to the model. These parameters influence the rate of change of the concentrations.
  - `p[1]`: Rate of NO3- reduction
  - `p[2]`: Rate of N2O production
  - `p[3]`: Gas exchange rate
- `t::Float64`: The current time step.
- `H::Float64`: Henry's law constant for N2O.

# Returns
- `du::Vector{Float64}`: A vector containing the rates of change of the concentrations of the substances at the current time step.

# Example
"""
function const_denitrification(du, u, p, t, H)
    r_no3 = p[1]
    r_n2o = p[2]
    r_g = p[3]
    Vw = u[5]
    Vg = u[6]
    c_g = u[4] # gas concentration of N2O in ppmv
    c_w = u[3] # concentration of N2O in water (mol/L)
    gas_rate = r_g*(c_g*H - c_w)
    r_n2o = IfElse.ifelse(c_w > 0., r_n2o, 0.)
    r_no3 = IfElse.ifelse(u[1] > 0., r_no3, 0.)
    du[1] = -r_no3
    du[3] = 1/2*r_no3 - r_n2o + gas_rate/H
    du[4] = -gas_rate
end