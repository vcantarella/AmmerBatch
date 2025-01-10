"""
    const_denitrification(du, u, p, t)

This function represents the right-hand side (RHS) of an ordinary differential equation (ODE) in the context of the SciML ecosystem. It models the rates of change of concentrations of various substances over time.

# Arguments
- `du::Vector{Float64}`: A vector that will be modified in place to contain the rates of change of the concentrations of the substances at the current time step.
- `u::Vector{Float64}`: A vector containing the concentrations of the substances at the current time step. Specifically, it includes:
  - `u[1]`: Concentration of NO3- (nitrate)
  - `u[2]`: Concentration of DOC (dissolved organic carbon)
  - `u[3]`: Concentration of N2O (nitrous oxide)
  - `u[4]`: gas concentration of N2O (nitrous oxide)
  - `u[5]`: Mass of N2O (nitrous oxide)
- `p::Vector{Float64}`: A vector of constant rates applied to the model. These parameters influence the rate of change of the concentrations.
- `t::Float64`: The current time step.

# Returns
- `du::Vector{Float64}`: A vector containing the rates of change of the concentrations of the substances at the current time step.

# Example
"""
function const_denitrification(du, u, p, t)
    r_no3 = p[1]
    r_n2o = p[2]
    H = p[3]
    Vw = u[5]
    Vg = u[6]
    assert((Vw + Vg) ≈ 0.08)
    du[1] = -r_no3
    du[4] = (Vw+H*Vg)*(1/2*r_no3-r_n2o)
    du[3] = du[4]/(Vw+H*Vg)
    du[5] = du[3]*H
end