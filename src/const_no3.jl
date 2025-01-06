"""
    const_denitrification(du, u, p, t)

This function represents the right-hand side (RHS) of an ordinary differential equation (ODE) in the context of the SciML ecosystem. It models the rates of change of concentrations of various substances over time.

# Arguments
- `du::Vector{Float64}`: A vector that will be modified in place to contain the rates of change of the concentrations of the substances at the current time step.
- `u::Vector{Float64}`: A vector containing the concentrations of the substances at the current time step. Specifically, it includes:
  - `u[1]`: Concentration of NO3- (nitrate)
  - `u[2]`: Concentration of DOC (dissolved organic carbon)
  - `u[3]`: Concentration of SO4-2 (sulfate)
- `p::Vector{Float64}`: A vector of constant rates applied to the model. These parameters influence the rate of change of the concentrations.
- `t::Float64`: The current time step.

# Returns
- `du::Vector{Float64}`: A vector containing the rates of change of the concentrations of the substances at the current time step.

# Example
"""
function const_denitrification(du, u, p, t)
    r_no3 = p[1]
    no3_ = u[1]
    du[1] = -r_no3
end