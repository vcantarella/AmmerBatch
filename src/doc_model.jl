"""
    doc_model(du, u, p, t)

Right-hand side function for an ODE system representing DOC and NO3- dynamics.

# Arguments
- `du::Vector{Float64}`: Derivative of the state vector `u`.
- `u::Vector{Float64}`: State vector containing concentrations of NO3- and DOC.
  - `u[1]`: Concentration of NO3- (mmol L⁻¹).
  - `u[2]`: Concentration of DOC (mmol L⁻¹).
  - `u[4]`: Concentration of sorbed DOC (mmol L⁻¹).
- `p::Vector{Float64}`: Parameter vector.
  - `p[1]`: Zero-order rate of NO3- reduction (mmol L⁻¹ day⁻¹).
  - `p[2]`: Maximum rate of DOC denitrification by Michaelis-Menten kinetics (mmol L⁻¹ day⁻¹).
  - `p[3]`: First-order rate of DOC transfer to sorbed phase (day⁻¹).
  - `p[4]`: Equilibrium concentration of labile DOC in water (mmol L⁻¹).
  - `p[5]`: Half-saturation constant for DOC (mmol L⁻¹).
  - `p[6]`: Half-saturation constant for NO3- (mmol L⁻¹).
- `t::Float64`: Time (days).

# Description
This function calculates the time derivatives of the concentrations of NO3-, DOC, and sorbed DOC based on the provided rate equations. It is intended to be used with an ODE solver from the SciML ecosystem.

# Equations
- `r_no3`: Zero-order rate of NO3- reduction.
- `r_doc`: Michaelis-Menten kinetics for DOC denitrification.
- `r_transfer`: First-order rate of DOC transfer to sorbed phase.

# Returns
- Updates the `du` vector with the calculated derivatives:
  - `du[1]`: Derivative of NO3- concentration.
  - `du[2]`: Derivative of DOC concentration.
  - `du[4]`: Derivative of sorbed DOC concentration.
"""
function doc_model(du, u, p, t)
    no3_, doc = u[[1, 2]]
    doc_s = u[4] # sorbed DOC
    # define the rate constants
    r_no3 = p[1] # zero-order rate of NO3- reduction (mmol L-1 day-1)
    r_doc = p[2] # max rate of DOC denitrification by Michaelis-Menten kinetics (mmol L-1 day-1)
    αˡ = p[3] # first-order rate of DOC transfer to sorbed phase (day-1)
    Kd = p[4] # equilibrium concentration of labile DOC in water (mmol L-1)
    K_doc = p[5] # half-saturation constant for DOC (mmol L-1)
    K_no3 = p[6] # half-saturation constant for NO3- (mmol L-1)
    c_eq = Kd*doc_s
    # define the rate equations
    r_doc = r_doc * doc / (K_doc + doc) * no3_ / (K_no3 + no3_)
    r_transfer = αˡ * (c_eq - doc) * tanh(doc_s)
    du[1] = -r_no3 - r_doc
    du[2] = -5/4 * r_doc + r_transfer
    du[3] = 0
    du[4] = -r_transfer
end