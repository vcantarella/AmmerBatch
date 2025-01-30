"""
    doc_model(du, u, p, t)

Right-hand side function for an ODE system representing DOC and NO3- dynamics.

# Arguments
- `du::Vector{Float64}`: Derivative of the state vector `u`.
- `u::Vector{Float64}`: State vector containing concentrations and volumes.
  - `u[1]`: Concentration of NO3- (mmol L⁻¹)
  - `u[2]`: Concentration of DOC (mmol L⁻¹)
  - `u[3]`: Concentration of N2O in water (mmol L⁻¹)
  - `u[4]`: Concentration of N2O in gas phase (atm)
  - `u[5]`: Volume of water phase (L)
  - `u[6]`: Volume of gas phase (L)
  - `u[7]`: Concentration of sorbed DOC (mmol L⁻¹)
  - `u[8]`: Mass of soil/sediment (g)
- `p::Vector{Float64}`: Parameter vector.
  - `p[1]`: Zero-order rate of NO3- reduction (mmol L⁻¹ day⁻¹)
  - `p[2]`: Maximum rate of DOC denitrification (mmol L⁻¹ day⁻¹)
  - `p[3]`: First-order rate of DOC transfer to sorbed phase (day⁻¹)
  - `p[4]`: Equilibrium concentration of labile DOC in water (mmol L⁻¹)
  - `p[5]`: Half-saturation constant for DOC (mmol L⁻¹)
  - `p[6]`: Half-saturation constant for NO3- (mmol L⁻¹)
  - `p[7]`: Rate of N2O reduction by matrix DOC (mmol L⁻¹ day⁻¹)
  - `p[8]`: Rate of N2O reduction by soluble DOC (mmol L⁻¹ day⁻¹)
  - `p[9]`: Half-saturation constant for N2O (mmol L⁻¹)
"""
function doc_model!(du, u, p, t, H)
    no3_, doc = u[[1, 2]]
    doc_s = u[7] # sorbed DOC
    c_w = u[3] # concentration of N2O in water (mol/L)
    c_g = u[4] # gas concentration of N2O in atm
    Vw = u[5]
    Vg = u[6]
    mg = u[8]
    # define the rate constants
    k_no3 = p[1] # zero-order rate of NO3- reduction (mmol L-1 day-1)
    k_doc = p[2] # max rate of DOC denitrification by Michaelis-Menten kinetics (mmol L-1 day-1)
    αˡ = p[3] # first-order rate of DOC transfer to sorbed phase (day-1)
    Kd = p[4] # equilibrium concentration of labile DOC in water (mmol L-1)
    K_doc = p[5] # half-saturation constant for DOC (mmol L-1)
    K_no3 = p[6] # half-saturation constant for NO3- (mmol L-1)
    k_n2o = p[7] # rate of N2O reduction by matrix born DOC (mmol L-1 day-1)
    k_n2o_doc = p[8] # rate of N2O reduction by solluble bioavail. DOC (mmol L-1 day-1)
    K_n2o = p[9] # half-saturation constant for N2O (mmol L-1)
    c_eq = Kd*doc_s
    # define the rate equations
    r_doc = k_doc * doc / (K_doc + doc) * no3_ / (K_no3 + no3_)
    r_transfer = αˡ * (c_eq - doc_s)
    r_g = 1e-1
    r_no3 = ifelse(no3_ > 0.0, k_no3, 0.0)
    gas_rate = r_g*(c_g*H - c_w)
    rate_n2o = k_n2o*c_w/(K_n2o + c_w)
    rate_n2o_doc = k_n2o_doc * doc/(K_doc + doc) *  c_w/(K_n2o + c_w)
    du[1] = -r_no3 - r_doc
    du[2] = -r_doc - r_transfer*mg/Vw - 1/2*rate_n2o_doc
    du[3] = 1/2*(r_no3+r_doc) - rate_n2o - rate_n2o_doc + gas_rate
    du[4] = -gas_rate/H
    du[7] = r_transfer
end