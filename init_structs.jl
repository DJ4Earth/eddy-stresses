# This script will contain the structures needed for a barotropic gyre model

mutable struct gyre_matrix
    u::Matrix{Float64}
    v::Matrix{Float64}
    η::Matrix{Float64}
    lastu_t::Matrix{Float64}
    lastv_t::Matrix{Float64}
    lastη_t::Matrix{Float64}
end

mutable struct gyre_vector
    u::Vector{Float64}
    v::Vector{Float64}
    eta::Vector{Float64}
end

# constants that appear in various places 
mutable struct Parameters 
    Tspinup::Int64          # how long to spin up from rest
    Tsteps::Int64           # how many time steps to take after spinning up the model
    dt::Float64             # timestep
    g::Float64              # gravity
    f0::Float64             # Coriolis parameter
    beta::Float64           # Coriolis parameter
    H::Float64              # ocean depth 
    A_h::Float64            # horizontal Laplacian viscosity 
    ρ_c::Float64            # reference density 
    eps_ab::Float64         # parameter pertaining to AB2
    bottom_drag::Float64    # bottom drag coefficient
    wind_stress::Vector{Float64}    # wind stress values on the grid 
    coriolis::Vector{Float64}       # coriolis values on the grid 
end


# parameters relating to the grid 
mutable struct Grid
    Lx::Int             # length of box in x direction [meters]
    Ly::Int             # length of box in y direction [meters]
    T::Int              # total number of steps to take after spinning up the model from rest
    nx::Int             # total number of grid cells in x direction
    ny::Int             # total number of grid cells in y direction
    NT::Int              # total number of cells on eta-grid 
    Nu::Int             # total number of cells on u-grid
    Nv::Int             # total number of cells on v-grid
    Nq::Int             # total number of cells on vorticity grid 
    dx::Float64 
    dy::Float64 
end 

# The gyre model written by Kloewer and used by Zanna leaves all the states as vectors stacked 
# row-wise. So, computing things like gradients, Laplacians, etc. become matrix operations, and 
# in this struct we store all of these operators for use in computing the RHS of the system 
# (namely the time derivatives)

mutable struct Derivatives
    GTx::Matrix{Float64}
    GTy::Matrix{Float64}
    Gux::Matrix{Float64}
    Guy::Matrix{Float64}
    Gvx::Matrix{Float64}
    Gvy::Matrix{Float64}
    Gqy::Matrix{Float64}
    Gqx::Matrix{Float64}
    Lu::Matrix{Float64}
    Lv::Matrix{Float64}
    LT::Matrix{Float64}
    Lq::Matrix{Float64}
end

mutable struct Interps 
    Ivu::Matrix{Float64}
    Iuv::Matrix{Float64}
    IqT::Matrix{Float64}
    IuT::Matrix{Float64}
    IvT::Matrix{Float64}
    ITu::Matrix{Float64}
    ITv::Matrix{Float64}
    Iqu::Matrix{Float64}
    Iqv::Matrix{Float64}
    Iuq::Matrix{Float64}
    Ivq::Matrix{Float64}
    ITq::Matrix{Float64}
end

mutable struct Advection
    AL1::Matrix{Float64}
    AL2::Matrix{Float64}
    index_av::Vector{Int64}
    index_bv::Vector{Int64}
    index_cv::Vector{Int64}
    index_dv::Vector{Int64}
    ALeur::Matrix{Float64}
    ALeul::Matrix{Float64}
    Seul::Matrix{Float64}
    Seur::Matrix{Float64}
    Sau::Matrix{Float64}
    Sbu::Matrix{Float64}
    Scu::Matrix{Float64}
    Sdu::Matrix{Float64}
    ALpvu::Matrix{Float64}
    ALpvd::Matrix{Float64}
    Spvu::Matrix{Float64}
    Spvd::Matrix{Float64}
    Sav::Matrix{Float64}
    Sbv::Matrix{Float64}
    Scv::Matrix{Float64}
    Sdv::Matrix{Float64}
end
