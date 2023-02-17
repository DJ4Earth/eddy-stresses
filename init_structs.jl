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
    dt::Float64             # timestep
    g::Float64              # gravity
    f0::Float64             # Coriolis parameter
    beta::Float64           # Coriolis parameter
    H::Float64              # ocean depth 
    A_h::Float64            # horizontal Laplacian viscosity 
    ρ_c::Float64            # reference density 
    bottom_drag::Float64    # bottom drag coefficient
    wind_stress::Vector{Float64}    # wind stress values on the grid 
    coriolis::Vector{Float64}       # coriolis values on the grid 
end


# parameters relating to the grid 
mutable struct Grid
    Lx::Int             # length of box in x direction [meters]
    Ly::Int             # length of box in y direction [meters]
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
    GTx::SparseMatrixCSC{Float64, Int64}
    GTy::SparseMatrixCSC{Float64, Int64}
    Gux::SparseMatrixCSC{Float64, Int64}
    Guy::SparseMatrixCSC{Float64, Int64}
    Gvx::SparseMatrixCSC{Float64, Int64}
    Gvy::SparseMatrixCSC{Float64, Int64}
    Gqy::SparseMatrixCSC{Float64, Int64}
    Gqx::SparseMatrixCSC{Float64, Int64}
    Lu::SparseMatrixCSC{Float64, Int64}
    Lv::SparseMatrixCSC{Float64, Int64}
    LT::SparseMatrixCSC{Float64, Int64}
    Lq::SparseMatrixCSC{Float64, Int64}
    LLu::SparseMatrixCSC{Float64, Int64}
    LLv::SparseMatrixCSC{Float64, Int64}
end

mutable struct Interps 
    Ivu::SparseMatrixCSC{Float64, Int64}
    Iuv::SparseMatrixCSC{Float64, Int64}
    IqT::SparseMatrixCSC{Float64, Int64}
    IuT::SparseMatrixCSC{Float64, Int64}
    IvT::SparseMatrixCSC{Float64, Int64}
    ITu::SparseMatrixCSC{Float64, Int64}
    ITv::SparseMatrixCSC{Float64, Int64}
    Iqu::SparseMatrixCSC{Float64, Int64}
    Iqv::SparseMatrixCSC{Float64, Int64}
    Iuq::SparseMatrixCSC{Float64, Int64}
    Ivq::SparseMatrixCSC{Float64, Int64}
    ITq::SparseMatrixCSC{Float64, Int64}
end

mutable struct Advection
    AL1::SparseMatrixCSC{Float64, Int64}
    AL2::SparseMatrixCSC{Float64, Int64}
    index_av::Vector{Int64}
    index_bv::Vector{Int64}
    index_cv::Vector{Int64}
    index_dv::Vector{Int64}
    ALeur::SparseMatrixCSC{Float64, Int64}
    ALeul::SparseMatrixCSC{Float64, Int64}
    Seul::SparseMatrixCSC{Float64, Int64}
    Seur::SparseMatrixCSC{Float64, Int64}
    Sau::SparseMatrixCSC{Float64, Int64}
    Sbu::SparseMatrixCSC{Float64, Int64}
    Scu::SparseMatrixCSC{Float64, Int64}
    Sdu::SparseMatrixCSC{Float64, Int64}
    ALpvu::SparseMatrixCSC{Float64, Int64}
    ALpvd::SparseMatrixCSC{Float64, Int64}
    Spvu::SparseMatrixCSC{Float64, Int64}
    Spvd::SparseMatrixCSC{Float64, Int64}
    Sav::SparseMatrixCSC{Float64, Int64}
    Sbv::SparseMatrixCSC{Float64, Int64}
    Scv::SparseMatrixCSC{Float64, Int64}
    Sdv::SparseMatrixCSC{Float64, Int64}
end

# per a suggestion, I'm creating a new structure that will pre-allocate space to operators that only appear 
# during the timestepping loop (mainly appear in the computation of the RHS equation)
mutable struct RHS_terms
    # Appear in advance
    u0::Vector{Float64}
    v0::Vector{Float64}
    eta0::Vector{Float64}
    u1::Vector{Float64}
    v1::Vector{Float64}
    eta1::Vector{Float64}
    # Appear in comp_u_v_eta_t
    h::Vector{Float64}          # height of water columns [meters]
    h_u::Vector{Float64}        # height of water columns interpolated to u-grid 
    h_v::Vector{Float64}        # height of water columns interpolated to v-grid 
    h_q::Vector{Float64}
    U::Vector{Float64}
    V::Vector{Float64}
    kinetic::Vector{Float64}
    kinetic_sq::Vector{Float64}
    q::Vector{Float64}
    p::Vector{Float64}
    bfric_u::Vector{Float64}
    bfric_v::Vector{Float64}
    Mu::Vector{Float64}
    Mv::Vector{Float64}
    rhs_u::Vector{Float64}
    rhs_v::Vector{Float64}
    rhs_eta::Vector{Float64}
    # Appear in comp_advection 
    AL1q::Vector{Float64}
    AL2q::Vector{Float64}
    AL1q_au::Vector{Float64}
    AL1q_du::Vector{Float64}
    AL1q_av::Vector{Float64}
    AL1q_dv::Vector{Float64}
    AL2q_bu::Vector{Float64}
    AL2q_cu::Vector{Float64}
    AL2q_bv::Vector{Float64}
    AL2q_cv::Vector{Float64}
end