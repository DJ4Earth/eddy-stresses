# This script will contain the structures needed for a barotropic gyre model

mutable struct gyre

    lastu::Matrix{Float64}
    lastv::Matrix{Float64}
    lastη::Matrix{Float64} 
    
    u::Matrix{Float64}
    v::Matrix{Float64}
    η::Matrix{Float64}

end

mutable struct gyre_Parameters 
    g::Float64              # gravity
    f0::Float64             # Coriolis parameter
    beta::Float64           # Coriolis parameter
    τ_0::Float64             # wind stress amplitude  
    H::Float64              # ocean depth 
    A_h::Float64            # horizontal Laplacian viscosity 
    ρ_c::Float64            # reference density 
    eps_ab::Float64           # parameter pertaining to AB2
    Nx::Int 
    Ny::Int 
    dx::Float64 
    dy::Float64 
    dt::Float64 
    τ::Function              # wind stress function 
    f::Function              # coriolis beta-plane approx 
    wind_stress::Matrix{Float64}    # wind stress values on the grid 
    coriolis::Matrix{Float64}       # coriolis values on the grid 
end
