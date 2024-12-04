"""
Placing the checkpointed integraton, loop, and all of the packages here, so that when running 
include("technical_paper.jl") we never re-include them. This should help with Enzyme compile times

When running the experiments for the technical paper only ever include("technical_paper_integration.jl") once
"""

include("../ShallowWaters.jl/src/ShallowWaters.jl")
using .ShallowWaters

using Enzyme
using Checkpointing, HDF5, Serialization
using NetCDF, JLD2, CairoMakie

Enzyme.API.looseTypeAnalysis!(true)

using Parameters
using Optim
using LaTeXStrings

function checkpointed_integration(S, scheme)
    h = S.Diag.VolumeFluxes.h
    H = S.forcing.H
    nu = S.Prog.η
    @inbounds for i in eachindex(nu)
        h[i] = nu[i]
    end
    # run integration loop with checkpointing
    loop(S, scheme)

    return S.parameters.J

end

function myremove_halo(   u::Array{T,2},
                        v::Array{T,2},
                        η::Array{T,2},
                        sst::Array{T,2},
                        S) where {T<:AbstractFloat}

    @unpack halo,haloη,halosstx,halossty = S.grid

    # undo scaling as well
    @views ucut = u[halo+1:end-halo,halo+1:end-halo]
    ηcut = η

    return ucut,ηcut
end

mutable struct MyPrognosticVars{T<:AbstractFloat}
    u::Array{T,2}           # u-velocity
    nu::Array{T,2}           # sea surface height / interface displacement
end

function loop(S,scheme)


    @checkpoint_struct scheme S for S.parameters.i = 1:S.grid.nt

        Diag = S.Diag
        Prog = S.Prog

        # if S.parameters.i in (S.grid.nt - 30*224):1:S.grid.nt

            temp = MyPrognosticVars{Float64}(myremove_halo(S.Prog.u,
            S.Prog.v,
            S.Prog.η,
            S.Prog.sst,S)...)

            energy_lr = first(temp.u.^2)

            S.parameters.J += energy_lr

        # end

    end

    ##### use if time-averaging the objective function #######
    ##########################################################

    ##### Energy objective function, not time averaged ###########
    # temp = ShallowWaters.PrognosticVars{Float32}(ShallowWaters.remove_halo(S.Prog.u,
    # S.Prog.v,
    # S.Prog.η,
    # S.Prog.sst,S)...)

    # energy_lr = (sum(temp.u.^2) + sum(temp.v.^2)) / (S.grid.nx * S.grid.ny)
    # S.parameters.J = energy_lr
    ###########################################

    return nothing

end
