"""
Placing the checkpointed integraton, loop, and all of the packages here, so that when running 
include("technical_paper.jl") we never re-include them. This should help with Enzyme compile times

When running the experiments for the technical paper only ever include("technical_paper_integration.jl") once
"""

include("../ShallowWaters.jl/src/ShallowWaters.jl")
using .ShallowWaters

using Enzyme
using Checkpointing

mutable struct MyPrognosticVars{T<:AbstractFloat}
    u::Array{T,2}           # u-velocity
    nu::Array{T,2}           # sea surface height / interface displacement
end

function loop(S,scheme)
    @checkpoint_struct scheme S for S.parameters.i = 1:S.grid.nt
        Diag = S.Diag
        Prog = S.Prog

        halo = S.grid.halo

        # undo scaling as well
        @views ucut = S.Prog.u[halo+1:end-halo,halo+1:end-halo]
        ηcut = S.Prog.η

        temp = MyPrognosticVars{Float64}(ucut, ηcut) 

        energy_lr = first(temp.u.^2)

        S.parameters.J += energy_lr
    end
    return nothing
end

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

