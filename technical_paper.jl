using Enzyme
using Checkpointing

mutable struct MyPrognosticVars{T<:AbstractFloat}
    u::Array{T,2}           # u-velocity
    nu::Array{T,2}           # sea surface height / interface displacement
end

mutable struct MyModelSetup
    i::Int
    J::Float64
    halo::Int
    nt::Int
    Prog::MyPrognosticVars{Float32}
    h::Matrix{Float64}
    t::Int                              # SW: I believe this has something to do with Checkpointing, need to verify
end

function checkpointed_integration(S, scheme)
    h = S.h
    nu = S.Prog.nu
      @inbounds  h[1] = 2.7
    
    @checkpoint_struct scheme S for S.i = 1:S.nt
        Prog = S.Prog

        halo = S.halo

        # undo scaling as well
        @views ucut = S.Prog.u[halo+1:end-halo,halo+1:end-halo]
        ηcut = S.Prog.nu

        temp = MyPrognosticVars{Float64}(ucut, ηcut) 

        energy_lr = first(temp.u.^2)

        S.J += energy_lr
    end

    return
end

function mymodel_setup()
    nt = 6733
    S = MyModelSetup(0, 0.0, 2, nt,MyPrognosticVars{Float32}(ones(Float32, 131, 132),ones(Float32, 130, 130)),ones(Float64, 130, 130),0)

    return S

end
function run_adjoint_plusfd()

    S = mymodel_setup()


    dS = Enzyme.make_zero(S)
    snaps = Int(floor(sqrt(S.nt)))
    revolve = Revolve{MyModelSetup}(S.nt,
        snaps;
        gc=true,
    )

    autodiff(Enzyme.Reverse, checkpointed_integration, Duplicated(S, dS), Const(revolve))

end

diffs, enzyme_deriv, S, dS = run_adjoint_plusfd()
