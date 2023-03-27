# returning an out of bounds error from something in checkpointing

using Parameters
using Enzyme, Checkpointing, Zygote

@with_kw mutable struct checkpoint_struct_debug

    Nu::Int
    Nv::Int
    NT::Int
    Nq::Int

    u::Vector{Float64} = zeros(Nu)
    v::Vector{Float64} = zeros(Nv)
    eta::Vector{Float64} = zeros(NT)
    u0::Vector{Float64} = zeros(Nu) 
    v0::Vector{Float64} = zeros(Nv)
    eta0::Vector{Float64} = zeros(NT)

end

function advance_debug(states_rhs)

    @assert all(x -> x < 7.0, states_rhs.u0)
    @assert all(x -> x < 7.0, states_rhs.v0)
    @assert all(x -> x < 7.0, states_rhs.eta0)

    return nothing 

end

function chkpt_debug(chkpt_struct::checkpoint_struct_debug, T::Int)
    @checkpoint_struct revolve chkpt_struct for j in 1:T
        advance_debug(chkpt_struct)
        copyto!(chkpt_struct.u, chkpt_struct.u0)
        copyto!(chkpt_struct.v, chkpt_struct.v0)
        copyto!(chkpt_struct.eta, chkpt_struct.eta0)
    end
    return chkpt_struct.u
end

nx = 20
ny = 20 

deriv_struct = checkpoint_struct_debug(Nu = (nx - 1) * ny, 
    Nv = (ny - 1) * nx, 
    NT = nx * ny, 
    Nq = (nx + 1) * (ny + 1)
)

Trun = 100
snaps = 30
verbose = 0
revolve = Revolve{checkpoint_struct_debug}(Trun, snaps; verbose=verbose)

g = Zygote.jacobian(chkpt_debug, deriv_struct, Trun)

# checking that the function advance_debug runs
# for j in 1:10
#     advance_debug(deriv_struct)#, grid, gyre_params, interp_ops, grad_ops, advec_ops)
# end

