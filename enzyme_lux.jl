using Enzyme
# Enzyme.Compiler.VERBOSE_ERRORS[] = true
using Checkpointing, HDF5, Serialization
using NetCDF, JLD2, CairoMakie
using Lux, Random

# Enzyme.API.looseTypeAnalysis!(true)

using Parameters
using Optim
using LaTeXStrings

include("../ShallowWaters.jl/src/ShallowWaters.jl")
using .ShallowWaters

mutable struct Chkp{T1,T2}
    S::ShallowWaters.ModelSetup{T1,T2}      # model structure
    data_steps::StepRange{Int, Int}         # location of data points temporally
    J::Float64                              # objective function value
    input1::Array{T1,1}
    input2::Array{T1,1}
    j::Int                                  # for keeping track of location in data
    i::Int                                  # timestep iterator
    t::Int64                                # model time
end

function luxfunction(chkp, scheme)

    input1 = chkp.input1
    input2 = chkp.input2

    weights_corner = chkp.S.Diag.NNVars.weights_corner
    weights_center = chkp.S.Diag.NNVars.weights_center

    corner_layers = Lux.Dense(22 => 2, relu)
    center_layers = Lux.Dense(17 => 1, relu)

    ps_corner, st_corner = Lux.setup(Random.default_rng(), corner_layers)
    ps_center, st_center = Lux.setup(Random.default_rng(), center_layers)

    corner_params = (weight=weights_corner, bias=ps_corner.bias)
    center_params = (weight=weights_center, bias=ps_center.bias)

    @checkpoint_struct scheme chkp for chkp.i = 1:300

        corner_model = StatefulLuxLayer{true}(corner_layers, corner_params, st_corner)
        center_model = StatefulLuxLayer{true}(center_layers, center_params, st_center)

        _ = corner_model(input1)
        b1 = center_model(input2)

        chkp.J += b1[1]

    end

    return chkp.J

end

function compute_gradient()

    # Type precision
    T = Float32

    P = ShallowWaters.Parameter(T=T;
        output=false,
        L_ratio=1,
        g=9.81,
        H=500,
        wind_forcing_x="double_gyre",
        Lx=3840e3,
        seasonal_wind_x=false,
        topography="flat",
        bc="nonperiodic",
        bottom_drag="quadratic",
        tracer_advection=false,
        tracer_relaxation=false,
        zb_forcing_momentum=false,
        zb_forcing_dissipation=false,
        nn_forcing_momentum=false,
        nn_forcing_dissipation=false,
        handwritten=false,
        zb_filtered=true,
        N=1,
        α=2,
        nx=128,
        Ndays=1,
        initial_cond="rest",
    )

    S = ShallowWaters.model_setup(P)

    input1 = randn(T,22)
    input2 = randn(T,17)
    weights_corner = randn(T,2,22)
    weights_center = randn(T,1,17)

    S.Diag.NNVars.weights_corner = weights_corner
    S.Diag.NNVars.weights_center = weights_center

    data_steps = 1:1:S.grid.nt

    snaps = Int(floor(sqrt(S.grid.nt)))
    revolve = Revolve{Chkp{T, T}}(S.grid.nt,
        snaps;
        verbose=1,
        gc=true,
        write_checkpoints=false,
        write_checkpoints_filename = "",
        write_checkpoints_period = 224
    )

    chkp = Chkp{T, T}(S,
        data_steps,
        0.0,
        input1,
        input2,
        1,
        1,
        0
    )
    dchkp = Enzyme.make_zero(chkp)

    snaps = 10
    revolve = Revolve{Chkp}(1,
        snaps;
        verbose=1,
        gc=true,
        write_checkpoints=false,
        write_checkpoints_filename = "",
        write_checkpoints_period = 50
    )

    J = autodiff(
        set_runtime_activity(Enzyme.ReverseWithPrimal),
        luxfunction,
        Active,
        Duplicated(chkp, make_zero(chkp)),
        Const(revolve)
    )[2]
    println("Cost with AD: $J")

    return nothing

end

# G = compute_gradient()