# this will attempt to make the weights in my NN match the ZB parameterization term. I will
# then use those weights as an initial guess for further optimization, assuming that it works
mutable struct InitWeightsModel{T, S} <: AbstractNLPModel{T, S}
    meta::NLPModelMeta{T,S}
    counters::Counters
    SZB::ShallowWaters.ModelSetup{T,T}      # model structure, ZB parameterization
    SNN::ShallowWaters.ModelSetup{T,T}      # model struct, NN parameterization
    snapshot::Array{Array{T,2}, 1}          # snapshot being used for offline learning
    J::Float64                              # objective value
end

function InitWeightsModel{T}() where {T<:AbstractFloat}

    SZB = ShallowWaters.model_setup(output=false,
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
        zb_forcing_dissipation=true,
        zb_filtered=true,
        nn_forcing_momentum=false,
        nn_forcing_dissipation=false,
        N=1,
        α=2,
        nx=128
    )

    SNN = ShallowWaters.model_setup(output=false,
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
        zb_filtered=true,
        nn_forcing_momentum=false,
        nn_forcing_dissipation=true,
        N=1,
        α=2,
        nx=128
    )

    ulr = ncread("./spinup_files/128_postspinup_noforcing_cginitcondition_oneyear_071825/u.nc", "u")
    vlr = ncread("./spinup_files/128_postspinup_noforcing_cginitcondition_oneyear_071825/v.nc", "v")
    etalr = ncread("./spinup_files/128_postspinup_noforcing_cginitcondition_oneyear_071825/eta.nc", "eta")

    param_guess = 1e-1.*randn(Lux.parameterlength(SNN.Diag.CNNVars.model_Su) + Lux.parameterlength(SNN.Diag.CNNVars.model_Sv))
    # param_guess = load_object("./result_workedupto21by21_32-32-32-4CNN_091825.jld2").solution
    # current = 1
    # for model in (SNN.Diag.CNNVars.model_Su, SNN.Diag.CNNVars.model_Sv)
    #     for layers in model[1]
    #         for array in layers
    #                 sz = prod(size(array))
    #                 param_guess[current:(current + sz - 1)] .= vec(array)
    #                 current += sz
    #         end
    #     end
    # end

    result = nothing
    j = 10

    u, v, _ = ShallowWaters.add_halo(ulr[:,:,j], vlr[:,:,j], etalr[:,:,j], SNN)
    snapshot = [u,v]
    meta = NLPModelMeta(Lux.parameterlength(SNN.Diag.CNNVars.model_Su) + Lux.parameterlength(SNN.Diag.CNNVars.model_Sv); ncon=0, nnzh=0,x0=param_guess)
    counters = Counters()

    return InitWeightsModel{T, typeof(param_guess)}(meta, Counters(), SZB, SNN, snapshot, 0.0)

end

function compute_init_weights_newoptimizer()

    nlp = InitWeightsModel{Float64}()
    qn_options = MadNLP.QuasiNewtonOptions(; max_history=200)
    results = madnlp(
        nlp;
        # linear_solver=LapackCPUSolver,
        hessian_approximation=MadNLP.CompactLBFGS,
        quasi_newton_options=qn_options,
    )

    return results

end

function for_enzyme(param_guess, state, SNN, SZB)

    # current = 1
    # for model in (SNN.Diag.NNVars.model_diag, SNN.Diag.NNVars.model_offdiag)
    #     for layers in model[1]
    #         for array in layers
    #             sz = prod(size(array))
    #             array .= reshape(param_guess[current:(current + sz - 1)], size(array)...)
    #             current += sz
    #         end
    #     end
    # end

    current = 1
    for model in (SNN.Diag.CNNVars.model_Su, SNN.Diag.CNNVars.model_Sv)
        for layers in model[1]
            for array in layers
                    sz = prod(size(array))
                    array .= reshape(param_guess[current:(current + sz - 1)], size(array)...)
                    current += sz
            end
        end
    end

    ShallowWaters.ZB_momentum(state[1], state[2], SZB, SZB.Diag)
    ShallowWaters.CNN_momentum(state[1], state[2], SNN)

    # return sum((SZB.Diag.ZBVars.S_u .- SNN.Diag.CNNVars.S_u).^2) ./ (128*127) + sum((SZB.Diag.ZBVars.S_v .- SNN.Diag.CNNVars.S_v).^2) ./ (128*127)
    return sum((SZB.Diag.ZBVars.S_u[25:85,25:85] - SNN.Diag.CNNVars.S_u[25:85,25:85]).^2 + (SZB.Diag.ZBVars.S_v[25:85,25:85] - SNN.Diag.CNNVars.S_v[25:85,25:85]).^2)
    # temp = reshape(collect(1:36), 6, 6)
    # return sum((temp - SNN.Diag.CNNVars.S_u[40:45,40:45]).^2)
end

function NLPModels.obj(model, param_guess)

    PZB = ShallowWaters.parameters(output=false,
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
        zb_forcing_dissipation=true,
        zb_filtered=true,
        nn_forcing_momentum=false,
        nn_forcing_dissipation=false,
        N=1,
        α=2,
        nx=128
    )

    PNN = ShallowWaters.parameters(output=false,
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
        zb_filtered=true,
        nn_forcing_momentum=false,
        nn_forcing_dissipation=true,
        N=1,
        α=2,
        nx=128
    )

    model.SZB = ShallowWaters.model_setup(PZB)
    model.SNN = ShallowWaters.model_setup(PNN)
    model.J = 0

    # current = 1
    # for model in (SNN.Diag.NNVars.model_diag, SNN.Diag.NNVars.model_offdiag)
    #     for layers in model[1]
    #         for array in layers
    #             sz = prod(size(array))
    #             array .= reshape(param_guess[current:(current + sz - 1)], size(array)...)
    #             current += sz
    #         end
    #     end
    # end

    current = 1
    for m in (model.SNN.Diag.CNNVars.model_Su, model.SNN.Diag.CNNVars.model_Sv)
        for layers in m[1]
            for array in layers
                    sz = prod(size(array))
                    array .= reshape(param_guess[current:(current + sz - 1)], size(array)...)
                    current += sz
            end
        end
    end

    ShallowWaters.ZB_momentum(state[1], state[2], SZB, SZB.Diag)
    ShallowWaters.CNN_momentum(state[1], state[2], SNN)

    # return sum((SZB.Diag.ZBVars.S_u .- SNN.Diag.CNNVars.S_u).^2) ./ (128*127) + sum((SZB.Diag.ZBVars.S_v .- SNN.Diag.CNNVars.S_v).^2) ./ (128*127)
    # return sum((SZB.Diag.ZBVars.S_u[45:55,45:55] - SNN.Diag.CNNVars.S_u[45:55,45:55]).^2)
    return sum((SZB.Diag.ZBVars.S_u[25:85,25:85] - SNN.Diag.CNNVars.S_u[25:85,25:85]).^2 + (SZB.Diag.ZBVars.S_v[25:85,25:85] - SNN.Diag.CNNVars.S_v[25:85,25:85]).^2)
    # temp = reshape(collect(1:36), 6, 6)
    # return sum((temp - SNN.Diag.CNNVars.S_u[40:45,40:45]).^2)

    # return SZB, SNN

end

function NLPModels.grad!(model, param_guess, G)

    PZB = ShallowWaters.parameters(output=false,
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
        zb_forcing_dissipation=true,
        zb_filtered=true,
        nn_forcing_momentum=false,
        nn_forcing_dissipation=false,
        N=1,
        α=2,
        nx=128
    )

    PNN = ShallowWaters.parameters(output=false,
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
        zb_filtered=true,
        nn_forcing_momentum=false,
        nn_forcing_dissipation=true,
        N=1,
        α=2,
        nx=128
    )

    model.SZB = ShallowWaters.model_setup(PZB)
    model.SNN = ShallowWaters.model_setup(PNN)
    model.J = 0

    dparam = Enzyme.make_zero(param_guess)
    dSZB = Enzyme.make_zero(model.SZB)
    dSNN = Enzyme.make_zero(model.SNN)

    J = autodiff(
        set_runtime_activity(Enzyme.ReverseWithPrimal),
        for_enzyme,
        Active,
        Duplicated(param_guess, dparam),
        Const(model.snapshot),
        Duplicated(model.SZB, dSZB),
        Duplicated(model.SNN, dSNN)
    )[2]

    G .= dparam

    return G

end

function ignore(result)

    ulr = ncread("./spinup_files/128_postspinup_noforcing_cginitcondition_oneyear_071825/u.nc", "u")
    vlr = ncread("./spinup_files/128_postspinup_noforcing_cginitcondition_oneyear_071825/v.nc", "v")
    etalr = ncread("./spinup_files/128_postspinup_noforcing_cginitcondition_oneyear_071825/eta.nc", "eta")

    param_guess = result.solution
    j = 10

    SZB = ShallowWaters.model_setup(output=false,
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
        zb_forcing_dissipation=true,
        zb_filtered=true,
        nn_forcing_momentum=false,
        nn_forcing_dissipation=false,
        N=1,
        α=2,
        nx=128
    )

    u, v, _ = ShallowWaters.add_halo(ulr[:,:,j], vlr[:,:,j], etalr[:,:,j], SZB)
    state = [u,v]

    SNN = ShallowWaters.model_setup(output=false,
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
        zb_filtered=true,
        nn_forcing_momentum=false,
        nn_forcing_dissipation=true,
        N=1,
        α=2,
        nx=128
    )

    # current = 1
    # for model in (SNN.Diag.NNVars.model_diag, SNN.Diag.NNVars.model_offdiag)
    #     for layers in model[1]
    #         for array in layers
    #             sz = prod(size(array))
    #             array .= reshape(param_guess[current:(current + sz - 1)], size(array)...)
    #             current += sz
    #         end
    #     end
    # end

    current = 1
    for model in (SNN.Diag.CNNVars.model_Su, SNN.Diag.CNNVars.model_Sv)
        for layers in model[1]
            for array in layers
                    sz = prod(size(array))
                    array .= reshape(param_guess[current:(current + sz - 1)], size(array)...)
                    current += sz
            end
        end
    end

    ShallowWaters.ZB_momentum(state[1], state[2], SZB, SZB.Diag)
    ShallowWaters.CNN_momentum(state[1], state[2], SNN)

    return SZB, SNN


end

function error()

    Slr = ShallowWaters.model_setup(output=false,
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
        zb_filtered=true,
        nn_forcing_momentum=false,
        nn_forcing_dissipation=true,
        N=1,
        α=2,
        nx=128
    )

    param_guess = 1000 .* randn(883)

    dparam = Enzyme.make_zero(param_guess)
    # dstate = Enzyme.make_zero(state)

    @time autodiff(
        set_runtime_activity(Enzyme.Reverse),
        initweights_compute_loss,
        Active,
        Duplicated(param_guess, dparam),
        Const([Slr.Prog.u, Slr.Prog.v])
    )



end