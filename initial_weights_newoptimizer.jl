# this will attempt to make the weights in my NN match the ZB parameterization term. I will
# then use those weights as an initial guess for further optimization, assuming that it works
using Statistics

mutable struct InitWeightsModel{T, S} <: AbstractNLPModel{T, S}
    meta::NLPModelMeta{T,S}
    counters::Counters
    SZB::ShallowWaters.ModelSetup{T,T}      # model structure, ZB parameterization
    SNN::ShallowWaters.ModelSetup{T,T}      # model struct, NN parameterization
    snapshot::Array{Array{T,3}, 1}          # snapshot(s) being used for offline learning
    J::Float64                              # objective value
    T11::Array{T, 3}                        # for storing T11 computed from *coarse-grained and filtered* hr states
    T22::Array{T, 3}                        # same as above
    T12::Array{T, 3}
end

function InitWeightsModel{T}() where {T<:AbstractFloat}

    PZB = ShallowWaters.Parameter(T=T;
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
        zb_forcing_dissipation=true,
        zb_filtered=true,
        nn_forcing_momentum=false,
        nn_forcing_dissipation=false,
        N=1,
        α=2,
        nx=128
    )
    SZB = ShallowWaters.model_setup(PZB)

    PNN = ShallowWaters.Parameter(T=T;
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
        zb_filtered=true,
        nn_forcing_momentum=false,
        nn_forcing_dissipation=true,
        N=1,
        α=2,
        nx=1024)
    SNN = ShallowWaters.model_setup(PNN)

    # param_guess = 1e-1.*randn(Lux.parameterlength(SNN.Diag.CNNVars.model_Su) + Lux.parameterlength(SNN.Diag.CNNVars.model_Sv))
    # param_guess = load_object("./offline_files/offlineresult_3-25-25_1e-3objective_relu_activation.jld2").solution
    param_guess = zeros(Lux.parameterlength(SNN.Diag.CNNVars.model_Su) + Lux.parameterlength(SNN.Diag.CNNVars.model_Sv))
    current = 1
    for model in (SNN.Diag.CNNVars.model_Su, SNN.Diag.CNNVars.model_Sv)
        for layers in model[1]
            for array in layers
                    sz = prod(size(array))
                    param_guess[current:(current + sz - 1)] .= vec(array)
                    current += sz
            end
        end
    end

    j = 1

    nx = SZB.grid.nx
    ny = SZB.grid.ny
    nux = SZB.grid.nux
    nuy = SZB.grid.nuy
    nvx = SZB.grid.nvx
    nvy = SZB.grid.nvy
    halo = SZB.grid.halo
    haloη = SZB.grid.haloη

    filteredstates = load_object("./offline_files/hrstates_filtered_downsized_hourly_tendays_uveta_beginsatonehour_102825.jld2")
    true_Ts = load_object("./offline_files/true_Ts_hourlysaves_filteredandcg_T11T22T12_beginsatonehour_111725.jld2")

    # lvar is by default -Inf * ones(Float64, nvar)
    # uvar is by default Inf * ones(Float64, nvar)

    # u, v, _ = ShallowWaters.add_halo(ulr[:,:,j], vlr[:,:,j], etalr[:,:,j], SNN)

    meta = NLPModelMeta(Lux.parameterlength(SNN.Diag.CNNVars.model_Su) + Lux.parameterlength(SNN.Diag.CNNVars.model_Sv);
        ncon=0,
        nnzh=0,
        x0=param_guess
    )

    counters = Counters()

    return InitWeightsModel{T, typeof(param_guess)}(meta, Counters(), SZB, SNN, filteredstates, 0.0, true_Ts[1], true_Ts[2], true_Ts[3])

end

function compute_init_weights_newoptimizer()

    nlp = InitWeightsModel{Float64}()
    qn_options = MadNLP.QuasiNewtonOptions(; max_history=200)
    results = madnlp(
        nlp;
        # linear_solver=LapackCPUSolver,
        hessian_approximation=MadNLP.CompactLBFGS,
        quasi_newton_options=qn_options,
        max_iter=150
    )

    return results

end

# function for_enzyme(param_guess, state, SNN, T11, T22, T12)

#     current = 1
#     for model in (SNN.Diag.CNNVars.model_Su, SNN.Diag.CNNVars.model_Sv)
#         for layers in model[1]
#             for array in layers
#                     sz = prod(size(array))
#                     array .= reshape(param_guess[current:(current + sz - 1)], size(array)...)
#                     current += sz
#             end
#         end
#     end

#     ShallowWaters.ZB_momentum(state[1], state[2], SZB, SZB.Diag)
#     ShallowWaters.CNN_momentum(state[1], state[2], SNN)

#     denom = SZB.grid.Δ^2 * SZB.grid.scale

#     # True T11 - CNN T11
#     J = 0.0
#     denom1 = sqrt(sum((SNN.Diag.CNNVars.T11 .- mean(SNN.Diag.CNNVars.T11)).^2))
#     J += sum((SNN.Diag.CNNVars.T11 - T11).^2) / 128^2 + ((sqrt(sum((SNN.Diag.CNNVars.T11 .- mean(SNN.Diag.CNNVars.T11)).^2)) - sqrt(sum((T11 .- mean(T11)).^2)))^2) / denom1

#     # True T22 - CNN T22
#     denom2 = sqrt(sum((SNN.Diag.CNNVars.T22 .- mean(SNN.Diag.CNNVars.T22)).^2))
#     J += sum((SNN.Diag.CNNVars.T22 - T22).^2) / 128^2 + ((sqrt(sum((SNN.Diag.CNNVars.T22 .- mean(SNN.Diag.CNNVars.T22)).^2)) - sqrt(sum((T22 .- mean(T22)).^2)))^2) / denom2

#     # True T12 - CNN T12
#     denom3 = sqrt(sum((SNN.Diag.CNNVars.T12 .- mean(SNN.Diag.CNNVars.T12)).^2))
#     J += sum((SNN.Diag.CNNVars.T12 - T12).^2) / 129^2 + ((sqrt(sum((SNN.Diag.CNNVars.T12 .- mean(SNN.Diag.CNNVars.T12)).^2)) - sqrt(sum((T12 .- mean(T12)).^2)))^2) / denom3

#     return J

# end

function for_enzyme(model, param_guess)

    PNN = ShallowWaters.Parameter(T = Float64;
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
        zb_filtered=true,
        nn_forcing_momentum=false,
        nn_forcing_dissipation=true,
        N=1,
        α=2,
        nx=128
    )

    model.SNN = ShallowWaters.model_setup(PNN)

    J = 0
    # adding up the difference for 5 different hourly snapshots of the coarse-grained, hr states true Ts and 
    # the output from the NN
    for j = 1:5

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

        u,v,_ = ShallowWaters.add_halo(model.snapshot[1][:,:,j], model.snapshot[2][:,:,j], model.snapshot[3][:,:,1],model.SNN)
        ShallowWaters.CNN_momentum(u, v, model.SNN)
        T11 = model.T11[:,:,j]
        T22 = model.T22[:,:,j]
        T12 = model.T12[:,:,j]

        denom = model.SNN.grid.Δ^2 * model.SNN.grid.scale

        # trying with the "true" S tensors
        # True T11 - CNN T11
        denom = model.SZB.grid.Δ^2 * model.SZB.grid.scale

        # trying with the "true" S tensors
        # True T11 - CNN T11
        denom1 = sqrt(sum((model.SNN.Diag.CNNVars.T11 .- mean(model.SNN.Diag.CNNVars.T11)).^2))
        model.J += sum((model.SNN.Diag.CNNVars.T11 - T11).^2) / 128^2 + ((sqrt(sum((model.SNN.Diag.CNNVars.T11 .- mean(model.SNN.Diag.CNNVars.T11)).^2)) - sqrt(sum((T11 .- mean(T11)).^2)))^2) / denom1

        # True T22 - CNN T22
        denom2 = sqrt(sum((model.SNN.Diag.CNNVars.T22 .- mean(model.SNN.Diag.CNNVars.T22)).^2))
        model.J += sum((model.SNN.Diag.CNNVars.T22 - T22).^2) / 128^2 + ((sqrt(sum((model.SNN.Diag.CNNVars.T22 .- mean(model.SNN.Diag.CNNVars.T22)).^2)) - sqrt(sum((T22 .- mean(T22)).^2)))^2) / denom2

        # True T12 - CNN T12
        denom3 = sqrt(sum((model.SNN.Diag.CNNVars.T12 .- mean(model.SNN.Diag.CNNVars.T12)).^2))
        model.J += sum((model.SNN.Diag.CNNVars.T12 - T12).^2) / 129^2 + ((sqrt(sum((model.SNN.Diag.CNNVars.T12 .- mean(model.SNN.Diag.CNNVars.T12)).^2)) - sqrt(sum((T12 .- mean(T12)).^2)))^2) / denom3

    end

    return model.J

end

function NLPModels.obj(model, param_guess)

    PNN = ShallowWaters.Parameter(T = Float64;
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
        zb_filtered=true,
        nn_forcing_momentum=false,
        nn_forcing_dissipation=true,
        N=1,
        α=2,
        nx=128
    )

    model.J = 0

    # adding up the difference for 5 different snapshots of the coarse-grained, hr states true Ts and 
    # the output from the NN
    for j = 1:5

        model.SNN = ShallowWaters.model_setup(PNN)

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

        u,v,_ = ShallowWaters.add_halo(model.snapshot[1][:,:,j], model.snapshot[2][:,:,j], model.snapshot[3][:,:,1],model.SNN)

        ShallowWaters.CNN_momentum(u, v, model.SNN)
        T11 = model.T11[:,:,j]
        T22 = model.T22[:,:,j]
        T12 = model.T12[:,:,j]

        denom = model.SNN.grid.Δ^2 * model.SNN.grid.scale

        # trying with the "true" S tensors
        # True T11 - CNN T11
        denom1 = sqrt(sum((model.SNN.Diag.CNNVars.T11 .- mean(model.SNN.Diag.CNNVars.T11)).^2))
        model.J += sum((model.SNN.Diag.CNNVars.T11 - T11).^2) / 128^2 + ((sqrt(sum((model.SNN.Diag.CNNVars.T11 .- mean(model.SNN.Diag.CNNVars.T11)).^2)) - sqrt(sum((T11 .- mean(T11)).^2)))^2) / denom1

        # True T22 - CNN T22
        denom2 = sqrt(sum((model.SNN.Diag.CNNVars.T22 .- mean(model.SNN.Diag.CNNVars.T22)).^2))
        model.J += sum((model.SNN.Diag.CNNVars.T22 - T22).^2) / 128^2 + ((sqrt(sum((model.SNN.Diag.CNNVars.T22 .- mean(model.SNN.Diag.CNNVars.T22)).^2)) - sqrt(sum((T22 .- mean(T22)).^2)))^2) / denom2

        # True T12 - CNN T12
        denom3 = sqrt(sum((model.SNN.Diag.CNNVars.T12 .- mean(model.SNN.Diag.CNNVars.T12)).^2))
        model.J += sum((model.SNN.Diag.CNNVars.T12 - T12).^2) / 129^2 + ((sqrt(sum((model.SNN.Diag.CNNVars.T12 .- mean(model.SNN.Diag.CNNVars.T12)).^2)) - sqrt(sum((T12 .- mean(T12)).^2)))^2) / denom3

    end

    return model.J

end

function NLPModels.grad!(model, param_guess, G)

    PNN = ShallowWaters.Parameter(T = Float64;
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
        zb_filtered=true,
        nn_forcing_momentum=false,
        nn_forcing_dissipation=true,
        N=1,
        α=2,
        nx=128
    )

    model.SNN = ShallowWaters.model_setup(PNN)
    model.J = 0

    dparam = Enzyme.make_zero(param_guess)
    dmodel = Enzyme.make_zero(model)

    J = autodiff(
        set_runtime_activity(Enzyme.ReverseWithPrimal),
        for_enzyme,
        Active,
        Duplicated(model, dmodel),
        Duplicated(param_guess, dparam),
    )[2]

    G .= dparam

    return nothing

end

function ignore(result)

    ulr = ncread("./spinup_files/128_postspinup_noforcing_cginitcondition_oneyear_071825/u.nc", "u")
    vlr = ncread("./spinup_files/128_postspinup_noforcing_cginitcondition_oneyear_071825/v.nc", "v")
    etalr = ncread("./spinup_files/128_postspinup_noforcing_cginitcondition_oneyear_071825/eta.nc", "eta")

    param_guess = result.solution
    j = 1

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

# function error()

#     Slr = ShallowWaters.model_setup(output=false,
#         L_ratio=1,
#         g=9.81,
#         H=500,
#         wind_forcing_x="double_gyre",
#         Lx=3840e3,
#         seasonal_wind_x=false,
#         topography="flat",
#         bc="nonperiodic",
#         bottom_drag="quadratic",
#         tracer_advection=false,
#         tracer_relaxation=false,
#         zb_forcing_momentum=false,
#         zb_forcing_dissipation=false,
#         zb_filtered=true,
#         nn_forcing_momentum=false,
#         nn_forcing_dissipation=true,
#         N=1,
#         α=2,
#         nx=128
#     )

#     param_guess = 1000 .* randn(883)

#     dparam = Enzyme.make_zero(param_guess)
#     # dstate = Enzyme.make_zero(state)

#     @time autodiff(
#         set_runtime_activity(Enzyme.Reverse),
#         initweights_compute_loss,
#         Active,
#         Duplicated(param_guess, dparam),
#         Const([Slr.Prog.u, Slr.Prog.v])
#     )

# end

# function checking()

#     nlp = InitWeightsModel{Float64}()

#     current = 1
#     l = 0
#     for model in (nlp.SNN.Diag.CNNVars.model_Su, nlp.SNN.Diag.CNNVars.model_Sv)
#         for layers in model[1]
#             for array in layers
#                     sz = prod(size(array))
#                     l += sz
#             end
#         end
#     end
#     G = zeros(l)
#     NLPModels.grad!(nlp, nlp.meta.x0, G)

#     println("norm of initial gradient ", norm(G))

#     G = zeros(l)
#     NLPModels.grad!(nlp, nlp.meta.x0, G)
#     println("should be the same norm ", norm(G))


# end