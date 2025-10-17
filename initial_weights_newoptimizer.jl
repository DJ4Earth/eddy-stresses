# this will attempt to make the weights in my NN match the ZB parameterization term. I will
# then use those weights as an initial guess for further optimization, assuming that it works
using Statistics

mutable struct InitWeightsModel{T, S} <: AbstractNLPModel{T, S}
    meta::NLPModelMeta{T,S}
    counters::Counters
    SZB::ShallowWaters.ModelSetup{T,T}      # model structure, ZB parameterization
    SNN::ShallowWaters.ModelSetup{T,T}      # model struct, NN parameterization
    snapshot::Array{Array{T,2}, 1}          # snapshot being used for offline learning
    J::Float64                              # objective value
    T11::Array{T, 2}
    T22::Array{T, 2}
    T12::Array{T, 2}
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
        nx=128)
    SNN = ShallowWaters.model_setup(PNN)

    ulr = ncread("./spinup_files/128_postspinup_noforcing_cginitcondition_oneyear_071825/u.nc", "u")
    vlr = ncread("./spinup_files/128_postspinup_noforcing_cginitcondition_oneyear_071825/v.nc", "v")
    etalr = ncread("./spinup_files/128_postspinup_noforcing_cginitcondition_oneyear_071825/eta.nc", "eta")

    # param_guess = 1e-1.*randn(Lux.parameterlength(SNN.Diag.CNNVars.model_Su) + Lux.parameterlength(SNN.Diag.CNNVars.model_Sv))
    # param_guess = load_object("./result_workedupto21by21_32-32-32-4CNN_091825.jld2").solution
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

    ucg = load_object("./coarsegrained_hr_ubar_ubarsq_t1_foroffline_101525.jld2")
    vcg = load_object("./coarsegrained_hr_vbar_vbarsq_t1_foroffline_101525.jld2")
    uvbar = load_object("./coarsegrained_hr_uvbar_t1_foroffline_101525.jld2")

    ubar = ucg[1]
    ubarsq = ucg[2]

    vbar = vcg[1]
    vbarsq = vcg[2]

    ubarh = cat(zeros(T,nux+2*halo,halo),cat(zeros(T,halo,nuy),ubar,zeros(T,halo,nuy),dims=1),zeros(T,nux+2*halo,halo),dims=2)
    ubarhsq = cat(zeros(T,nux+2*halo,halo),cat(zeros(T,halo,nuy),ubarsq,zeros(T,halo,nuy),dims=1),zeros(T,nux+2*halo,halo),dims=2)

    vbarh = cat(zeros(T,nvx+2*halo,halo),cat(zeros(T,halo,nvy),vbar,zeros(T,halo,nvy),dims=1),zeros(T,nvx+2*halo,halo),dims=2)
    vbarhsq = cat(zeros(T,nvx+2*halo,halo),cat(zeros(T,halo,nvy),vbarsq,zeros(T,halo,nvy),dims=1),zeros(T,nvx+2*halo,halo),dims=2)

    uvbarh = cat(zeros(T,nx+2*haloη,haloη),cat(zeros(T,haloη,ny),uvbar,zeros(T,haloη,ny),dims=1),zeros(T,nx+2*haloη,haloη),dims=2)

    T12_true = ShallowWaters.Iy(ubarh)[2:end-1,2:end-1] .* ShallowWaters.Ix(vbarh)[2:end-1,2:end-1] - ShallowWaters.Ixy(uvbarh)

    T11_true = (ShallowWaters.Ixy(ShallowWaters.Iy(ubarh)[2:end-1,2:end-1])).^2 - ShallowWaters.Ixy(ShallowWaters.Iy(ubarhsq)[2:end-1,2:end-1])
    T22_true = ShallowWaters.Ixy((ShallowWaters.Ix(vbarh)[2:end-1,2:end-1])).^2 - ShallowWaters.Ixy(ShallowWaters.Ix(vbarhsq)[2:end-1,2:end-1])

    u, v, _ = ShallowWaters.add_halo(ulr[:,:,j], vlr[:,:,j], etalr[:,:,j], SNN)
    snapshot = [u,v]
    meta = NLPModelMeta(Lux.parameterlength(SNN.Diag.CNNVars.model_Su) + Lux.parameterlength(SNN.Diag.CNNVars.model_Sv); ncon=0, nnzh=0,x0=param_guess)
    counters = Counters()

    return InitWeightsModel{T, typeof(param_guess)}(meta, Counters(), SZB, SNN, snapshot, 0.0, T11_true, T22_true, T12_true)

end

function compute_init_weights_newoptimizer()

    nlp = InitWeightsModel{Float64}()
    qn_options = MadNLP.QuasiNewtonOptions(; max_history=200)
    results = madnlp(
        nlp;
        # linear_solver=LapackCPUSolver,
        hessian_approximation=MadNLP.CompactLBFGS,
        quasi_newton_options=qn_options,
        max_iter=500
    )

    return results

end

function for_enzyme(param_guess, state, SNN, SZB, T11, T22, T12)

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

    denom = SZB.grid.Δ^2 * SZB.grid.scale

    # ZB T11 - CNN T11
    J = 0.0
    denom1 = sqrt(sum((SNN.Diag.CNNVars.T11 .- mean(SNN.Diag.CNNVars.T11)).^2))
    J += sum((SNN.Diag.CNNVars.T11 - T11).^2) / 128^2 + ((sqrt(sum((SNN.Diag.CNNVars.T11 .- mean(SNN.Diag.CNNVars.T11)).^2)) - sqrt(sum((T11 .- mean(T11)).^2)))^2) / denom1

    #ZB T22 - CNN T22
    denom2 = sqrt(sum((SNN.Diag.CNNVars.T22 .- mean(SNN.Diag.CNNVars.T22)).^2))
    J += sum((SNN.Diag.CNNVars.T22 - T22).^2) / 128^2 + ((sqrt(sum((SNN.Diag.CNNVars.T22 .- mean(SNN.Diag.CNNVars.T22)).^2)) - sqrt(sum((T22 .- mean(T22)).^2)))^2) / denom2

    #ZB T12 - CNN T12
    denom3 = sqrt(sum((SNN.Diag.CNNVars.T12 .- mean(SNN.Diag.CNNVars.T12)).^2))
    J += sum((SNN.Diag.CNNVars.T12 - T12).^2) / 129^2 + ((sqrt(sum((SNN.Diag.CNNVars.T12 .- mean(SNN.Diag.CNNVars.T12)).^2)) - sqrt(sum((T12 .- mean(T12)).^2)))^2) / denom3


    # ZB T11 - CNN T11
    # J = 0.0
    # J += sum((SNN.Diag.CNNVars.T11 - (SZB.Diag.ZBVars.trace_filtered - SZB.Diag.ZBVars.ζD_filtered)./denom).^2) / (128^2)

    # #ZB T22 - CNN T22
    # J += sum((SNN.Diag.CNNVars.T22 - (SZB.Diag.ZBVars.trace_filtered + SZB.Diag.ZBVars.ζD_filtered)./denom).^2) / (128^2)

    # #ZB T12 - CNN T22
    # J += sum((SNN.Diag.CNNVars.T12 - SZB.Diag.ZBVars.ζDhat_filtered./denom).^2) / (129^2)

    # attempting with a variance regularization term, might need a coefficient here
    # J += sum(SZB.Diag.ZBVars.S_u .- mean(SZB.Diag.ZBVars.S_u).^2) + sum((SZB.Diag.ZBVars.S_v .- mean(SZB.Diag.ZBVars.S_v)).^2)

    return J

end

function NLPModels.obj(model, param_guess)

    PZB = ShallowWaters.Parameter(T = Float64;
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

    ShallowWaters.ZB_momentum(model.snapshot[1], model.snapshot[2], model.SZB, model.SZB.Diag)
    ShallowWaters.CNN_momentum(model.snapshot[1], model.snapshot[2], model.SNN)

    denom = model.SZB.grid.Δ^2 * model.SZB.grid.scale

    ucg = load_object("./coarsegrained_hr_ubar_ubarsq_t1_foroffline_101525.jld2")
    vcg = load_object("./coarsegrained_hr_vbar_vbarsq_t1_foroffline_101525.jld2")

    # trying with the "true" S tensors
    # ZB T11 - CNN T11
    model.J = 0.0
    denom1 = sqrt(sum((model.SNN.Diag.CNNVars.T11 .- mean(model.SNN.Diag.CNNVars.T11)).^2))
    model.J += sum((model.SNN.Diag.CNNVars.T11 - model.T11).^2 ) / 128^2 + ((sqrt(sum((model.SNN.Diag.CNNVars.T11 .- mean(model.SNN.Diag.CNNVars.T11)).^2)) - sqrt(sum((model.T11 .- mean(model.T11)).^2)))^2) / denom1

    #ZB T22 - CNN T22
    denom2 = sqrt(sum((model.SNN.Diag.CNNVars.T22 .- mean(model.SNN.Diag.CNNVars.T22)).^2))
    model.J += sum((model.SNN.Diag.CNNVars.T22 - model.T22).^2) / 128^2  + ((sqrt(sum((model.SNN.Diag.CNNVars.T22 .- mean(model.SNN.Diag.CNNVars.T22)).^2)) - sqrt(sum((model.T22 .- mean(model.T22)).^2)))^2) / denom2

    #ZB T12 - CNN T12
    denom3 = sqrt(sum((model.SNN.Diag.CNNVars.T12 .- mean(model.SNN.Diag.CNNVars.T12)).^2))
    model.J += sum((model.SNN.Diag.CNNVars.T12 - model.T12).^2) / 129^2  + ((sqrt(sum((model.SNN.Diag.CNNVars.T12 .- mean(model.SNN.Diag.CNNVars.T12)).^2)) - sqrt(sum((model.T12 .- mean(model.T12)).^2)))^2) / denom3

    # # ZB T11 - CNN T11
    # model.J += sum((model.SNN.Diag.CNNVars.T11 - (model.SZB.Diag.ZBVars.trace_filtered - model.SZB.Diag.ZBVars.ζD_filtered)./denom).^2) / (128^2)

    # #ZB T22 - CNN T22
    # model.J += sum((model.SNN.Diag.CNNVars.T22 - (model.SZB.Diag.ZBVars.trace_filtered + model.SZB.Diag.ZBVars.ζD_filtered)./denom).^2) / (128^2)

    # #ZB T12 - CNN T12
    # model.J += sum((model.SNN.Diag.CNNVars.T12 - model.SZB.Diag.ZBVars.ζDhat_filtered./denom).^2) / (129^2)

    # attempting with a variance regularization term, might need a coefficient here
    # model.J += sum(model.SZB.Diag.ZBVars.S_u .- mean(model.SZB.Diag.ZBVars.S_u).^2) + sum((model.SZB.Diag.ZBVars.S_v .- mean(model.SZB.Diag.ZBVars.S_v)).^2)

    return model.J

end

function NLPModels.grad!(model, param_guess, G)

    PZB = ShallowWaters.Parameter(T=Float64;
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
        Duplicated(model.SNN, dSNN),
        Const(model.T11),
        Const(model.T22),
        Const(model.T12)
    )[2]

    G .= dparam

    return G

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

function checking()

    nlp = InitWeightsModel{Float64}()

    current = 1
    l = 0
    for model in (nlp.SNN.Diag.CNNVars.model_Su, nlp.SNN.Diag.CNNVars.model_Sv)
        for layers in model[1]
            for array in layers
                    sz = prod(size(array))
                    l += sz
            end
        end
    end
    G = zeros(l)
    NLPModels.grad!(nlp, nlp.meta.x0, G)

    println("norm of initial gradient ", norm(G))

    G = zeros(l)
    NLPModels.grad!(nlp, nlp.meta.x0, G)
    println("should be the same norm ", norm(G))


end