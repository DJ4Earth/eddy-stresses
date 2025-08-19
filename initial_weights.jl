# this will attempt to make the weights in my NN match the ZB parameterization term. I will
# then use those weights as an initial guess for further optimization, assuming that it works

function for_enzyme(param_guess, state, SNN, SZB)

    current = 1
    for model in (SNN.Diag.NNVars.model_diag, SNN.Diag.NNVars.model_offdiag)
        for layers in model[1]
            for array in layers
                sz = prod(size(array))
                array .= reshape(param_guess[current:(current + sz - 1)], size(array)...)
                current += sz
            end
        end
    end

    ShallowWaters.ZB_momentum(state[1], state[2], SZB, SZB.Diag)
    ShallowWaters.NN_momentum(state[1], state[2], SNN)

    # return sum((SZB.Diag.ZBVars.S_u .- SNN.Diag.NNVars.S_u).^2) ./ (128*127) #+ sum((SZB.Diag.ZBVars.S_v .- SNN.Diag.NNVars.S_v).^2) ./ (128*127)
    # return sum((SZB.Diag.ZBVars.S_u[50:52,50] - SNN.Diag.NNVars.S_u[50:52,50]).^2)
    temp = reshape(collect(1:36), 6, 6)
    return sum((temp - SNN.Diag.NNVars.T11[40:45,40:45]).^2)
end

function initweights_compute_loss(param_guess, state)

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

    current = 1
    for model in (SNN.Diag.NNVars.model_diag, SNN.Diag.NNVars.model_offdiag)
        for layers in model[1]
            for array in layers
                sz = prod(size(array))
                array .= reshape(param_guess[current:(current + sz - 1)], size(array)...)
                current += sz
            end
        end
    end

    ShallowWaters.ZB_momentum(state[1], state[2], SZB, SZB.Diag)
    ShallowWaters.NN_momentum(state[1], state[2], SNN)

    # return sum((SZB.Diag.ZBVars.S_u .- SNN.Diag.NNVars.S_u).^2) ./ (128*127) #+ sum((SZB.Diag.ZBVars.S_v .- SNN.Diag.NNVars.S_v).^2) ./ (128*127)
    # return sum((SZB.Diag.ZBVars.S_u[50:52,50] - SNN.Diag.NNVars.S_u[50:52,50]).^2)
    temp = reshape(collect(1:36), 6, 6)
    return sum((temp - SNN.Diag.NNVars.T11[40:45,40:45]).^2)

end

function initweights_compute_gradient(G, param_guess, state)

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

    dparam = Enzyme.make_zero(param_guess)
    dSZB = Enzyme.make_zero(SZB)
    dSNN = Enzyme.make_zero(SNN)

    J = autodiff(
        set_runtime_activity(Enzyme.ReverseWithPrimal),
        for_enzyme,
        Active,
        Duplicated(param_guess, dparam),
        Const(state),
        Duplicated(SZB, dSZB),
        Duplicated(SNN, dSNN)
    )[2]

    G .= dparam

    return nothing

end

function initweights_FG(F, G, param_guess, state)

    G === nothing || initweights_compute_gradient(G, param_guess, state)
    F === nothing || return initweights_compute_loss(param_guess, state)

end

function compute_init_weights()

    S = ShallowWaters.model_setup(output=false,
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

    ulr = ncread("./spinup_files/128_postspinup_noforcing_cginitcondition_oneyear_071825/u.nc", "u")
    vlr = ncread("./spinup_files/128_postspinup_noforcing_cginitcondition_oneyear_071825/v.nc", "v")
    etalr = ncread("./spinup_files/128_postspinup_noforcing_cginitcondition_oneyear_071825/eta.nc", "eta")

    param_guess = randn(5442+5521)

    result = nothing
    for j in [10]

        u, v, _ = ShallowWaters.add_halo(ulr[:,:,j], vlr[:,:,j], etalr[:,:,j], S)
        fg!_closure(F, G, param_guess) = initweights_FG(F, G, param_guess, [u, v])
        obj_fg = Optim.only_fg!(fg!_closure)
        result = Optim.optimize(obj_fg, param_guess, Optim.Adam(; alpha = .1), Optim.Options(show_trace=true, store_trace=true, iterations=1500))
        param_guess = result.minimizer

    end

    return result

end

function ignore(result)

    ulr = ncread("./spinup_files/128_postspinup_noforcing_cginitcondition_oneyear_071825/u.nc", "u")
    vlr = ncread("./spinup_files/128_postspinup_noforcing_cginitcondition_oneyear_071825/v.nc", "v")
    etalr = ncread("./spinup_files/128_postspinup_noforcing_cginitcondition_oneyear_071825/eta.nc", "eta")

    param_guess = result.minimizer
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

    current = 1
    for model in (SNN.Diag.NNVars.model_diag, SNN.Diag.NNVars.model_offdiag)
        for layers in model[1]
            for array in layers
                sz = prod(size(array))
                array .= reshape(param_guess[current:(current + sz - 1)], size(array)...)
                current += sz
            end
        end
    end

    ShallowWaters.ZB_momentum(state[1], state[2], SZB, SZB.Diag)
    ShallowWaters.NN_momentum(state[1], state[2], SNN)

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