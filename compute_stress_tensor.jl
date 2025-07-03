# for computing the stress tensor from high-resolution model states

function compute_stress_tensor()

    Ndays = 1
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
        handwritten=false,
        N=1,
        α=2,
        nx=128,
        Ndays=Ndays
    )

    Shr = ShallowWaters.model_setup(output=false,
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
        N=1,
        α=2,
        nx=1024,
        Ndays=Ndays
    )

    nuxlr, nuylr = size(Slr.Prog.u)
    nuxhr, nuyhr = size(Shr.Prog.u)

    nvxlr, nvylr = size(Slr.Prog.v)
    nvxhr, nvyhr = size(Shr.Prog.v)

    dudxhr = Shr.Diag.NNVars.dudx
    dvdyhr = Shr.Diag.NNVars.dvdy
    Suhr = Shr.Diag.NNVars.S_u
    Svhr = Shr.Diag.NNVars.S_v

    dudxlr = Slr.Diag.NNVars.dudx
    dvdylr = Slr.Diag.NNVars.dvdy
    Sulr = Slr.Diag.NNVars.S_u
    Svlr = Slr.Diag.NNVars.S_v

    uhrall = ncread("./spinup_files/1024_postspinup_noslip_5years_061824/u.nc", "u")
    vhrall = ncread("./spinup_files/1024_postspinup_noslip_5years_061824/v.nc", "v")
    etahrall = ncread("./spinup_files/1024_postspinup_noslip_5years_061824/eta.nc", "eta")

    S_hr = []

    ucg = zeros(127,128)
    vcg = zeros(128,127)

    for j = 1:10

        # computing u ⋅ ∇u and then coarse-graining
        uhr = uhrall[:,:,j]
        vhr = vhrall[:,:,j]

        # uhrh,vhrh,_ = ShallowWaters.add_halo(Float32.(uhr), Float32.(vhr), etahrall[:,:,1], Shr)
        # dudxhr = ShallowWaters.∂x(uhrh, Shr.grid.Δ)
        # dvdyhr = ShallowWaters.∂y(vhrh, Shr.grid.Δ)
        # dot(uhrh, cat(zeros(Float32, 1, nuyhr), dudxhr, dims=1))
        # dot(cat(zeros(Float32, nvxhr, 1), dvdyhr, dims=2), vhrh)

        # coarse-graining and then computing u ⋅ ∇u
        # ucg= ShallowWaters.coarse_grain_u(uhr[:,:,j], Shr.grid.nx, Slr)
        # vcg= ShallowWaters.coarse_grain_v(vhr[:,:,j], Shr.grid.nx, Slr)
        # ucgh,vcgh,_ = ShallowWaters.add_halo(Float32.(ucg), Float32.(vcg), Shr.Prog.η, Slr)
        # dudxlr = ShallowWaters.∂x(ucgh, Slr.grid.Δ)
        # dvdylr = ShallowWaters.∂y(vcgh, Slr.grid.Δ)
        # dot(ulrh, cat(zeros(Float32, 1, nuylr), dudxlr, dims=1))
        # dot(cat(zeros(Float32, nvxlr, 1), dvdylr, dims=2), vlrh)

    end

end