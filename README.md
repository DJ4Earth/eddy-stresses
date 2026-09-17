This is the supporting repo for "An online-learned neural network-based parameterization in a differentiable shallow water model". It contains scripts for the different experiments that were run, functions used for creating figures, and the tuned weights. In addition to this repo, one can find the synthetic data files at the Zenodo source **need to add link**. Some of the larger data files used, including the fine-resolution runs, are only found in the Zenodo folder. 

These online experiments are designed to work alongside [ShallowWaters.jl](https://github.com/milankl/ShallowWaters.jl). All experiments were run with the following ShallowWaters model configuration:

```julia
ShallowWaters.Parameter(T=T,
        output=false,
        L_ratio=1,
        g=9.81,
        H=500,
        cfl=.898,
        wind_forcing_x="double_gyre",
        Lx=3840e3,
        seasonal_wind_x=false,
        topography="flat",
        adv_scheme="Sadourny",
        bc="nonperiodic",
        bottom_drag="quadratic",
        tracer_advection=false,
        tracer_relaxation=false,
        N=1,
        α=2,
        nx=128,
        Ndays=Ndays
    );
```
This sets the model up to be a flat bottomed, non-periodic domain with wind-driven double barotropic gyre circulation. The forcing functions to run with either the baseline ZB20 (Eq. 5 [here](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2020GL088376)) or the convolutional neural network source terms are within ShallowWaters.jl. 

## Model setup for experiments

The model setup that was (generally) used for all coarse-resolution, neural network experiments is 

```julia

S = ShallowWaters.model_setup(T=Float64,
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
zb_forcing_dissipation=false,
zb_filtered=true,
nn_forcing_momentum=false,
nn_forcing_dissipation=true,
N=1,
α=2,
nx=128,
Ndays=Ndays
);

```

where `Ndays` varies across experiments (in all models this is the parameter determining how many days to integrate). For the high-resolution data, the only change is to set `nx = 1024`, a factor of 8 increase to the number of spatial cells. 

The neural network forcing term used for the online experiments has been added to ShallowWaters as an extension. To activate the extension, one needs to use Lux.jl after using ShallowWaters.jl

```julia

using ShallowWaters

using Lux

```

Before `using Lux`, ShallowWaters will effectively ignore the neural network forcing term. In early experiments, such as the initial offline problem, the CFL condition was left with `cfl = 0.9`, the default value. In later online learning experiments, the CFL condition was modified to be `cfl = .898`, as `cfl = .9` led to a *slight* mismatch with the high-resolution states used for data. 

### Neural network setup

The neural network in all experiments is initialized as 

```julia

Su_dims = [3,25,25,1]
Sv_dims = [3,25,25,2]

Su_layers = Lux.Chain(
(
Lux.Conv((5,5), Su_dims[i] => Su_dims[i+1], (i == (length(Su_dims)-1) ? identity : gelu); pad=SamePad(),use_bias=Lux.False())
for i in 1:(length(Su_dims)-1)
)...
)


Sv_layers = Lux.Chain(
(
Lux.Conv((5,5), Sv_dims[i] => Sv_dims[i+1], (i == (length(Sv_dims)-1) ? identity : gelu); pad=SamePad(),use_bias=Lux.False())
for i in 1:(length(Sv_dims)-1)
)...
)


model_Su = Lux.setup(Random.default_rng(), Su_layers)
model_Sv = Lux.setup(Random.default_rng(), Sv_layers)

```

within [`ShallowWaters.jl/ext/`](https://github.com/milankl/ShallowWaters.jl/blob/282fa603701d600949c15f5f1333d2543b5cb1e4/ext/ShallowWatersNNForcing.jl#L8). This setup which gives a total of 36875 tunable weights. When inserting tuned weights into `model_Su` and `model_Sv` we use the following for loop, 

```julia

current = 1
for m in (S.Diag.CNNVars.model_Su, S.Diag.CNNVars.model_Sv)
	for layers in m[1]
		for array in layers
			sz = prod(size(array))
			array .= reshape(param_guess[current:(current + sz - 1)], size(array)...)
			current += sz
		end
	end
end

```

which first inserts weights into `model_Su` and then inserts weights into `model_Sv`. If the size of `model_Su` and `model_Sv` are changed by changing `Su_dims` and/or `Sv_dims`, then the tuned weights provided here will no longer "fit" the neural networks and new weights are needed. 

**Caution regarding the use of the neural network forcing term**: the default weights chosen by Lux.jl will lead to model divergence (in my experience this happened every time). How long it takes for ShallowWaters to diverge varies, but generally a few days is enough. In order to recreate stable weights one needs to run a new optimization problem.

## Quick start

To run ShallowWaters with some of the tuned weights from the online learning experiments, one can run 

```julia
	using ShallowWaters#sw/nn_forcing
	using Lux

	Ndays = 5
    S = ShallowWaters.model_setup(T=Float64,
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
        nx=128,
        Ndays=Ndays,
        initial_cond="rest"
    );

    onlineweights = load_object("./dissipation_constant/tuned_weights/ensemble_initial_condition/result_multistate_1-4-8-13-18-23-28-33-38-41-44-48-53-58-63-68-73-78-83-86daystart_3dayoptimization_initialweightsmulti3daystate_20iterations.jld2").solution

    current = 1
    for m in (S.Diag.CNNVars.model_Su, S.Diag.CNNVars.model_Sv)
        for layers in m[1]
            for array in layers
                    sz = prod(size(array))
                    array .= reshape(onlineweights[current:(current + sz - 1)], size(array)...)
                    current += sz
            end
        end
    end

    P = ShallowWaters.time_integration(Sonline);

```
