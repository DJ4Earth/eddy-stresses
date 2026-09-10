This is the supporting repo for "An online-learned neural network-based parameterization in a differentiable shallow water model". It contains scripts for the different experiments that were run, functions used for creating figures, and the tuned weights. In addition to this repo, one can find the synthetic data files at the Zenodo source **need to add link**. The relevant files included here are:

1. eddy_paper.jl
2. ./dissipation_constant/states_madnlp.jl
3. ./dissipation_constant/ensemble_states_madnlp.jl

These online experiments are designed to work alongside [ShallowWaters.jl](https://github.com/milankl/ShallowWaters.jl), a single layer ocean model. All experiments were run with the following ShallowWaters model configuration:

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
        zb_forcing_momentum=false,
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
This sets the model up to be a flat bottomed, non-periodic domain with wind-driven double barotropic gyre circulation. The forcing functions to run with either the baseline ZB20 (Eq. 5 [here](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2020GL088376)) or the convolutional neural network source terms are within ShallowWaters.jl. Important to note that the NN forcing term is setup as an extension: to use it one needs to run 
```julia
using ShallowWaters#main
using Lux
```
All tuned weights from experiments are provided here in 

`./dissipation_constant/tuned_weights`

and are split into folders depending on if they were the result of a single initial condition experiment or an ensemble initial condition experiment. 

The single initial condition experiments were all run in the script 

`./dissipation_constant/states_madnlp.jl`

and the ensemble initial condition experiments were run in

`./dissipation_constant/ensemble_states_madnlp.jl`.

The initial offline weights were learned in

`./dissipation_constant/initial_weights_offlineproblem.jl`

To run ShallowWaters with a NN forcing term, one can do

```julia
	using ShallowWaters#main
	using Lux
 	Ndays = 10
    Ponline = ShallowWaters.Parameter(T=Float64,
        output=false,
        output_vars=["u", "v", "η", "ζ"],
        L_ratio=1,
        g=9.81,
        H=500,
        cfl=.898,
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

    Sonline = ShallowWaters.model_setup(Ponline);

    onlineweights = load_object("./dissipation_constant/tuned_weights/ensemble_initial_condition/result_multistate_1-4-6-8-10-13-15-18-23-28-33-38-41-44-48-51-53-58-63-65-68-73-78-83-86-88daystart_2dayoptimization_initialweightsmulti3daystate_20iterations.jld2").solution;


	# after loading the weights into variable onlineweights, this for loop places them into the NN structure created when initializing the model. Important to note
	# that all of the tuned weights are specific to a certain NN size -- changing the size will cause the weights to be incompatible. 
  	current = 1
	for m in (Sonline.Diag.CNNVars.model_Su, Sonline.Diag.CNNVars.model_Sv)
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
