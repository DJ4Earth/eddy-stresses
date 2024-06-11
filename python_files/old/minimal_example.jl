# now going to try a sample adjoint calculation. just going to look at the sensitivity of the final 
# displacement in center of grid w.r.t. initial displacement. 

using Enzyme

mutable struct gyre
    lastη::Matrix{Float64} 
    η::Matrix{Float64}
end

mutable struct gyre_Parameters 
    dt::Float64 
end

function advance(u_v_eta, parameters) 

    Nx = 62
    Ny = 62

    dt = parameters.dt

    η = copy(u_v_eta.η)

    nextη = zeros(Nx, Ny)

    # update fields 

    for j = 2:Nx - 1
        for k = 2:Ny - 1 

            nextη[j,k] = η[j,k] #+ dt * 0.0 

        end
    end

    copyto!(u_v_eta.η, nextη)
    copyto!(u_v_eta.lastη, η)

    return nothing 

end 

dt = (0.5 * 1e3 * 20.0) / (2 * sqrt(9.81 * 5000.0))

gyre_parameters = gyre_Parameters(
    dt
)

u_v_eta = gyre(
    zeros(62, 62),
    zeros(62, 62),
)

# structure for adjoint variables 
initial_adjoint = deepcopy(u_v_eta)

init_adj_eta = 0.0 .* u_v_eta.η
init_adj_eta[31, 31] = 1.0

initial_adjoint.lastη = zeros(62, 62)
initial_adjoint.η = copy(init_adj_eta)

function ad_calc(ad_struct, state, params)

        @show ad_struct.η[31, 31]
        @show sum(abs.(ad_struct.η)), state.η[31, 31]

        autodiff(advance, Const, Duplicated(state, ad_struct), Const(params))

        @show ad_struct.η[31,31]
        @show sum(abs.(ad_struct.η)), state.η[31, 31]

        return ad_struct
    
end 
    
@show sum(abs.(initial_adjoint.η))
ad_u_v_eta = ad_calc(initial_adjoint, u_v_eta, gyre_parameters);