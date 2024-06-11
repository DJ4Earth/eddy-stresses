using Oceananigans
using Oceananigans: time_step!
using Printf

Nx = 128
Ny = 128

# A spherical domain
# grid = LatitudeLongitudeGrid(size = (Nx, Ny, 1),
#                              longitude = (-30, 30),
#                              latitude = (30, 35),
#                              z = (-500, 0))
grid = RectilinearGrid(size = (Nx, Ny, 1),
                    topology = (Periodic, Periodic, Bounded),
                    extent = (3840e3, 3840e3, 500)              
)

@show surface_wind_stress_parameters = (τ₀ = 0.12 / 1000,
                                        L = grid.Ly
)

@inline surface_wind_stress(y, p) = p.τ₀ * (cos(2 * pi * (y/p.L - 0.5)) + 2 * sin(2 * pi *(y/p.L - 0.5)))

# @inline surface_wind_stress(λ, φ, t, p) = p.τ₀ * cos(2π * (φ - p.φ₀) / p.L)

surface_wind_stress_bc = FluxBoundaryCondition(surface_wind_stress,
    parameters = surface_wind_stress_parameters
)

μ = 1/(60 * 24 * 60 * 60)#60days
# @inline u_bottom_drag(i, j, grid, clock, fields, μ) = @inbounds - μ * fields.u[i, j, 1]
# @inline v_bottom_drag(i, j, grid, clock, fields, μ) = @inbounds - μ * fields.v[i, j, 1]
@inline u_quadratic_drag(x, y, t, u, v) = - 1e-5 * u * sqrt(u^2 + v^2)
@inline v_quadratic_drag(x, y, t, u, v) = - 1e-5 * v * sqrt(u^2 + v^2)

u_bottom_drag_bc = FluxBoundaryCondition(u_quadratic_drag, 
    field_dependencies = (:u, :v)
)
v_bottom_drag_bc = FluxBoundaryCondition(v_quadratic_drag,
    field_dependencies = (:u, :v)
)

u_bcs = FieldBoundaryConditions(top = surface_wind_stress_bc, bottom = u_bottom_drag_bc)
v_bcs = FieldBoundaryConditions(bottom = v_bottom_drag_bc)

νh₀ = 5e3 * (60 / grid.Nx)^2
constant_horizontal_diffusivity = HorizontalScalarDiffusivity(ν = νh₀)

model = HydrostaticFreeSurfaceModel(; grid,
                                    momentum_advection = VectorInvariant(),
                                    free_surface = ImplicitFreeSurface(gravitational_acceleration=0.1),
                                    coriolis = BetaPlane(rotation_rate=7.27e-5, latitude=35),
                                    boundary_conditions = (u=u_bcs, v=v_bcs),
                                    closure = constant_horizontal_diffusivity,
                                    tracers = nothing,
                                    buoyancy = nothing
)

time_step!(model, 6*60*60)#6hours)

@show maximum(model.velocities.u)
@show maximum(model.velocities.v)

# u, v = model.velocities

# set!(model, u=u.data, v=v.data)

time_step!(model, 6*60*60)

@show maximum(model.velocities.u)
@show maximum(model.velocities.v)

# u, v = model.velocities

# set!(model, u=u.data, v=v.data)

# time_step!(model, 6*60*60)

# @show maximum(model.velocities.u)
# @show maximum(model.velocities.v)
