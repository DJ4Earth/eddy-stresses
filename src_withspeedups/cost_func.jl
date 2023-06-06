# this will have just one function, which will compute the cost function given
# a data value and the computed state at that given timestep

function data_misfit(d, u, v, eta)

    @show size(u)
    @show size(v)
    @show size(eta)
    return dot((d - [u; v; eta]), (d - [u; v; eta]))

end

function energy(grid, u, v) 

    return sum(u.^2 .+ v.^2) / (grid.nx * grid.ny)
    
end