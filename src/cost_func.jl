# this will have just one function, which will compute the cost function given
# a data value and the computed state at that given timestep

function data_misfit(t, d, u, v, eta)

    return dot((d[:, t] - [u; v; eta]), (d[:, t] - [u; v; eta]))

end