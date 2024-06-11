using Serialization, Plots
using HDF5, LaTeXStrings
using Parameters, UnPack, SparseArrays
using ColorSchemes

include("/Users/swilliamson/Documents/GitHub/eddy-stresses/src_withspeedups/init_structs.jl")

function deserialize(x)
    s = IOBuffer(x)
    Serialization.deserialize(s)
end

# if !HDF5.ishdf5("primal_chkp.h5")
#     error("File not found primal_chkp.h5")
# end
# if !HDF5.ishdf5("adjoint_chkp.h5")
#     error("File not found adjoint_chkp.h5")
# end

# prim_fid = h5open("energyex_primal_chkp_30days.h5", "r")
# adj_fid = h5open("energyex_adjoint_chkp_30days.h5", "r")

# # Read first adjoint field
# blob = read(adj_fid["1"])
# adj_chkp = deserialize(blob)

function create_adjoint_gif()
    
    # prim_fid = h5open("nxny128_365_primal_chkp.h5", "r")
    adj_fid = h5open("nxny128_365_adjoint_chkp.h5", "r")

    # eta_anim = Animation()
    # u_anim = Animation()
    # v_anim = Animation()

    # for j = 15000:-100:1

        blob = read(adj_fid[string(1)])
        adj_chkp = deserialize(blob)
        
        # frame(eta_anim, heatmap(reshape(adj_chkp.eta, 128, 128)', title=L"\partial \mathcal{E}(t_f)/\partial \eta(%$j)" , 
        #     clim=(-150, 150), xlabel=L"x", ylabel=L"y", c=:balance, dpi=300))
        # frame(u_anim, heatmap(reshape(adj_chkp.u, 127, 128)', title=L"\partial \mathcal{E}(t_f)/\partial u(%$j)", 
        #     clim=(-10000, 10000), xlabel=L"x", ylabel=L"y", c=:balance, dpi=300))
        # frame(v_anim, heatmap(reshape(adj_chkp.v, 128, 127)', title="Energy sensitivity w.r.t v"))

        for_billy = heatmap(reshape(adj_chkp.eta, 128, 128)', title=L"\partial J/\partial \eta(%$j)" , 
             clim=(-150, 150), xlabel=L"x", ylabel=L"y", c=:balance, dpi=300)
    


    # end

    savefig("for_poster.png", for_billy)
    # gif(eta_anim, "deta_integration_071323_newcolor_neg150pos150_fps7.gif", fps = 7)
    # gif(u_anim, "du_integration_071323_newcolor_fps7.gif", fps = 7)
    # gif(v_anim, "dv_integration_days=365-061523.gif", fps = 30)

end


create_adjoint_gif()