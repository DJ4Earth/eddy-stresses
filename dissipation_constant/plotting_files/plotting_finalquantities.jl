using Statistics

function compute_energy(energy_true, u, v)

    energy = zeros(522-24)
    for j = 1:(522-24)
        energy[j] = (sum(abs2, u[:,:,j+24]) + sum(abs2, v[:,:,j+24])) / (128*127)
    end

    return mean(energy_true .- energy) / (522 - 24)

end

"""
Computing the difference in the spatially averaged KE. I'm going to compute this is a difference of the spatially
averaged KE fields at each timestep, and then time-average those values
"""
function KE_difference()

    uhrcg10 = cat(uhrcgall[:,:,1:7:1096], uhrcgall[:,:,1097:end];dims=3)
    vhrcg10 = cat(vhrcgall[:,:,1:7:1096], vhrcgall[:,:,1097:end];dims=3)
    energy_hrcg = zeros(522-24)
    for j = 1:(522-24) energy_hrcg[j] = (sum(abs2, uhrcg10[:,:,j+24]) + sum(abs2, vhrcg10[:,:,j+24])) / (127 * 128) end

    energy_zb10 = compute_energy(energy_hrcg, uzb10, vzb10)
    energy_noparam10 = compute_energy(energy_hrcg, unoparam10, vnoparam10)
    energy_tenday10 = compute_energy(energy_hrcg, u10s10, v10s10)
    energy_twentyday10 = compute_energy(energy_hrcg, u20s10, v20s10)
    energy_thirtyday10 = compute_energy(energy_hrcg, u30s10, v30s10)
    energy_multi210 = compute_energy(energy_hrcg, umulti210, vmulti210)
    energy_multi3more10 = compute_energy(energy_hrcg, umulti3more10, vmulti3more10)
    energy_multi1010 = compute_energy(energy_hrcg, umulti1010, vmulti1010)
    energy_multi2010 = compute_energy(energy_hrcg, umulti2010, vmulti2010)

end

"""
RMSE of the sea surface height in the models. This is just a spatial sum of the difference in the time-averaged ssh fields
"""
function compute_ssh_rmse(etatrue, etamodel)

    etamodel_mean = mean(etamodel; dims=3)
    error = sqrt( sum( (etatrue .- etamodel_mean).^2 ./ 128^2 ) )

    return error
end

"""
computing the error in the time-averaged ssh between coarse-resolution models and the 
coarse-grained data
"""
function timeavg_ssh_rmse()

    etahrcg10 = cat(etahrcgall[:,:,1:7:1096], etahrcgall[:,:,1097:end];dims=3)
    etahrcgmean = mean(etahrcg10; dims=3)[:,:,1]

    denom = sqrt(sum( (etahrcgmean).^2 ./ (128^2) ))
    etadiff_zb20 = compute_ssh_rmse(etahrcgmean, etazb10) / denom
    etadiff_noparam = compute_ssh_rmse(etahrcgmean, etanoparam10) / denom
    etadiff_tenday = compute_ssh_rmse(etahrcgmean, eta10s10) / denom
    etadiff_twentyday = compute_ssh_rmse(etahrcgmean, eta20s10) / denom
    etadiff_thirtyday = compute_ssh_rmse(etahrcgmean, eta30s10) / denom
    etadiff_multi2 = compute_ssh_rmse(etahrcgmean, etamulti210) / denom
    etadiff_multi3 = compute_ssh_rmse(etahrcgmean, etamulti3more10) / denom
    etadiff_multi10 = compute_ssh_rmse(etahrcgmean, etamulti1010) / denom
    etadiff_multi20 = compute_ssh_rmse(etahrcgmean, etamulti2010) / denom

end