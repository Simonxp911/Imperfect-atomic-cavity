

# ================================================
#   Calculate steady state atomic expectation values
# ================================================
function calc_σ_ss(Δ, Gnm, drivemode)
    # The transmission amplitude for some drive and identical detection 
    # is found by a matrix inversion to solve the EoMs in the steady state
    return -((Δ + 1im)*I + Gnm)\drivemode
end


# ================================================
#   Calculate E-field for finite array
# ================================================
function calc_atomic_Efield_fin(r, array, σ_ss, dipoleMoment, Gmat_rn=nothing)
    # Get the GF matrix (evaluated at r - r_n for each atom n)
    if isnothing(Gmat_rn) Gmat_rn = get_Gmat_rn(r, array) end
    
    # # Calculate and return the E-field (in units of d, i.e. return Ed)
    # return sum(Gmat_rn.*σ_ss) * dipoleMoment
    
    # Calculate and return the E-field (in units of d, i.e. return Ed) [optimized calculation]
    atomic_Efield_fin = zeros(ComplexF64, 3)
    for i in 1:3, j in eachindex(σ_ss)
        atomic_Efield_fin[i] += (Gmat_rn[j][i, 1]*dipoleMoment[1] + Gmat_rn[j][i, 2]*dipoleMoment[2] + Gmat_rn[j][i, 3]*dipoleMoment[3])*σ_ss[j]
    end
    return atomic_Efield_fin
end


function calc_total_Efield_fin(r, array, σ_ss, drive_type, w0, dipoleMoment, Gmat_rn=nothing)
    # Incoming E-field (assuming the drive to have dipoleMoment as its polarization)
    Ed_in = get_drivemode(drive_type, r, w0)*dipoleMoment
    
    # Get the atomic contribution to the E-field
    Ed_at = calc_atomic_Efield_fin(r, array, σ_ss, dipoleMoment, Gmat_rn)
    
    # Sum the contributions and return (in units of 1/d, i.e. return Ed)
    return Ed_in + Ed_at
end


# ================================================
#   Calculate transmission amplitudes
# ================================================
function calc_transAmpl_fin(Δ, array, Gnm, drivemode, SP)
    # Get the steady state
    σ_ss = calc_σ_ss(Δ, Gnm, drivemode)
    
    if SP.DeP.detec_type == "drive_mode"
        # The transmission amplitude for a paraxial drive and detection in that mode
        return 1 + 3π*1im/ωa^2*drivemode'*σ_ss
    
    # The transmission amplitude for some drive and different choices of detection 
    # as calculated by a direct integration of the E-field (the usual analytic expression is only for paraxial detection modes)
    elseif SP.DeP.detec_type == "integrated_drive_mode"
        # Calculate E-field and drive on integration plane, as well as the detection mode
        Ed = calc_total_Efield_fin.(SP.DeP.integration_plane, Ref(array), Ref(σ_ss), SP.DrP.drive_type, SP.DrP.w0, Ref(SP.EP.dipoleMoment))
        drive = get_drivemode.(SP.DrP.drive_type, SP.DeP.integration_plane, SP.DrP.w0) .* Ref(SP.EP.dipoleMoment)
        detec_type = drive
        
    elseif SP.DeP.detec_type == "flat_mode_on_detection_plane"
        # Calculate E-field and drive on detection plane, as well as the detection mode
        Ed = calc_total_Efield_fin.(SP.DeP.detection_plane, Ref(array), Ref(σ_ss), SP.DrP.drive_type, SP.DrP.w0, Ref(SP.EP.dipoleMoment))
        drive = get_drivemode.(SP.DrP.drive_type, SP.DeP.detection_plane, SP.DrP.w0) .* Ref(SP.EP.dipoleMoment)
        detec_type = ones(length(SP.DeP.detection_plane)) .* Ref(SP.EP.dipoleMoment) / sqrt(length(SP.DeP.detection_plane)*SP.DeP.dx*SP.DeP.dy)
        
    elseif SP.DeP.detec_type == "incoming_mode_on_detection_plane"
        # Calculate E-field and drive on detection plane, as well as the detection mode
        Ed = calc_total_Efield_fin.(SP.DeP.detection_plane, Ref(array), Ref(σ_ss), SP.DrP.drive_type, SP.DrP.w0, Ref(SP.EP.dipoleMoment))
        drive = get_drivemode.(SP.DrP.drive_type, SP.DeP.detection_plane, SP.DrP.w0) .* Ref(SP.EP.dipoleMoment)
        # detec_type = Ed./norm.(Ed)
        detec_type = Ed./norm.(Ed) / sqrt(length(SP.DeP.detection_plane)*SP.DeP.dx*SP.DeP.dy)
        
    else
        throw(DomainError(SP.DeP.detec_type, "This detec_type has not been implemented in calc_transAmpl_fin"))
    end
    
    # The transmission is then calculated as the integral of the product of detec_type^\dagger and Ed 
    # (assuming detec_type to be normalized such that integral of detec_type^\dagger*detec_type is unity)
    normalization = sum(adjoint.(detec_type).*drive)*SP.DeP.dx*SP.DeP.dy
    return sum(adjoint.(detec_type).*Ed)*SP.DeP.dx*SP.DeP.dy/normalization
end


function calc_transCoef_fin(Δ, array, Gnm, drivemode, SP, Gmat_rn_plane=nothing)
    if SP.DeP.detec_type ∈ ("drive_mode", "integrated_drive_mode", "flat_mode_on_detection_plane", "incoming_mode_on_detection_plane")
        return abs2(calc_transAmpl_fin(Δ, array, Gnm, drivemode, SP))
        
    elseif SP.DeP.detec_type == "intensity_on_detection_plane"
        # Get the steady state
        σ_ss = calc_σ_ss(Δ, Gnm, drivemode)
        
        # Calculate E-field and drive on detection plane
        Ed = calc_total_Efield_fin.(SP.DeP.detection_plane, Ref(array), Ref(σ_ss), SP.DrP.drive_type, SP.DrP.w0, Ref(SP.EP.dipoleMoment), Gmat_rn_plane)
        drive = get_drivemode.(SP.DrP.drive_type, SP.DeP.detection_plane, SP.DrP.w0) .* Ref(SP.EP.dipoleMoment)
        
        # The transmission is then calculated as the integral of the product of Ed^\dagger and Ed 
        normalization = sum(adjoint.(drive).*drive)*SP.DeP.dx*SP.DeP.dy
        return real(sum(adjoint.(Ed).*Ed)*SP.DeP.dx*SP.DeP.dy/normalization)
    
    else
        throw(DomainError(SP.DeP.detec_type, "This detec_type has not been implemented in calc_transCoef_fin"))
    end
    
end


function calc_transAmpl_inf(lattice_type, N_sheets, a, L, Δ, dipoleMoment)
    # The transmission amplitude for a normal-incidence, plane-wave 
    # i.e. the system is driven by a k=0 plane-wave and only the 
    # k=0 plane-wave component of the emitted light is detected
    
    # Calculate tildeΓ0 directly and get the self-energy
    tildeΓ0 = 3π/(ωa*a)^2
    Σ0 = -ana_FT_GF(lattice_type, a, dipoleMoment, dipoleMoment, 0.0, 0.0, 0.0, 1e-5)
        
    if N_sheets == 1
        return 1 - 1im*tildeΓ0/(Δ - Σ0)
        
    elseif N_sheets == 2
        # Get the self-energy for L
        ΣL = -ana_FT_GF(lattice_type, a, dipoleMoment, dipoleMoment, 0.0, 0.0, L, 1e-5)
        
        # Get the trigonometric factors
        CL_k = cos(ωa*L)
        
        # Calculate and return the transmission amplitude
        return 1 - 1im*tildeΓ0*((1 + CL_k)/(Δ - (Σ0 + ΣL)) + (1 - CL_k)/(Δ - (Σ0 - ΣL)))
        
    else
        throw(ArgumentError("N_sheets = $N_sheets has not been implemented in calc_transAmpl_inf"))
    end
end


function calc_transAmpl_inf(lattice_type, N_sheets, a, L, Δ, dipoleMoment, dipoleMoment_label, drive_type, w0, k_n)
    # The transmission amplitude for the infinite system with some drive and identical detection
    # It is assumed that the drive has the same symmetries as the FT GF, 
    # such that the integration can go over the first octant only
    
    if N_sheets ∉ [1, 2] throw(ArgumentError("N_sheets = $N_sheets has not been implemented in calc_transAmpl_inf")) end
    
    # Calculate tildeΓ0 directly and prepare Σ0 and ΣL
    tildeΓ0 = 3π/(ωa*a)^2
    Σ0 = prepare_Σ(lattice_type, a, 0.0, k_n, dipoleMoment, dipoleMoment_label)
    if N_sheets == 2 ΣL = prepare_Σ(lattice_type, a, L, k_n, dipoleMoment, dipoleMoment_label) end
    
    # Set the k-space/BZ resolution
    dk = (π/a)/(k_n - 1)
    k_range = (0:k_n-1)*dk
    
    # Integrate over the BZ
    integrals = zeros(ComplexF64, N_sheets)
    for (i, kx) in enumerate(k_range), (j, ky) in enumerate(k_range)
        # Only consider the first octant
        if ky > kx continue end
        
        # Get the squared kz and skip if negative
        kz2 = ωa^2 - kx^2 - ky^2
        if kz2 < 0 continue end
        
        # Get the drive and skip if it's is small
        drive = get_drive_kspace(drive_type, w0, dk, [kx, ky])
        if drive < 1e-3*maximum(drive) continue end
        
        # Define multiplicities of points in the octant
        if     kx == ky == 0       mult = 1
        elseif ky == 0 || kx == ky mult = 4
        else                       mult = 8 end
        
        # Get kz
        kz = sqrt(kz2)
        
        # Get the polarization overlap factor 
        pol_overlap = 1 - (kx^2 + ky^2)/(2*ωa^2)
        
        if N_sheets == 1
            det = Δ - Σ0[i, j]
            integrals[1] += mult * drive^2*ωa/kz*pol_overlap^2/det
            
        elseif N_sheets == 2
            # Get the trigonometric factors
            CL_k = cos(kz*L)
            
            # Complex detuning factors
            det_p = Δ - (Σ0[i, j] + ΣL[i, j])
            det_m = Δ - (Σ0[i, j] - ΣL[i, j])
            
            integrals[1] += mult * drive^2*ωa/kz*pol_overlap^2*(1 + CL_k)/det_p
            integrals[2] += mult * drive^2*ωa/kz*pol_overlap^2*(1 - CL_k)/det_m
        end
    end
    
    # Normalize integrals
    integrals .*= dk^2/(2π)^2
    
    # Calculate and return the transmission amplitude
    return 1.0 - 1im*tildeΓ0*sum(integrals)
end


# ================================================
#   Make scans of the transmission amplitudes and do statistics for them
# ================================================
function scan_transCoef_fin(SP)
    # Check if the scan has already been performed
    postfix = get_postfix_Tscan(SP.AP.lattice_type, SP.AP.N_sheets, SP.AP.radius, SP.AP.cut_corners, SP.AP.a, SP.AP.L, SP.AP.ff, SP.pos_unc_ratio, SP.AP.N_inst, SP.DrP.drive_type, SP.w0_ratio, SP.EP.dipoleMoment_label, SP.DeP.detec_type, SP.DeP.detec_radius, SP.DeP.detec_z, SP.Delta_specs)
    filename_ts = "Tscan_" * postfix
    data = check_if_already_calculated(save_dir, [filename_ts])
    if length(data) == 1 return data[1] end
    
    # Prepare coupling matrices (to avoid calculation these for each value of detuning)
    Gnm = get_Gnm.(SP.AP.array, SP.AP.N, Ref(SP.EP.dipoleMoment))
    if SP.DeP.detec_type ∈ ("integrated_drive_mode",)
        Gmat_rn_plane = [get_Gmat_rn.(SP.DeP.integration_plane, Ref(arr)) for arr in SP.AP.array]
    elseif SP.DeP.detec_type ∈ ("flat_mode_on_detection_plane", "incoming_mode_on_detection_plane", "intensity_on_detection_plane")
        Gmat_rn_plane = [get_Gmat_rn.(SP.DeP.detection_plane, Ref(arr)) for arr in SP.AP.array]
    else 
        Gmat_rn_plane = fill(nothing, SP.AP.N_inst)
    end
    
    # Perform scan
    Tscan = calc_transCoef_fin.(reshape(SP.Delta_range, 1, SP.Delta_specs[3]), 
                                reshape(SP.AP.array, SP.AP.N_inst, 1), 
                                reshape(Gnm, SP.AP.N_inst, 1), 
                                reshape(SP.DrP.drivemode, SP.AP.N_inst, 1),
                                Ref(SP),
                                reshape(Gmat_rn_plane, SP.AP.N_inst, 1))
    
    # Save the scan
    save_as_txt(Tscan, save_dir, filename_ts)
    
    return Tscan
end


function scan_transAmpl_inf(SP)
    # Check if the scan has already been performed
    postfix = get_postfix_tscan_inf(SP.AP.lattice_type, SP.AP.N_sheets, SP.DrP.drive_type, SP.w0_ratio, SP.k_n ,SP.EP.dipoleMoment_label, SP.Delta_specs)
    filename_ts_k0 = "tscan_inf_k0_" * postfix
    filename_ts_k  = "tscan_inf_k_"  * postfix
    data = check_if_already_calculated(save_dir, [filename_ts_k0, filename_ts_k], ComplexF64)
    if length(data) == 2 return unpack_tscan_inf.(data) end
    
    # Perform scan
    t_inf_k0 = calc_transAmpl_inf.(SP.AP.lattice_type, SP.AP.N_sheets, SP.AP.a, SP.AP.L, SP.Delta_range, Ref(SP.EP.dipoleMoment))
    t_inf_k  = calc_transAmpl_inf.(SP.AP.lattice_type, SP.AP.N_sheets, SP.AP.a, SP.AP.L, SP.Delta_range, Ref(SP.EP.dipoleMoment), SP.EP.dipoleMoment_label, SP.DrP.drive_type, SP.DrP.w0, SP.k_n)
    
    
    # # Save the scan
    # data_k0 = pack_tscan_inf(t_inf_k0)
    # data_k  = pack_tscan_inf(t_inf_k)
    # save_as_txt(data_k0, save_dir, filename_ts_k0)
    # save_as_txt(data_k , save_dir, filename_ts_k)
    
    return t_inf_k0, t_inf_k
end


function scan_statistics(scan)
    # Get a slice for each value of Delta on which to do statistics
    Delta_slices = eachslice(scan, dims=2)
    
    # Find mean and standard deviation for each slice
    means = mean.(Delta_slices)
    stds  = std.(Delta_slices)
    
    return means, stds
end