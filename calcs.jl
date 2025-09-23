

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
function calc_atomic_Efield_fin(r, array, σ_ss, e1, Gmat_rn=nothing)
    # Get the GF matrix (evaluated at r - r_n for each atom n)
    if isnothing(Gmat_rn) Gmat_rn = get_Gmat_rn(r, array) end
    
    # Calculate and return the E-field (in units of d, i.e. return Ed)
    return sum(Gmat_rn.*σ_ss) * e1
end


function calc_total_Efield_fin(r, array, σ_ss, drive_type, w0, e1, Gmat_rn=nothing)
    # Incoming E-field (assuming the drive to have e1 as its polarization)
    Ed_in = get_drivemode(drive_type, r, w0)*e1
    
    # Get the atomic contribution to the E-field
    Ed_at = calc_atomic_Efield_fin(r, array, σ_ss, e1, Gmat_rn)
    
    # Sum the contributions and return (in units of 1/d, i.e. return Ed)
    return Ed_in + Ed_at
end


# ================================================
#   Calculate transmission amplitudes
# ================================================
function calc_transAmpl_fin(Δ, array, Gnm, drivemode, SP)
    # Get the steady state
    σ_ss = calc_σ_ss(Δ, Gnm, drivemode)
    
    if SP.detec_mode == "drive_mode"
        # The transmission amplitude for a paraxial drive and detection in that mode
        return 1 + 3π*1im/wa^2*drivemode'*σ_ss
    
    # The transmission amplitude for some drive and different choices of detection 
    # as calculated by a direct integration of the E-field (the usual analytic expression is only for paraxial detection modes)
    elseif SP.detec_mode == "integrated_drive_mode"
        # Calculate E-field and drive on integration plane, as well as the detection mode
        Ed = calc_total_Efield_fin.(SP.integration_plane, Ref(array), Ref(σ_ss), SP.drive_type, SP.w0, Ref(SP.e1))
        drive = get_drivemode.(SP.drive_type, SP.integration_plane, SP.w0) .* Ref(SP.e1)
        detec_mode = drive
        
    elseif SP.detec_mode == "flat_mode_on_detection_plane"
        # Calculate E-field and drive on detection plane, as well as the detection mode
        Ed = calc_total_Efield_fin.(SP.detection_plane, Ref(array), Ref(σ_ss), SP.drive_type, SP.w0, Ref(SP.e1))
        drive = get_drivemode.(SP.drive_type, SP.detection_plane, SP.w0) .* Ref(SP.e1)
        detec_mode = ones(length(SP.detection_plane)) .* Ref(SP.e1) / sqrt(length(SP.detection_plane)*SP.dx*SP.dy)
        
    elseif SP.detec_mode == "incoming_mode_on_detection_plane"
        # Calculate E-field and drive on detection plane, as well as the detection mode
        Ed = calc_total_Efield_fin.(SP.detection_plane, Ref(array), Ref(σ_ss), SP.drive_type, SP.w0, Ref(SP.e1))
        drive = get_drivemode.(SP.drive_type, SP.detection_plane, SP.w0) .* Ref(SP.e1)
        # detec_mode = Ed./norm.(Ed)
        detec_mode = Ed./norm.(Ed) / sqrt(length(SP.detection_plane)*SP.dx*SP.dy)
        
    else
        throw(DomainError(SP.detec_mode, "This detec_mode has not been implemented in calc_transAmpl_fin"))
    end
    
    # The transmission is then calculated as the integral of the product of detec_mode^\dagger and Ed 
    # (assuming detec_mode to be normalized such that integral of detec_mode^\dagger*detec_mode is unity)
    normalization = sum(adjoint.(detec_mode).*drive)*SP.dx*SP.dy
    return sum(adjoint.(detec_mode).*Ed)*SP.dx*SP.dy/normalization
end


function calc_transCoef_fin(Δ, array, Gnm, drivemode, SP, Gmat_rn_plane=nothing)
    if SP.detec_mode ∈ ("drive_mode", "integrated_drive_mode", "flat_mode_on_detection_plane", "incoming_mode_on_detection_plane")
        return abs2(calc_transAmpl_fin(Δ, array, Gnm, drivemode, SP))
        
    elseif SP.detec_mode == "intensity_on_detection_plane"
        # Get the steady state
        σ_ss = calc_σ_ss(Δ, Gnm, drivemode)
        
        # Calculate E-field and drive on detection plane
        Ed = calc_total_Efield_fin.(SP.detection_plane, Ref(array), Ref(σ_ss), SP.drive_type, SP.w0, Ref(SP.e1), Gmat_rn_plane)
        drive = get_drivemode.(SP.drive_type, SP.detection_plane, SP.w0) .* Ref(SP.e1)
        
        # The transmission is then calculated as the integral of the product of Ed^\dagger and Ed 
        normalization = sum(adjoint.(drive).*drive)*SP.dx*SP.dy
        return real(sum(adjoint.(Ed).*Ed)*SP.dx*SP.dy/normalization)
              
    else
        throw(DomainError(SP.detec_mode, "This detec_mode has not been implemented in calc_transCoef_fin"))
    end
    
end


function calc_transAmpl_inf(lattice_type, N_sheets, a, L, Δ, e1)
    # The transmission amplitude for a normal-incidence, plane-wave 
    # i.e. the system is driven by a k=0 plane-wave and only the 
    # k=0 plane-wave component of the emitted light is detected
    
    # Calculate tildeΓ0 directly and get the self-energy
    tildeΓ0 = 3π/(wa*a)^2
    Σ0 = -ana_FT_GF(lattice_type, a, e1, e1, 0.0, 0.0, 0.0, 1e-5)
        
    if N_sheets == 1
        return 1 - 1im*tildeΓ0/(Δ - Σ0)
        
    elseif N_sheets == 2
        # Get the self-energy for L
        ΣL = -ana_FT_GF(lattice_type, a, e1, e1, 0.0, 0.0, L, 1e-5)
        
        # Get the trigonometric factors
        CL_k = cos(wa*L)
        
        # Calculate and return the transmission amplitude
        return 1 - 1im*tildeΓ0*((1 + CL_k)/(Δ - (Σ0 + ΣL)) + (1 - CL_k)/(Δ - (Σ0 - ΣL)))
        
    else
        throw(ArgumentError("N_sheets = $N_sheets has not been implemented in calc_transAmpl_inf"))
    end
end


function calc_transAmpl_inf(lattice_type, N_sheets, a, L, Δ, e1, e1_label, drive_type, w0, k_n)
    # The transmission amplitude for the infinite system with some drive and identical detection
    # It is assumed that the drive has the same symmetries as the FT GF, 
    # such that the integration can go over the first octant only
    
    if N_sheets ∉ [1, 2] throw(ArgumentError("N_sheets = $N_sheets has not been implemented in calc_transAmpl_inf")) end
    
    # Calculate tildeΓ0 directly and prepare Σ0 and ΣL
    tildeΓ0 = 3π/(wa*a)^2
    Σ0 = prepare_Σ(lattice_type, a, 0.0, k_n, e1, e1_label)
    if N_sheets == 2 ΣL = prepare_Σ(lattice_type, a, L, k_n, e1, e1_label) end
    
    # Set the k-space/BZ resolution
    dk = (π/a)/(k_n - 1)
    k_range = (0:k_n-1)*dk
    
    # Integrate over the BZ
    integrals = zeros(ComplexF64, N_sheets)
    for (i, kx) in enumerate(k_range), (j, ky) in enumerate(k_range)
        # Only consider the first octant
        if ky > kx continue end
        
        # Get the squared kz and skip if negative
        kz2 = wa^2 - kx^2 - ky^2
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
        pol_overlap = 1 - (kx^2 + ky^2)/(2*wa^2)
        
        if N_sheets == 1
            det = Δ - Σ0[i, j]
            integrals[1] += mult * drive^2*wa/kz*pol_overlap^2/det
            
        elseif N_sheets == 2
            # Get the trigonometric factors
            CL_k = cos(kz*L)
            
            # Complex detuning factors
            det_p = Δ - (Σ0[i, j] + ΣL[i, j])
            det_m = Δ - (Σ0[i, j] - ΣL[i, j])
            
            integrals[1] += mult * drive^2*wa/kz*pol_overlap^2*(1 + CL_k)/det_p
            integrals[2] += mult * drive^2*wa/kz*pol_overlap^2*(1 - CL_k)/det_m
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
    println("Runnning scan_transCoef_fin")
    
    # Check if the scan has already been performed
    postfix = get_postfix(SP.lattice_type, SP.N_sheets, SP.radius, SP.cut_corners, SP.a, SP.L, SP.ff, SP.pos_unc_ratio, SP.N_inst, SP.drive_type, SP.w0_ratio, SP.e1_label, SP.detec_mode, SP.detec_radius, SP.detec_z, SP.Delta_specs)
    filename_ts = "Tscan" * postfix
    data = check_if_already_calculated(save_dir, [filename_ts], ComplexF64)
    if length(data) == 1 return data[1] end
    
    # Perform scan
    println("Performing scan")
    Tscan = calc_transCoef_fin.(reshape(SP.Delta_range, 1, SP.Delta_specs[3]), 
                                reshape(SP.array, SP.N_inst, 1), 
                                reshape(SP.Gnm, SP.N_inst, 1), 
                                reshape(SP.drivemode, SP.N_inst, 1),
                                Ref(SP),
                                reshape(SP.Gmat_rn_plane, SP.N_inst, 1))
    
    # Save the scan
    save_as_txt(Tscan, save_dir, filename_ts)
    
    return Tscan
end


function scan_transAmpl_inf(SP)
    println("Runnning scan_transAmpl_inf")
    
    # Check if the scan has already been performed
    postfix = get_postfix(SP.lattice_type, SP.N_sheets, SP.drive_type, SP.w0_ratio, SP.k_n ,SP.e1_label, SP.Delta_specs)
    filename_ts_k0 = "tscan_inf_k0" * postfix
    filename_ts_k  = "tscan_inf_k"  * postfix
    data = check_if_already_calculated(save_dir, [filename_ts_k0, filename_ts_k], ComplexF64)
    if length(data) == 2 return unpack_tscan_inf.(data) end
    
    # Perform scan
    println("Calculating t_inf_k0")
    t_inf_k0 = calc_transAmpl_inf.(SP.lattice_type, SP.N_sheets, SP.a, SP.L, SP.Delta_range, Ref(SP.e1))
    println("Calculating t_inf_k")
    t_inf_k  = calc_transAmpl_inf.(SP.lattice_type, SP.N_sheets, SP.a, SP.L, SP.Delta_range, Ref(SP.e1), SP.e1_label, SP.drive_type, SP.w0, SP.k_n)
    
    
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