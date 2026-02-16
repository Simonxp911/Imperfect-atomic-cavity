

# ================================================
#   Functions related to the real space atomic array
# ================================================
function get_array(lattice_type, N_sheets, a, L, radius, ff, pos_unc, N_inst, cut_corners=true)
    # Get the perfect array (no missing atoms, no randomness in position)
    perfect_array = get_perfect_array(lattice_type, N_sheets, a, L, radius, cut_corners)
    
    # Remove atoms to reach desired filling fraction
    holed_array = remove_atoms_from_array(perfect_array, ff, N_inst)
    
    # Introduce randomness to the atomic positions
    uncertain_holed_array = introduce_position_uncertainty_to_array_sites(holed_array, pos_unc)
    
    return uncertain_holed_array
end


function get_perfect_array(lattice_type, N_sheets, a, L, radius, cut_corners=true)
    # Get lattice vectors
    lattice_vectors = get_lattice_vectors(lattice_type, a)
    
    # Construct single sheet, cutting off corners
    sheet = get_single_sheet(lattice_vectors, a, radius, cut_corners)
    
    # Repeat sheet N_sheets times
    array = Vector{Vector{Float64}}(undef, 0)
    for i in 0:N_sheets-1
        z_shift = i*L - (N_sheets - 1)/2*L
        push!.(Ref(array), sheet .+ Ref([0, 0, z_shift]))
    end
    
    # return array[sortperm(norm.(array))]
    return array
end


function get_single_sheet(lattice_vectors, a, radius, cut_corners)
    # To ensure that all sites within the radius are included, an oversize sheet 
    # with radius loop_radius = 1.5*radius is constructed, and its corners cut off
    if cut_corners
        loop_radius = ceil(1.5*radius)
    else
        loop_radius = radius
    end
    
    # Construc sheet
    sheet = Vector{Vector{Float64}}(undef, 0)
    for site_indices in Iterators.product((-loop_radius:loop_radius for x in eachindex(lattice_vectors))...)
        site_vec = sum(site_indices .* lattice_vectors)
        if cut_corners
            if norm(site_vec) <= a*radius + 1e-3
                push!(sheet, site_vec)
            end
        else
            push!(sheet, site_vec)
        end
    end
    return sheet
end


function get_lattice_vectors(lattice_type, a)
    if lattice_type == "square"
        return a*[1, 0, 0], a*[0, 1, 0]
    elseif lattice_type == "triangular"
        return a*[1, 0, 0], a*[1/2, sqrt(3)/2, 0]
    elseif lattice_type == "linear"
        return (a*[1, 0, 0],)
    end
end


function remove_atoms_from_array(array, ff, N_inst)
    if ff == 1
        return fill(deepcopy(array), N_inst)
    end
    
    # Find total number of atoms and the number of atoms to be kept to match the desired filling fraction
    N = length(array)
    N_to_be_kept = Int(floor(N*ff))
    
    # Return N_to_be_kept of the original array sites
    return [array[randperm(N)[1:N_to_be_kept]] for _ in 1:N_inst]
end


function introduce_position_uncertainty_to_array_sites(array, pos_unc)
    if pos_unc == 0
        return array
    end
    
    # Generate normally-distributed random numbers for each coordinate of each atom
    random_shift = [pos_unc*randn.(fill(3, N)) for N in length.(array)]
    
    # Return the randomly shifted array sites
    return array + random_shift
end


# ================================================
#   Functions related to driving
# ================================================
function get_drivemode(drive_type, r, w0)
    # Get the driving mode evaluated at r
    if drive_type == "homogenous"
        return exp(1im*wa*r[3]) 
    elseif drive_type == "Gaussian"
        return Gaussian(w0, r, true)
    else
        throw(DomainError(drive_type, "This drive_type has not been implemented in get_drivemode"))
    end
end


function Gaussian(w0, rvec, normalize_or_not)
    # Extract coordinates (squared transverse distance and z)
    r2 = rvec[1]^2 + rvec[2]^2
    z  = rvec[3]
    
    # Calculate components
    zR = π*w0^2
    wz = w0*sqrt(1 + (z/zR)^2)
    prop_phase = wa*z
    tran_phase = wa*r2*z/(2*(zR^2 + z^2))
    Gouy_phase = atan(z/zR)
    
    # Put together components
    width_factor = w0/wz
    Gaussian     = exp(-r2/wz^2)
    phase        = exp(1im*(prop_phase + tran_phase - Gouy_phase))
    
    # The mode can be normalized or not
    if normalize_or_not
        N = π*w0^2/2
    else
        N = 1.0
    end
    
    # Finally return the full Gaussian mode
    return width_factor*Gaussian*phase/sqrt(N)
end


function get_drive_kspace(drive_type, w0, dk, kperp)
    if drive_type == "homogenous"
        return all(kperp .== 0.0) * sqrt((2π)^2/dk^2)
    elseif drive_type == "Gaussian"
        return Gaussian_k(w0, kperp, true)
    else
        throw(DomainError(drive_type, "This drive_type has not been implemented in get_drive_kspace"))
    end
end


function Gaussian_k(w0, kperp, normalize_or_not)
    # The mode can be normalized or not
    if normalize_or_not
        N = 1/(2π*w0^2)
    else
        N = 1.0
    end
    
    # Finally return the full Gaussian mode
    return exp(-w0^2*norm(kperp)^2/4.0)/sqrt(N)
end


# ================================================
#   Functions related to the Green's function
# ================================================
function get_Gnm(array, N, dipoleMoment)
    return realspace_GF.(reshape(array, N, 1), 
                         reshape(array, 1, N), 
                         Ref(dipoleMoment), Ref(dipoleMoment))
end


function get_Gmat_rn(r, array)
    return realspace_GF_matrix.(Ref(r), array)
end


function prepare_Σ(lattice_type, a, L, k_n, dipoleMoment, dipoleMoment_label)
    # Check if the Σ0 has already been calculated
    postfix = get_postfix_Sigma0(lattice_type, a, L, k_n, dipoleMoment_label)
    filename_Σ0 = "Sigma0_" * postfix
    data = check_if_already_calculated(save_dir, [filename_Σ0], ComplexF64)
    if length(data) == 1 return data[1] end
    
    # Set the k-space/BZ resolution
    dk = (π/a)/(k_n - 1)
    k_range = (0:k_n-1)*dk
    Σ0 = zeros(ComplexF64, k_n, k_n)
    for (i, kx) in enumerate(k_range), (j, ky) in enumerate(k_range)
        # Only consider the first octant
        if ky > kx continue end
        
        # Calculate and assign self-energy
        Σ0[i, j] = -ana_FT_GF(lattice_type, a, dipoleMoment, dipoleMoment, kx, ky, L, 1e-5)
        Σ0[j, i] = Σ0[i, j]
    end
    
    # Save Σ0
    save_as_txt(Σ0, save_dir, filename_Σ0)
    
    return Σ0
end


