

# ================================================
#   Functions related to packing/unpacking data
# ================================================
function pack_tscan(tscan)
    return tscan
end


function unpack_tscan(data)
    return data
end


function pack_tscan_inf(tscan_inf)
    return tscan_inf
end


function unpack_tscan_inf(data)
    return data
end


# ================================================
#   Functions related to filenames
# ================================================
function format_Complex_to_String(z)
    # Return a zero if z = 0
    if z == 0 return "0" end
    
    # Format the real and imaginary part
    r_string = real(z) != 0 ? @sprintf("%.3f", real(z))        : ""
    i_string = imag(z) != 0 ? @sprintf("%.3f", imag(z)) * "im" : ""
    
    # If the imaginary part is negative is has a sign, otherwise we add a plus, but only if the real part is nonzero
    if real(z) != 0 && imag(z) > 0 i_string = "+" * i_string end
    
    return r_string * i_string
end


function get_postfix(lattice_type, N_sheets, radius, cut_corners, a, L, ff, pos_unc_ratio, N_inst, drive_type, w0_ratio, e1_label, detec_mode, detec_radius, detec_z, Delta_specs)
    postfix_components = []
    
    push!(postfix_components, "_$(lattice_type[1:2])_Nsh_$(N_sheets)_r_$(radius)_cc_$(cut_corners)")
    push!(postfix_components, "_a_$(round(a, sigdigits=4))_L_$(round(L, sigdigits=4))")
    push!(postfix_components, "_ff_$(round(ff, sigdigits=3))_pur_$(round(pos_unc_ratio, sigdigits=3))_Ni_$(N_inst)")
    push!(postfix_components, "_dr_$(drive_type[1:4])_w0r_$(w0_ratio)")
    push!(postfix_components, "_e1_" * e1_label)
    push!(postfix_components, "_dm_$(detec_mode[1:5])_dr_$(round(detec_radius, sigdigits=4))_dz_$(round(detec_z, sigdigits=4))")
    push!(postfix_components, "_D_$(join(Delta_specs, ","))")
    
    return join(postfix_components)
end


function get_postfix(lattice_type, N_sheets, drive_type, w0_ratio, k_n, e1_label, Delta_specs)
    postfix_components = []
    
    push!(postfix_components, "_$(lattice_type[1:2])_Nsh_$(N_sheets)")
    push!(postfix_components, "_dr_$(drive_type[1:4])_w0r_$(w0_ratio)_kn_$(k_n)")
    push!(postfix_components, "_e1_" * e1_label)
    push!(postfix_components, "_D_$(join(Delta_specs, ","))")
    
    return join(postfix_components)
end


function get_postfix(lattice_type, a, L, k_n, e1_label)
    postfix_components = []
    
    push!(postfix_components, "_$(lattice_type[1:2])")
    push!(postfix_components, "_a_$(round(a, sigdigits=4))_L_$(round(L, sigdigits=4))")
    push!(postfix_components, "_kn_$(k_n)_e1_" * e1_label)
    
    return join(postfix_components)
end


# ================================================
#   Functions related to saving and loading
# ================================================
function check_if_already_calculated(save_dir, filenames, entry_type=Float64)
    data = []
    for filename in filenames
        if isfile(save_dir * filename * ".jld2")
            push!(data, load_as_jld2(save_dir, filename))
        elseif isfile(save_dir * filename * ".txt")
            push!(data, load_as_txt(save_dir, filename, entry_type))
        end
    end
    return tuple(data...)
end


function save_as_jld2(data, save_dir, filename)
    JLD2.save(save_dir * filename * ".jld2", "data", data)
end


function load_as_jld2(save_dir, filename)
    return load(save_dir * filename * ".jld2", "data")
end


function save_as_txt(data, save_dir, filename)
    writedlm(save_dir * filename * ".txt", data, '\t')
end


function load_as_txt(save_dir, filename, entry_type=Float64)
    return readdlm(save_dir * filename * ".txt", '\t', entry_type, '\n')
end
