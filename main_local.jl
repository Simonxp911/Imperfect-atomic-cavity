

using GLMakie               #for plotting
using Colors                #for generating distinguishable colors
using LaTeXStrings          #LaTeX formatting in string in plots

const save_dir = "C:/Users/Simon/Forskning/Data/imperfect_atomic_cavity_data/"

include("preamble.jl")
include("figures.jl")


# ================================================
#   Main functions
# ================================================
function define_SP()
    # Experimental parameters for lattice spacing = 532 nm
    EP = EP_a532
    
    # Number of sheets
    N_sheets = 2
    
    # Lattice spacing 
    a = NaN
    
    # Inter-sheet distance
    L = 1.321
    L_ratio = 1
    
    # Filling fraction 
    ff = 1.0 - 0.05
    
    # Gaussian position distribution width
    pos_unc = NaN
    pos_unc_ratio = 0.05
        
    # Set radius of sheets (in units of a) and whether to cut of corners (making the sheet rounded)
    radius = 6.0
    cut_corners = true
    
    # Number of array instantiations to calculate
    N_inst = 2
    
    # Set array parameters
    AP = AP_Square(N_sheets, a, L, ff, pos_unc, radius, cut_corners, N_inst, EP; L_ratio=L_ratio, pos_unc_ratio=pos_unc_ratio)
    
    # Specifications for detuning range
    Delta_specs = (-2.0, 3.0, 500)
    
    # Number of k-space points along the positive first axis (for integration over the first octant of the BZ, when calculating transmission of finite beam on infinite array)
    k_n = 100
    
    # Beam waist for Gaussian drive
    w0 = NaN
    w0_ratio = 5.0
    
    # Set drive parameters
    DrP = DrP_Gaussian(w0, AP.radius, AP.a, AP.array; w0_ratio=w0_ratio)
    
    # Set detection parameters
    DeP = DeP_IntensityDefault(AP.radius, AP.a, AP.N_sheets, AP.L)
    
                     
    return SystemPar(EP,
                     L_ratio, pos_unc_ratio, AP,
                     Delta_specs,
                     k_n,
                     w0_ratio, DrP,
                     DeP)
end


function define_ScP()
    
    N_sheets_specs      = (2, 3)
    L_ratio_specs       = (1, 3)
    L_specs             = (1.0, 2.0, 2)
    ff_specs            = (0.95, 0.95, 1)
    pos_unc_ratio_specs = (0.05, 0.1, 2)
    
    return ScanPar(N_sheets_specs, L_ratio_specs, L_specs, ff_specs, pos_unc_ratio_specs)
end


function main()
    # Define system parameters
    SP  = define_SP()
    ScP = define_ScP()
    # show(SP)
    # show(ScP)
    
    
    # Make figures
    # scan_transCoef_fin(SP)
    make_Tscan_fig(SP)
    # make_Tscan_comparison_fig(SP, ScP)
    # make_Efield_intensity_fig(SP)
    # make_Efield_intensity_3D_fig(SP)
        
    return nothing
end


# ================================================
#   Generate figures
# ================================================
function make_Tscan_fig(SP)
    # Perform the scan
    Tscan = scan_transCoef_fin(SP)
    
    # Do statistics on Tscan
    T_means, T_stds = scan_statistics(Tscan)
    
    # Get infinite system, normal-incidence, plane-wave transmission and with the chosen drive
    # t_inf_k0, t_inf_k = scan_transAmpl_inf(SP)
    # T_inf_k0 = abs2.(t_inf_k0)
    # T_inf_k  = abs2.(t_inf_k)
    T_inf_k0 = false
    T_inf_k  = false
    
    
    # Write a title and plot the Tscan
    # fig_Delta_scan(SP.Delta_range, Tscan, SP)
    fig_Delta_scan_stats(SP.Delta_range, T_means, T_stds, T_inf_k0, T_inf_k, SP)
end


function make_Tscan_comparison_fig(SP, ScP)
    # Collect pre-calculated scans and perform statistics on each of them
    means = Array{AbstractArray}(undef, len(ScP.ff_range), len(ScP.pos_unc_ratio_range))
    stds  = Array{AbstractArray}(undef, len(ScP.ff_range), len(ScP.pos_unc_ratio_range))
    missing_indices = []
    for (i, ff) in enumerate(ScP.ff_range), (j, pos_unc_ratio) in enumerate(ScP.pos_unc_ratio_range)
        # Check if the scan is available and load it
        postfix = get_postfix_Tscan(SP.AP.lattice_type, SP.AP.N_sheets, SP.AP.radius, SP.AP.cut_corners, SP.AP.a, SP.AP.L, ff, pos_unc_ratio, SP.AP.N_inst, SP.DrP.drive_type, SP.w0_ratio, SP.EP.dipoleMoment_label, SP.DeP.detec_type, SP.DeP.detec_radius, SP.DeP.detec_z, SP.Delta_specs)
        filename_ts = "Tscan_" * postfix
        
        data = check_if_already_calculated(save_dir, [filename_ts])
        if length(data) == 1 
            Tscan = data[1]
            means[i, j], stds[i, j] = scan_statistics(Tscan)
        else 
            push!(missing_indices, (i, j))
            means[i, j], stds[i, j] = [false], [false]
        end
    end
    if length(missing_indices) > 0 println("make_Tscan_comparison_fig missing_indices ($(length(missing_indices))) : ", join(missing_indices, ", ")) end
    
    # Make the figure
    # fig_Delta_scan_stats_comparison(SP.Delta_range, means, stds, SP)
    # fig_Delta_scan_stats_comparison_fixed_pos_unc(SP.Delta_range, means, stds, SP)
    fig_Delta_scan_stats_comparison_fixed_ff(SP.Delta_range, means, stds, SP)
end


function make_Efield_intensity_fig(SP)
    # Get collective energies and choose the detuning of perfect transmission
    Gk = ana_FT_GF(SP.AP.lattice_type, SP.AP.a, SP.EP.dipoleMoment, SP.EP.dipoleMoment)
    tildeDelta = -real(Gk)
    tildeGamma =  imag(Gk)
    Δ = tildeDelta - tildeGamma*tan(wa*SP.AP.L)
    
    # Find the steady state coherences
    Gnm = get_Gnm.(SP.AP.array, SP.AP.N, Ref(SP.EP.dipoleMoment))
    σ_ss = calc_σ_ss(Δ, Gnm[1], SP.DrP.drivemode[1])
    
    # Define x, y, and z ranges for the plot
    n = 101
    x_range = range(-3*SP.AP.radius*SP.AP.a, 3*SP.AP.radius*SP.AP.a, n)
    y_range = deepcopy(x_range)
    z_range = range(-3*SP.AP.L, 3*SP.AP.L, n)
    
    # Calculate the E-field intensity (in the xz and the zy planes)
    intensity_xz = zeros(length(x_range), length(z_range))
    for (i, x) in enumerate(x_range), (j, z) in enumerate(z_range)
        r = [x, 0.0, z]
        
        # We calculate E-field multiplied by d and divided by incoming amplitude
        Ed = calc_total_Efield_fin(r, SP.AP.array[1], σ_ss, SP.DrP.drive_type, SP.DrP.w0, SP.EP.dipoleMoment)
        
        intensity_xz[i, j] = Ed'*Ed
    end
    
    intensity_xy = zeros(length(x_range), length(y_range))
    for (i, x) in enumerate(x_range), (j, y) in enumerate(y_range)
        r = [x, y, maximum(z_range)]
        
        # We calculate E-field multiplied by d and divided by incoming amplitude
        Ed = calc_total_Efield_fin(r, SP.AP.array[1], σ_ss, SP.DrP.drive_type, SP.DrP.w0, SP.EP.dipoleMoment)
        
        intensity_xy[i, j] = Ed'*Ed
    end
    
    fig_Efield_intensity(x_range, y_range, z_range, intensity_xz, intensity_xy, SP.AP.array[1])
end


function make_Efield_intensity_3D_fig(SP)
    # Get collective energies and choose the detuning of perfect transmission
    Gk = ana_FT_GF(SP.AP.lattice_type, SP.AP.a, SP.EP.dipoleMoment, SP.EP.dipoleMoment)
    tildeDelta = -real(Gk)
    tildeGamma =  imag(Gk)
    Δ = tildeDelta - tildeGamma*tan(wa*SP.AP.L)
    
    # Find the steady state coherences
    Gnm = get_Gnm.(SP.AP.array, SP.AP.N, Ref(SP.EP.dipoleMoment))
    σ_ss = calc_σ_ss(Δ, Gnm[1], SP.DrP.drivemode[1])
    
    # Define x, y, and z ranges for the plot
    n = 101
    x_range = range(-3*SP.AP.radius*SP.AP.a, 3*SP.AP.radius*SP.AP.a, n)
    y_range = deepcopy(x_range)
    z_range = range(-5*SP.AP.L, 10*SP.AP.L, n)
    
    # Calculate the E-field intensity (in the xz, yz, and xy planes)
    intensities = [zeros(n, n) for i in 1:4]
    for i in 1:n, j in 1:n
        for (r, intensity) in zip(([x_range[i], 0.0, z_range[j]],
                                   [0.0, y_range[i], z_range[j]],
                                   [x_range[i], y_range[j], SP.DeP.detec_z],
                                   [x_range[i], y_range[j], z_range[end]]),
                                   intensities)
            # We calculate E-field multiplied by d and divided by incoming amplitude
            Ed = calc_total_Efield_fin(r, SP.AP.array[1], σ_ss, SP.DrP.drive_type, SP.DrP.w0, SP.EP.dipoleMoment)
            
            intensity[i, j] = Ed'*Ed
        end
    end
    
    fig_Efield_intensity_3D(x_range, y_range, z_range, intensities, SP.AP.array[1], SP.DeP.detec_z, SP.DeP.detec_radius)
end


println("\n -- Running main() -- \n")
@time main()



# TODO list:
# Implement phonon formalism
    # The appropriate functions can be copy-pasted from fiber_array
    # But because γ_a >> ν_α we would need to include phonons (?) and simulations would be limited in the number of atoms
# Read up on multiple layers (Shahmoon, Chang, Ruostekoski?)

# Change format of main_cluster, such that "input" files give define_SP
    # Use MPI
    # That will make it easier to do scans over any parameter
    # Consider making "subclasses" (subtypes? substructs?) to implement having a "standard" instantiation of the SP
    # and only having to define the parameters which are different from the "standard"
    # or only being forced to define some of the parameters while others take on standard values unless specifically changed
    # (Thus, define_SP would become much smaller to write and would exploit OOP or, in this case, being able to make subtypes

# Calculate reflection

# Implement comparison figures where L or N_sheets is the variable (instead of ff or pos_unc)
    # Plot both transmission and reflection

# Consider
    # Updating/cleaning up figures
    # Updating names (like Delta_ vs Δ_)

# Runs:
    # See David's summary
    # Make larger comparison of N_sheets = 2, 3, 4, 5
    # Find out when the bump/peak disappears (as a function of ff and pos_unc) for L = 2.05
    