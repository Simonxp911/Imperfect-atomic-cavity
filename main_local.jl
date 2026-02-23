

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
    EP = EP_a370
    
    # Number of sheets
    N_sheets = 2
    
    # Lattice spacing 
    a = NaN
    
    # Inter-sheet distance
    L_ratio = 1
    L = NaN
    
    # Filling fraction 
    ff = 1.0 - 0.0
    
    # Gaussian position distribution width
    pos_unc_ratio = 0.0
    pos_unc = NaN
        
    # Set radius of sheets (in units of a) and whether to cut of corners (making the sheet rounded)
    radius = 8.5
    cut_corners = true
    
    # Number of array instantiations to calculate
    N_inst = 100
    
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
    DrP = DrP_Gaussian(w0, "forward", AP.radius, AP.a, AP.array; w0_ratio=w0_ratio)
    
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
    
    N_sheets_specs      = (2, 5)
    L_ratio_specs       = [1]
    L_specs             = [NaN]
    ff_specs            = (0.9, 0.95, 2)
    pos_unc_ratio_specs = (0.035, 0.05, 2)
    
    return ScanPar(N_sheets_specs, L_ratio_specs, L_specs, ff_specs, pos_unc_ratio_specs)
end


function define_SP(scanParams)
    # N_sheets, L_ratio, L, ff, pos_unc_ratio = scanParams
    
    # Experimental parameters for lattice spacing = 532 nm
    EP = EP_a532
    
    # Number of sheets
    N_sheets = scanParams[1]
    
    # Lattice spacing 
    a = NaN
    
    # Inter-sheet distance
    L_ratio = scanParams[2]
    L = scanParams[3]
    
    # Filling fraction 
    ff = scanParams[4]
    
    # Gaussian position distribution width
    pos_unc_ratio = scanParams[5]
    pos_unc = NaN
        
    # Set radius of sheets (in units of a) and whether to cut of corners (making the sheet rounded)
    radius = 6.0
    cut_corners = true
    
    # Number of array instantiations to calculate
    N_inst = 100
    
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
    DrP = DrP_Gaussian(w0, "forward", AP.radius, AP.a, AP.array; w0_ratio=w0_ratio)
    
    # Set detection parameters
    DeP = DeP_IntensityDefault(AP.radius, AP.a, AP.N_sheets, AP.L)
    
                     
    return SystemPar(EP,
                     L_ratio, pos_unc_ratio, AP,
                     Delta_specs,
                     k_n,
                     w0_ratio, DrP,
                     DeP)
end




function main()
    # Define system parameters
    SP  = define_SP()
    ScP = define_ScP()
    show(SP)
    # show(ScP)
    
    
    # Make figures
    # scan_TRCoef_fin(SP)
    # make_Tscan_fig(SP)
    # make_Tscan_comparison_fig(ScP)
    # make_Efield_intensity_fig(SP)
    # make_Efield_intensity_3D_fig(SP)
        
    return nothing
end


# ================================================
#   Generate figures
# ================================================
function make_Tscan_fig(SP)
    # Perform the scan
    Tscan, Rscan = scan_TRCoef_fin(SP)
    
    # Do statistics on Tscan
    T_means, T_stds = scan_statistics(Tscan)
    R_means, R_stds = scan_statistics(Rscan)
    
    # Get infinite system, normal-incidence, plane-wave transmission and with the chosen drive
    # t_inf_k0, t_inf_k = scan_transAmpl_inf(SP)
    # T_inf_k0 = abs2.(t_inf_k0)
    # T_inf_k  = abs2.(t_inf_k)
    # T_inf_k0 = false
    # T_inf_k  = false
    
    
    # Write a title and plot the Tscan
    # fig_Delta_scan(SP.Delta_range, Tscan, SP)
    # fig_Delta_scan_stats(SP.Delta_range, T_means, T_stds, T_inf_k0, T_inf_k, SP)
    fig_Delta_TRscan_stats(SP.Delta_range, T_means, T_stds, R_means, R_stds, SP)
end


function make_Tscan_comparison_fig(ScP)
    # Collect pre-calculated scans and perform statistics on each of them
    scanProd = scanProduct(ScP)
    SP_rep = nothing
    T_means = Array{Any}(undef, size(scanProd))
    T_stds  = deepcopy(T_means)
    R_means = deepcopy(T_means)
    R_stds  = deepcopy(T_means)
    for (i, scanParams) in enumerate(scanProd)
        SP = define_SP(scanParams)
        if i == 1 SP_rep = SP end
        
        Tscan, Rscan = scan_TRCoef_fin(SP)
        T_means[i], T_stds[i] = scan_statistics(Tscan)
        R_means[i], R_stds[i] = scan_statistics(Rscan)
    end
    
    titl, labels = scanTitlAndLabels(scanProd)
    titl = "lattice_type, radius, cc, N_inst = $(SP_rep.AP.lattice_type), $(SP_rep.AP.radius), $(SP_rep.AP.cut_corners), $(SP_rep.AP.N_inst) \n" *
           "drive, w0_ratio, dipoleMoment, detec_type = $(SP_rep.DrP.drive_type), $(SP_rep.w0_ratio), $(SP_rep.EP.dipoleMoment_label), $(SP_rep.DeP.detec_type) \n" *
           titl
    
    # Choose which slice of the data to plot
    # for i in 1:4
        slice_raw = (1:4, :, :, 2:2, 2:2)
        slice = CartesianIndices(to_indices(labels, slice_raw))
        fig_scanStatsComparison(SP_rep.Delta_range, T_means[slice], T_stds[slice], R_means[slice], R_stds[slice], labels[slice], titl)
    # end
end


function make_Efield_intensity_fig(SP)
    # Get collective energies and choose the detuning of perfect transmission
    Gk = ana_FT_GF(SP.AP.lattice_type, SP.AP.a, SP.EP.dipoleMoment, SP.EP.dipoleMoment)
    tildeDelta = -real(Gk)
    tildeGamma =  imag(Gk)
    Δ = tildeDelta - tildeGamma*tan(ωa*SP.AP.L)
    
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
    Δ = tildeDelta - tildeGamma*tan(ωa*SP.AP.L)
    
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

# Ask David how they normalize their transmission 

# Consider cleaning up code (comments, names, )

# Runs:
    # See David's summary
    # Make larger comparison of N_sheets = 2, 3, 4, 5
    # Find out when the bump/peak disappears (as a function of ff and pos_unc) for L = 2.05
    