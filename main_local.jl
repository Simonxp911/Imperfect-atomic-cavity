
using LinearAlgebra #norm of vectors and other standard linear algebra
using JLD2 #saving and loading
using DelimitedFiles #read/write simple text data files
using Plots; pythonplot() #plot using Python-Matplotlib as backend
using Colors #for generating distinguishable colors
using LaTeXStrings #LaTeX formatting in string in plots
using Random #for randomly making imperfect lattices
using Printf #for formatting strings
using Statistics #for calculating mean, standard deviation, etc.

const wa = 2π
const save_dir = "C:/Users/Simon/Forskning/Data/imperfect_atomic_cavity_data/"

include("calcs.jl")
include("utility.jl")
include("calc_Greens.jl")
include("figures.jl")
include("save_load.jl")


# ================================================
#   Main functions
# ================================================
function define_system_parameters()
    # Experiment-specific parameters
    a_dimensionfull = 532                #nm  (lattice spacing)
    # a_dimensionfull = 370                #nm  (lattice spacing)
    # λ_dimensionfull = 780                #nm  (transition wavelength)
    λ_dimensionfull = 795                #nm  (transition wavelength)
    γ_dimensionfull = 2π*6.065           #MHz (excited state decay rate) 
    c_dimensionfull = 299792458*1e9*1e-6 #nm*MHz (speed of light)
    L_ratio = 3.0
    L_dimensionfull = L_ratio*a_dimensionfull  #nm  (inter-sheet distance = lattice spacing times any integer)
    
    # Dimesionless equivalents
    a = a_dimensionfull/λ_dimensionfull #0.682
    L = L_dimensionfull/λ_dimensionfull #2.046 for L_ratio = 3.0
    
    # Set L manually
    L = 2.1
    
    # Specifications for ranges of parameters
    Delta_specs         = (-2.0, 3.0, 100)
    ff_specs            = (1.0, 1.0, 1)
    pos_unc_ratio_specs = (0.0, 0.0, 1)
    
    # Define ranges
    Delta_range         = range(Delta_specs...)
    ff_range            = range(ff_specs...)
    pos_unc_ratio_range = range(pos_unc_ratio_specs...)
    
    # Polarization vector of atomic transition
    e1 = [1, 1im, 0]/sqrt(2)
    e1_label = "rc"
    
    # Define lattice_type ("square", "triangular", "linear")
    lattice_type = "square"
    
    # Number of sheets
    N_sheets = 2
    
    # Filling fraction 
    ff = 1.0 - 0.0
    
    # Gaussian position distribution width
    pos_unc_ratio = 0.0
    pos_unc = pos_unc_ratio*a
        
    # Set radius of sheets (in units of a) and whether to cut of corners (making the sheet rounded)
    # With radius = 7 and cut_corners = true, we match the experiment (for a single sheet)
    # With radius = 6.0-6.5 and cut_corners = true, we match the experiment (for two sheets)
    radius = 6.0
    cut_corners = true
    
    # Number of array instantiations to calculate
    N_inst = 1
    
    # Get array and determine number of atoms
    array = get_array(lattice_type, N_sheets, a, L, radius, ff, pos_unc, N_inst, cut_corners)
    # array = [0.5*radius*a*randn.(fill(3, N)) for N in length.(array)]
    # array = [[]]
    N = length(array[1])
    # fig_array(array[1])
    # fig_array.(array)
    # println(N)
    
    # Set up interaction matrix
    Gnm = get_Gnm.(array, N, Ref(e1))
    
    # The driving is assumed to have polarization e1, as it would only be the corresponding component that contributed to the driving anyway
    # Set type of driving for the case of a finite array ("homogenous" [is not properly normalized for the finite case], "Gaussian")
    drive_type = "Gaussian"
    
    # Beam waist for Gaussian drive
    w0_ratio = 5.0
    w0       = w0_ratio*radius*a
    
    # Number of k-space points along the positive first axis (for integration over the first octant of the BZ, when calculating transmission of finite beam on infinite array)
    k_n = 100
    
    # Vector of (normalized) driving amplitude at array sites 
    drivemode = prepare_drivemode.(drive_type, array, w0)
    
    # Set the detection mode when calculating the finite array transmission by direct integration ("drive_mode", "integrated_drive_mode", "flat_mode_on_detection_plane", "incoming_mode_on_detection_plane", "intensity_on_detection_plane")
    detec_mode = "integrated_drive_mode"
    
    # Set the radius of the detection plane (for detec_mode = "flat_mode_on_detection_plane", "incoming_mode_on_detection_plane", "intensity_on_detection_plane")
    detec_radius = radius*a/2
    
    # Set z-position of detector (for detec_mode = "integrated_drive_mode", "flat_mode_on_detection_plane", "incoming_mode_on_detection_plane", "intensity_on_detection_plane")
    if N_sheets == 1
        detec_z = 1.0
    elseif N_sheets == 2
        detec_z = L/2 + 1.0
        # detec_z *= -1
    else
        detec_z = 5.0
    end
    
    # Set the the plane at which we calculate the transmission by integration
    x_range = range(-2*radius*a, 2*radius*a, 31)
    y_range = deepcopy(x_range)
    integration_plane = [[x, y, detec_z] for x in x_range, y in y_range]
    dx = x_range[2] - x_range[1]
    dy = y_range[2] - y_range[1]
    
    # Set the detection plane (for detec_mode = "flat_mode_on_detection_plane", "incoming_mode_on_detection_plane", "intensity_on_detection_plane")
    detection_plane = [r for r in integration_plane if r[1]^2 + r[2]^2 <= detec_radius^2]
    
    
    return (a_dimensionfull=a_dimensionfull, L_dimensionfull=L_dimensionfull,
            λ_dimensionfull=λ_dimensionfull, γ_dimensionfull=γ_dimensionfull,
            c_dimensionfull=c_dimensionfull,
            L_ratio=L_ratio,
            a=a, L=L,
            Delta_specs=Delta_specs, Delta_range=Delta_range,
            ff_specs=ff_specs, ff_range=ff_range,
            pos_unc_ratio_specs=pos_unc_ratio_specs, pos_unc_ratio_range=pos_unc_ratio_range,
            e1=e1, e1_label=e1_label,
            lattice_type=lattice_type, N_sheets=N_sheets, ff=ff, pos_unc_ratio=pos_unc_ratio, pos_unc=pos_unc, 
            radius=radius, cut_corners=cut_corners, N_inst=N_inst, array=array,
            N=N, Gnm=Gnm,
            drive_type=drive_type, w0_ratio=w0_ratio, w0=w0,
            k_n=k_n,
            drivemode=drivemode, detec_mode=detec_mode, detec_radius=detec_radius, detec_z=detec_z,
            integration_plane=integration_plane, dx=dx, dy=dy, detection_plane=detection_plane)
end


function main()
    # Define system parameters
    SP = define_system_parameters()
    
    
    # Make figures
    scan_transCoef_fin(SP)
    # make_Tscan_fig(SP)
    make_Tscan_comparison_fig(SP)
    # fig_Efield_intensity(SP)
    # fig_Efield_intensity_3D(SP)
        
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


function make_Tscan_comparison_fig(SP)
    # Collect pre-calculated scans and perform statistics on each of them
    means = Array{AbstractArray}(undef, SP.ff_specs[3], SP.pos_unc_ratio_specs[3])
    stds  = Array{AbstractArray}(undef, SP.ff_specs[3], SP.pos_unc_ratio_specs[3])
    missing_indices = []
    for (i, ff) in enumerate(SP.ff_range), (j, pos_unc_ratio) in enumerate(SP.pos_unc_ratio_range)
        # Check if the scan is available and load it
        postfix = get_postfix(SP.lattice_type, SP.N_sheets, SP.radius, SP.cut_corners, SP.a, SP.L, ff, pos_unc_ratio, SP.N_inst, SP.drive_type, SP.w0_ratio, SP.e1_label, SP.detec_mode, SP.detec_radius, SP.detec_z, SP.Delta_specs)
        filename_ts = "Tscan" * postfix
        
        data = check_if_already_calculated(save_dir, [filename_ts])
        if length(data) == 1 
            Tscan = data[1]
            means[i, j], stds[i, j] = scan_statistics(Tscan)
        else 
            push!(missing_indices, i + SP.ff_specs[3]*(j - 1))
            means[i, j], stds[i, j] = [false], [false]
        end
    end
    if length(missing_indices) > 0 println("make_Tscan_comparison_fig missing_indices ($(length(missing_indices))) : ", join(missing_indices, ",")) end
    
    # Make the figure
    # fig_Delta_scan_stats_comparison(SP.Delta_range, means, stds, SP)
    # fig_Delta_scan_stats_comparison_fixed_pos_unc(SP.Delta_range, means, stds, SP)
    fig_Delta_scan_stats_comparison_fixed_ff(SP.Delta_range, means, stds, SP)
end





println("\n -- Running main() -- \n")
@time main()



