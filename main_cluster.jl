
using LinearAlgebra #norm of vectors and other standard linear algebra
using JLD2 #saving and loading
using DelimitedFiles #read/write simple text data files
using Random #for randomly making imperfect lattices
using Printf #for formatting strings
using Statistics #for calculating mean, standard deviation, etc.

const wa = 2π
const save_dir = ARGS[2] * "/imperfect_atomic_cavity_data/"

include("calcs.jl")
include("utility.jl")
include("calc_Greens.jl")
include("save_load.jl")


# ================================================
#   Main functions
# ================================================
function define_system_parameters()
    # Read input file
    input = readlines(ARGS[1])
    
    # Slurm task ID
    task_id = parse(Int, ARGS[3])
    printlnX("TASK_ID = $task_id")
    
    
    # Experiment-specific parameters
    a_dimensionfull = 532                #nm  (lattice spacing)
    # λ_dimensionfull = 780                #nm  (transition wavelength)
    λ_dimensionfull = 795                #nm  (transition wavelength)
    γ_dimensionfull = 2π*6.065           #MHz (excited state decay rate) 
    c_dimensionfull = 299792458*1e9*1e-6 #nm*MHz (speed of light)
    L_ratio                     = read_input(input, "L_ratio", "Float64")
    L_dimensionfull = L_ratio*a_dimensionfull  #nm  (inter-sheet distance = lattice spacing times any integer)
    
    # Dimesionless equivalents
    a = a_dimensionfull/λ_dimensionfull #0.682
    L = L_dimensionfull/λ_dimensionfull #2.046
    
    # Read L directly from the input
    L                           = read_input(input, "L", "Float64")
    
    # Specifications for ranges of parameters
    Delta_specs                 = read_input(input, "Delta_specs", "TupleFlFlInt")
    ff_specs                    = read_input(input, "ff_specs", "TupleFlFlInt")
    pos_unc_ratio_specs         = read_input(input, "pos_unc_ratio_specs", "TupleFlFlInt")
    
    # Define ranges
    Delta_range                 = range(Delta_specs...)
    ff_range                    = range(ff_specs...)
    pos_unc_ratio_range         = range(pos_unc_ratio_specs...)
    
    # Polarization vector of atomic transition
    e1                          = read_input(input, "e1", "VectorComplexF64")
    e1 /= norm(e1)
    e1_label                    = read_input(input, "e1_label", "String")
    
    # Define lattice_type ("square", "triangular", "linear")
    lattice_type                = read_input(input, "lattice_type", "String")
    
    # Number of sheets
    N_sheets                    = read_input(input, "N_sheets", "Int")
        
    # Filling fraction 
    # ff                          = read_input(input, "ff", "Float64")
    ff_ind = ((task_id .- 1) .% ff_specs[3]) .+ 1
    ff = ff_range[ff_ind]
    
    # Gaussian position distribution width
    # pos_unc_ratio               = read_input(input, "pos_unc_ratio", "Float64")
    pos_unc_ratio_ind = ((task_id .- 1) .÷ ff_specs[3]) .+ 1
    pos_unc_ratio = pos_unc_ratio_range[pos_unc_ratio_ind]
    pos_unc = pos_unc_ratio*a
        
    # Set radius of sheets and whether to cut of corners (making the sheet rounded)
    # With radius = 7 and cut_corners = true, we match the experiment
    radius                      = read_input(input, "radius", "Float64")
    cut_corners                 = read_input(input, "cut_corners", "Bool")
    
    # Number of array instantiations to calculate
    N_inst                      = read_input(input, "N_inst", "Int")
    
    # Get array and determine number of atoms
    array = get_array(lattice_type, N_sheets, a, L, radius, ff, pos_unc, N_inst, cut_corners)
    # array = [0.5*radius*a*randn.(fill(3, N)) for N in length.(array)]
    N = length(array[1])
    
    # Vector of interaction matrices
    Gnm = get_Gnm.(array, N, e1)
    
    # The driving is assumed to have polarization e1, as it would only be the corresponding component that contributed to the driving anyway
    # Set type of driving for the case of a finite array ("homogenous" [is not properly normalized for the finite case], "Gaussian")
    drive_type                  = read_input(input, "drive_type", "String")
    
    # Beam waist for Gaussian drive
    w0_ratio                    = read_input(input, "w0_ratio", "Float64")
    w0 = w0_ratio*radius*a
    
    # Number of k-space points along the positive first axis (for integration over the first octant of the BZ, when calculating transmission of finite beam on infinite array)
    k_n                         = read_input(input, "k_n", "Int") 
    
    # Vector of (normalized) driving amplitude at array sites 
    drivemode = prepare_drivemode.(drive_type, array, w0)
    
    # Set the detection mode when calculating the finite array transmission by direct integration ("drive_mode", "intensity_on_detection_plane")
    detec_mode                  = read_input(input, "detec_mode", "String")
    
    # Set the radius of the detection plane for detec_mode = "intensity_on_detection_plane"
    detec_radius = radius*a/2
    
    # Set z-position of detector for the case of detec_mode = "intensity_on_detection_plane"
    if N_sheets == 1
        detec_z = 1.0
    elseif N_sheets == 2
        detec_z = L/2 + 1.0
    else
        detec_z = 5.0
    end
    
    x_range = range(-2*SP.radius*SP.a, 2*SP.radius*SP.a, 31)
        y_range = deepcopy(x_range)
        integration_plane = [[x, y, SP.detec_z] for x in x_range, y in y_range]
        dx = x_range[2] - x_range[1]
        dy = y_range[2] - y_range[1]
    
    
    return (a_dimensionfull=a_dimensionfull, L_dimensionfull=L_dimensionfull,
            λ_dimensionfull=λ_dimensionfull, γ_dimensionfull=γ_dimensionfull,
            c_dimensionfull=c_dimensionfull,
            L_ratio=L_ratio,
            a=a, L=L,
            Delta_specs=Delta_specs, Delta_range=Delta_range,
            e1=e1, e1_label=e1_label,
            lattice_type=lattice_type, N_sheets=N_sheets, ff=ff, pos_unc_ratio=pos_unc_ratio, pos_unc=pos_unc, 
            radius=radius, cut_corners=cut_corners, N_inst=N_inst, array=array,
            N=N, Gnm=Gnm,
            drive_type=drive_type, w0_ratio=w0_ratio, w0=w0,
            k_n=k_n,
            drivemode=drivemode, detec_mode=detec_mode, detec_radius=detec_radius, detec_z=detec_z)
end


function main()
    # Define system parameters
    SP = define_system_parameters()
    print_SP(SP)
    
    # Make a scan of the transmission of the finite system
    scan_transmission_fin(SP)
        
    return nothing
end





printlnX("\n -- Running main() -- \n")
@time main()


