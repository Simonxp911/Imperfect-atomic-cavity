

# ================================================
#   Julia libraries
# ================================================
using LinearAlgebra         #norm of vectors and other standard linear algebra
BLAS.set_num_threads(1)     #needed to avoid LinearAlgebra becoming slow (BLAS.get_num_threads() = number of cores as default)
using JLD2                  #saving and loading
using DelimitedFiles        #read/write simple text data files
using Random                #for randomly making imperfect lattices
using Printf                #for formatting strings
using Statistics            #for calculating mean, standard deviation, etc.
using StatProfilerHTML      #profiling the code to see which parts take the most time to run


# ================================================
#   Files
# ================================================
include("calcs.jl")
include("utility.jl")
include("calc_Greens.jl")
include("save_load.jl")


# ================================================
#   Printing function
# ================================================
function showStruct(io::IO, structureTitle, structure, keys)
    longestKeyLength = maximum(length.(string.(keys)))
    
    println(io, "--- $structureTitle ---")
    for key in keys
        println(io, rpad("$key:", longestKeyLength + 3, " "), getfield(structure, key))
    end
    println(io, "---  ---")
end


# ================================================
#   Experimental parameters
# ================================================
struct ExperimentalPar
    latticeSpacing_inPlane::Real                # [nm] Lattice spacing within each sheet
    latticeSpacing_outOfPlane::Real             # [nm] Lattice spacing between sheets (i.e. the inter-sheet distance is a multiple of this)
    transitionWavelength::Real                  # [nm] Wavelength of the atomic transition
    excitedStateDecayRate::Real                 # [MHz] Excited state decay rate    
    dipoleMoment::Vector                        # Dipole moment vector
    dipoleMoment_label::String                  # Dipole moment vector label
    
    function ExperimentalPar(latticeSpacing_inPlane::Real, latticeSpacing_outOfPlane::Real, 
                             transitionWavelength::Real, excitedStateDecayRate::Real,
                             dipoleMoment, dipoleMoment_label)
        
        dipoleMoment /= norm(dipoleMoment)
        
        return new(latticeSpacing_inPlane, latticeSpacing_outOfPlane,
                   transitionWavelength, excitedStateDecayRate,
                   dipoleMoment, dipoleMoment_label)
    end
end


function Base.show(io::IO, EP::ExperimentalPar)
    keys = (:latticeSpacing_inPlane, :latticeSpacing_outOfPlane, :transitionWavelength, :excitedStateDecayRate, :dipoleMoment, :dipoleMoment_label)
    showStruct(io, "Experimental Parameters", EP, keys)
end


# ================================================
#   Array parameters
# ================================================
struct ArrayPar
    lattice_type::String                        # Which type of lattice 
    N_sheets::Int                               # Number of atomic sheets in the full array
    a::Real                                     # Dimensionless (in-plane) lattice spacing
    L::Real                                     # Dimensionless inter-lattice spacing
    ff::Real                                    # Filling fraction
    pos_unc::Real                               # Position uncertainty
    radius::Real                                # Radius of sheets in units of lattice spacing (i.e. number of atoms from the central atom to the sheet edge)
    cut_corners::Bool                           # Whether to cut corners of the array (to make it circular)
    N_inst::Int                                 # Number of instantiations (for array with random elements due to imperfect filling and position uncertainty)
    array::Vector                               # Vector of (N_inst) arrays
    N::Int                                      # Number of atoms in each array
    
    
    function ArrayPar(lattice_type::String, N_sheets::Int,
                      a::Real, L::Real,
                      ff::Real, pos_unc::Real,
                      radius::Real, cut_corners::Bool,
                      N_inst::Int)
        
        array = get_array(lattice_type, N_sheets, a, L, radius, ff, pos_unc, N_inst, cut_corners)
        N = length(array[1])
        
        return new(lattice_type, N_sheets,
                   a, L, ff, pos_unc,
                   radius, cut_corners,
                   N_inst,
                   array, N)
    end
end


function Base.show(io::IO, AP::ArrayPar)
    keys = (:lattice_type, :N_sheets, :a, :L, :ff, :pos_unc, :radius, :cut_corners, :N_inst, :N)
    showStruct(io, "Array Parameters", AP, keys)
end


function AP_Square(N_sheets::Int,
                   a::Real, L::Real,
                   ff::Real, pos_unc::Real,
                   radius::Real, cut_corners::Bool,
                   N_inst::Int,
                   EP::ExperimentalPar; 
                   L_ratio=nothing, pos_unc_ratio=nothing)
    
    lattice_type = "square"
    if isnan(a) a = EP.latticeSpacing_inPlane/EP.transitionWavelength end
    if isnan(L) 
        if !isnothing(L_ratio) 
            L = L_ratio*EP.latticeSpacing_outOfPlane/EP.transitionWavelength 
        else
            throw(ArgumentError("AP_Square received L = NaN but no L_ratio!"))
        end
    end
    if isnan(pos_unc)
        if !isnothing(pos_unc_ratio) 
            pos_unc = pos_unc_ratio*a 
        else
            throw(ArgumentError("AP_Square received pos_unc = NaN but no pos_unc_ratio!"))
        end
    end
                 
    return ArrayPar(lattice_type, N_sheets,
                    a, L,
                    ff, pos_unc,
                    radius, cut_corners,
                    N_inst)
end


# ================================================
#   Drive parameters
# ================================================
struct DrivePar
    drive_type::String                          # Which type of driving 
    w0::Real                                    # Driving width (or beam waist)
    drivemode::Vector                           # Vector holding the actual drive amplitudes
    
    
    function DrivePar(drive_type::String, w0::Real, array::Vector)
        drivemode = prepare_drivemode.(drive_type, array, w0)
        return new(drive_type, w0, drivemode)
    end
    
end


function Base.show(io::IO, DrP::DrivePar)
    keys = (:drive_type, :w0)
    showStruct(io, "Drive Parameters", DrP, keys)
end


function DrP_Gaussian(w0::Real, radius::Real, a::Real, array::Vector; 
                      w0_ratio=nothing)
    
    drive_type = "Gaussian"
    if isnan(w0)
        if !isnothing(w0_ratio)
            w0 = w0_ratio*radius*a 
        else
            throw(ArgumentError("DrP_Gaussian received w0 = NaN but no w0_ratio"))
        end
    end
    
    return DrivePar(drive_type, w0, array)
end


# ================================================
#   Detection parameters
# ================================================
struct DetectionPar
    detec_type::String                          # Which type of detection
    detec_z::Real                               # The distance at which the integration/detection planes are positioned
    integration_plane_n::Int                    # Resolution for the integration plane
    integration_plane_radius::Real              # Radius of the integration plane
    integration_plane::Matrix                   # Vector of position vectors of the detection plane
    dx::Real                                    # x difference for integration plane
    dy::Real                                    # y difference for integration plane
    detec_radius::Real                          # Radius of the detection plane
    detection_plane::Vector                     # Vector of position vectors of the detection plane
    
    
    function DetectionPar(detec_type::String, detec_z::Real,
                          integration_plane_n::Int, integration_plane_radius::Real, 
                          detec_radius::Real)
        
        x_range = range(-integration_plane_radius, integration_plane_radius, integration_plane_n)
        integration_plane = [[x, y, detec_z] for x in x_range, y in x_range]
        dx = x_range[2] - x_range[1]
        dy = dx
        detection_plane = [r for r in integration_plane if r[1]^2 + r[2]^2 <= detec_radius^2]
        
        return new(detec_type, detec_z,
                   integration_plane_n, integration_plane_radius, integration_plane,
                   dx, dy,
                   detec_radius, detection_plane)
    end
end


function Base.show(io::IO, DeP::DetectionPar)
    keys = (:detec_type, :detec_z, :integration_plane_n, :integration_plane_radius, :dx, :dy, :detec_radius)
    showStruct(io, "Detection Parameters", DeP, keys)
end


function DeP_IntensityDefault(radius::Real, a::Real, N_sheets::Int, L::Real)
    
    detec_type = "intensity_on_detection_plane"
    detec_z = (N_sheets - 1)*L/2 + 1.0
    integration_plane_n = 51
    integration_plane_radius = 2*radius*a
    detec_radius = radius*a/2

    return DetectionPar(detec_type, detec_z,
                        integration_plane_n, integration_plane_radius, 
                        detec_radius)
end


# ================================================
#   System parameters
# ================================================
struct SystemPar
    EP::ExperimentalPar                         # Experimental parameters struct
    
    L_ratio::Int                                # Inter-sheet distance in units EP.latticeSpacing_outOfPlane
    pos_unc_ratio::Real                         # Position uncertainty in units of lattice spacing
    AP::ArrayPar                                # Array parameters struct
    
    Delta_specs::Tuple{Real, Real, Int}         # Specs for detuning
    Delta_range::AbstractRange                  # Range of detuning
    
    k_n::Int                                    # Number of k-space points along the positive first axis (for integration over the first octant of the BZ, when calculating transmission of finite beam on infinite array)
    
    w0_ratio::Real                              # Driving width (or beam waist) in units of the array radius
    DrP::DrivePar                               # Driving parameters struct
    
    DeP::DetectionPar                           # Detection parameters struct
    
    
    function SystemPar(EP::ExperimentalPar,
                       L_ratio::Real, pos_unc_ratio::Real, AP::ArrayPar,
                       Delta_specs::Tuple{Real, Real, Int},
                       k_n::Int,
                       w0_ratio::Real, DrP::DrivePar,
                       DeP::DetectionPar)
        
        Delta_range = range(Delta_specs...)
        
        return new(EP,
                   L_ratio, pos_unc_ratio, AP,
                   Delta_specs, Delta_range,
                   k_n,
                   w0_ratio, DrP,
                   DeP)
    end
end


function Base.show(io::IO, SP::SystemPar)
    println(io, "--- System Parameters ---")
    println(io, "")
    show(SP.EP)
    println(io, "")
    println(io, rpad("L_ratio:", 16, " "), SP.L_ratio)
    println(io, rpad("pos_unc_ratio:", 16, " "), SP.pos_unc_ratio)
    println(io, "")
    show(SP.AP)
    println(io, "")
    println(io, rpad("Delta_specs:", 16, " "), SP.Delta_specs)
    println(io, rpad("k_n: :", 16, " "), SP.k_n)
    println(io, rpad("w0_ratio: :", 16, " "), SP.w0_ratio)
    println(io, "")
    show(SP.DrP)
    println(io, "")
    show(SP.DeP)
    println(io, "")
    println(io, "---  ---")
    println(io, "")
end


# ================================================
#   Scan parameters
# ================================================
struct ScanPar
    N_sheets_range::AbstractVector              # Range of number of atomic sheets
    L_ratio_range::AbstractVector               # Range of L_ratio
    L_range::AbstractVector                     # Range of inter-sheet distance
    ff_range::AbstractVector                    # Range of filling fraction
    pos_unc_ratio_range::AbstractVector         # Range of position uncertainty ratio
    
    
    function ScanPar(N_sheets_specs::Union{Tuple{Int, Int}, AbstractVector}, 
                     L_ratio_specs::Union{Tuple{Int, Int}, AbstractVector},
                     L_specs::Union{Tuple{Real, Real, Int}, AbstractVector},
                     ff_specs::Union{Tuple{Real, Real, Int}, AbstractVector},
                     pos_unc_ratio_specs::Union{Tuple{Real, Real, Int}, AbstractVector})
        
        if typeof(N_sheets_specs) <: Tuple{Int, Int}
            N_sheets_range = N_sheets_specs[1]:N_sheets_specs[2]
        else
            N_sheets_range = N_sheets_specs
        end
        
        if typeof(L_ratio_specs) <: Tuple{Int, Int}
            L_ratio_range = L_ratio_specs[1]:L_ratio_specs[2]
        else
            L_ratio_range = L_ratio_specs
        end
        
        if typeof(L_specs) <: Tuple{Real, Real, Int}
            L_range = range(L_specs...)
        else
            L_range = L_specs
        end
        
        if typeof(ff_specs) <: Tuple{Real, Real, Int}
            ff_range = range(ff_specs...)
        else
            ff_range = ff_specs
        end
        
        if typeof(pos_unc_ratio_specs) <: Tuple{Real, Real, Int}
            pos_unc_ratio_range = range(pos_unc_ratio_specs...)
        else
            pos_unc_ratio_range = pos_unc_ratio_specs
        end
        
        return new(N_sheets_range, L_ratio_range, L_range, ff_range, pos_unc_ratio_range)
    end
end


function Base.show(io::IO, ScP::ScanPar)
    keys = (:N_sheets_range, :L_ratio_range, :L_range, :ff_range, :pos_unc_ratio_range)
    showStruct(io, "Scan Parameters", ScP, keys)
end


# ================================================
#   Constants
# ================================================
const wa = 2π
const EP_a532 = ExperimentalPar(532, 532, 780, 2π*6.065, [1, 1im, 0], "rc")
const EP_a370 = ExperimentalPar(370, 532, 780, 2π*6.065, [1, 1im, 0], "rc")





