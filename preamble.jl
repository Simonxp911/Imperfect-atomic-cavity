

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
#   Structures and constants
# ================================================
struct ExperimentalPar
    latticeSpacing_inPlane::Real                # [nm] Lattice spacing within each sheet
    latticeSpacing_outOfPlane::Real             # [nm] Lattice spacing between sheets (i.e. the inter-sheet distance is a multiple of this)
    transitionWavelength::Real                  # [nm] Wavelength of the atomic transition
    excitedStateDecayRate::Real                 # [MHz] Excited state decay rate    
    dipoleMoment::Vector                                  # Dipole moment vector
    dipoleMoment_label::String                            # Dipole moment vector label
    
    function ExperimentalPar(latticeSpacing_inPlane::Real, latticeSpacing_outOfPlane::Real, 
                             transitionWavelength::Real, excitedStateDecayRate::Real,
                             dipoleMoment, dipoleMoment_label)
        
        dipoleMoment /= norm(dipoleMoment)
        
        return new(latticeSpacing_inPlane, latticeSpacing_outOfPlane,
                   transitionWavelength, excitedStateDecayRate,
                   dipoleMoment, dipoleMoment_label)
    end
end


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


struct DrivePar
    drive_type::String                          # Which type of driving 
    w0::Real                                    # Driving width (or beam waist)
    drivemode::Vector                           # Vector holding the actual drive amplitudes
    
    
    function DrivePar(drive_type::String, w0::Real, array::Vector)
        drivemode = prepare_drivemode.(drive_type, array, w0)
        return new(drive_type, w0, drivemode)
    end
    
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


struct DetectionPar
    detec_type::String                          # Which type of detection
    detec_z::Real                               # The distance at which the integration/detection planes are positioned
    integration_plane_n::Int                    # Resolution for the integration plane
    integration_plane_radius::Real              # Radius of the integration plane
    integration_plane::Vector                   # Vector of position vectors of the detection plane
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


struct SystemPar
    EP::ExperimentalPar                                  # Experimental parameters struct
    
    L_ratio::Int                                # Inter-sheet distance in units EP.latticeSpacing_outOfPlane
    pos_unc_ratio::Real                         # Position uncertainty in units of lattice spacing
    AP::ArrayPar                                # Array parameters struct
    
    Delta_specs::Tuple{Real, Real, Int}         # Specs for detuning
    Delta_range::AbstractRange                  # Range of detuning
    
    k_n::Int                                    # Number of k-space points along the positive first axis (for integration over the first octant of the BZ, when calculating transmission of finite beam on infinite array)
    
    w0_ratio::Real                              # Driving width (or beam waist) in units of the array radius
    DrP::DrivePar                               # Driving parameters struct
    
    DeP::DetectionPar                            # Detection parameters struct
    
    
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


struct ScanPar
    N_sheets_specs::Tuple{Real, Real, Int}      # Specs for number of atomic sheets
    L_ratio_specs::Tuple{Real, Real, Int}       # Specs for L_ratio
    L_specs::Tuple{Real, Real, Int}             # Specs for inter-sheet distance
    ff_specs::Tuple{Real, Real, Int}            # Specs for filling fraction
    pos_unc_ratio_specs::Tuple{Real, Real, Int} # Specs for position uncertainty ratio

    N_sheets_range::AbstractRange               # Range of number of atomic sheets
    L_ratio_range::AbstractRange                # Range of L_ratio
    L_range::AbstractRange                      # Range of inter-sheet distance
    ff_range::AbstractRange                     # Range of filling fraction
    pos_unc_ratio_range::AbstractRange          # Range of position uncertainty ratio
    
    
    function ScanPar(N_sheets_specs::Tuple{Real, Real, Int}, 
                     L_ratio_specs::Tuple{Real, Real, Int},
                     L_specs::Tuple{Real, Real, Int},
                     ff_specs::Tuple{Real, Real, Int},
                     pos_unc_ratio_specs::Tuple{Real, Real, Int})
        
        N_sheets_range      = range(N_sheets_specs...)
        L_ratio_range       = range(L_ratio_specs...)
        L_range             = range(L_specs...)
        ff_range            = range(ff_specs...)
        pos_unc_ratio_range = range(pos_unc_ratio_specs...)
                     
        return new(N_sheets_specs, L_ratio_specs, L_specs, ff_specs, pos_unc_ratio_specs,
                   N_sheets_range, L_ratio_range, L_range, ff_range, pos_unc_ratio_range)
    end
end


function Base.show(io::IO, EP::ExperimentalPar)
    println(io, "--- Experimental Parameters ---")
    
    println(io, "Dimensionfull experimental parameters")
    println(io, "latticeSpacing_inPlane: ", EP.latticeSpacing_inPlane)
    println(io, "latticeSpacing_outOfPlane: ", EP.latticeSpacing_outOfPlane)
    println(io, "transitionWavelength: ", EP.transitionWavelength)
    println(io, "excitedStateDecayRate: ", EP.excitedStateDecayRate)
    println(io, "dipoleMoment: ", EP.dipoleMoment)
    println(io, "dipoleMoment_label: ", EP.dipoleMoment_label)
    println(io, "")
    
    println(io, "---  ---")
end


function Base.show(io::IO, AP::ArrayPar)
    println(io, "--- Array Parameters ---")
    
    println(io, "Dimensionfull experimental parameters")
    println(io, "lattice_type: ", AP.lattice_type)
    println(io, "N_sheets: ", AP.N_sheets)
    println(io, "a: ", AP.a)
    println(io, "L: ", AP.L)
    println(io, "ff: ", AP.ff)
    println(io, "pos_unc: ", AP.pos_unc)
    println(io, "radius: ", AP.radius)
    println(io, "cut_corners: ", AP.cut_corners)
    println(io, "N_inst: ", AP.N_inst)
    println(io, "N: ", AP.N)
    println(io, "")
    
    println(io, "---  ---")
end


function Base.show(io::IO, DrP::DrivePar)
    println(io, "--- Drive Parameters ---")
    
    println(io, "Dimensionfull experimental parameters")
    println(io, "drive_type: ", DrP.drive_type)
    println(io, "w0: ", DrP.w0)
    println(io, "drivemode", DrP.drivemode)
    println(io, "")
    
    println(io, "---  ---")
end


function Base.show(io::IO, DeP::DetectionPar)
    println(io, "--- Detection Parameters ---")
    
    println(io, "Dimensionfull experimental parameters")
    println(io, "detec_type: ", DeP.detec_type)
    println(io, "detec_z: ", DeP.detec_z)
    println(io, "integration_plane_n: ", DeP.integration_plane_n)
    println(io, "integration_plane_radius: ", DeP.integration_plane_radius)
    println(io, "dx: ", DeP.dx)
    println(io, "dy: ", DeP.dy)
    println(io, "detec_radius: ", DeP.detec_radius)
    println(io, "")
    
    println(io, "---  ---")
end


function Base.show(io::IO, SP::SystemPar)
    println(io, "--- System Parameters ---")
    
    println(io, "")
    show(SP.EP)
    println(io, "")
    
    println(io, "L_ratio: ", SP.L_ratio)
    println(io, "pos_unc_ratio: ", SP.pos_unc_ratio)
    
    println(io, "")
    show(SP.AP)
    println(io, "")
    
    println(io, "Delta_specs: ", SP.Delta_specs)
    println(io, "k_n: ", SP.k_n)
    println(io, "w0_ratio: ", SP.w0_ratio)
    
    println(io, "")
    show(SP.DrP)
    println(io, "")
    
    println(io, "")
    show(SP.DeP)
    println(io, "")
    
    println(io, "---  ---")
end


function Base.show(io::IO, ScP::ScanPar)
    println(io, "--- Scan Parameters ---")

    println(io, "Dimensionfull experimental parameters")
    println(io, "N_sheets_specs: ", ScP.N_sheets_specs)
    println(io, "L_ratio_specs: ", ScP.L_ratio_specs)
    println(io, "L_specs: ", ScP.L_specs)
    println(io, "ff_specs: ", ScP.ff_specs)
    println(io, "pos_unc_ratio_specs: ", ScP.pos_unc_ratio_specs)
    println(io, "")
    
    println(io, "---  ---")
end


function Base.show(io::IO, ScP::ScanPar)
    println(io, "--- Scan Parameters ---")
    
    println(io, "Range specs")
    println(io, "N_sheets_specs: ", ScP.N_sheets_specs)
    println(io, "L_ratio_specs: ", ScP.L_ratio_specs)
    println(io, "L_specs: ", ScP.L_specs)
    println(io, "ff_specs: ", ScP.ff_specs)
    println(io, "pos_unc_ratio_specs: ", ScP.pos_unc_ratio_specs)
    println(io, "")
    
    println(io, "---  ---")
end

const wa = 2π
const EP_a532 = ExperimentalPar(532, 532, 780, 2π*6.065, [1, 1im, 0], "rc")
const EP_a370 = ExperimentalPar(370, 532, 780, 2π*6.065, [1, 1im, 0], "rc")


