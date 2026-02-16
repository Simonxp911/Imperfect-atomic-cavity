

const save_dir = ARGS[1] * "/imperfect_atomic_cavity_data/"

include("preamble.jl")

using MPI
MPI.Init()
const comm     = MPI.COMM_WORLD
const root     = 0
const myRank   = MPI.Comm_rank(comm)
const commSize = MPI.Comm_size(comm)

# Include the input file (which defines define_ScP() and define_SP())
include(ARGS[2])



# ================================================
#   Main functions
# ================================================
function main()
    # Define scan    
    ScP = define_ScP()
    MPI.Barrier(comm)
    if myRank == root show(ScP) end
    MPI.Barrier(comm)
    
    performScan(ScP)
        
    return nothing
end


function performScan(ScP)
    scanProd = scanProduct(ScP)
    totalNumberOfJobs = scanLength(ScP)
    myStartIndex, myEndIndex = myStartIndex_and_myEndIndex(totalNumberOfJobs)
    for index in myStartIndex:myEndIndex
        println("My rank is $myRank, I am working on index = $index")
        scanParams = collect(scanProd)[index]
        SP = define_SP(scanParams)
        scan_transCoef_fin(SP)
    end
end


function myStartIndex_and_myEndIndex(totalNumberOfJobs)
    # First divide the totalNumberOfJobs evenly among the processes, rounded down
    # Then divide the remainingJobs among the lowest-rank processes 
    # In case totalNumberOfJobs < commSize, some processes will have zero jobs to do
    
    jobsPerProcessFloor = totalNumberOfJobs ÷ commSize
    remainingJobs = totalNumberOfJobs - jobsPerProcessFloor*commSize
    if myRank < remainingJobs
        myNumberOfJobs = jobsPerProcessFloor + 1
        myStartIndex   = 1 + myRank*myNumberOfJobs
    elseif jobsPerProcessFloor > 0
        myNumberOfJobs = jobsPerProcessFloor
        myStartIndex   = 1 + remainingJobs + myRank*myNumberOfJobs
    else
        myNumberOfJobs = 0
        myStartIndex   = 0
    end
    myEndIndex = myStartIndex + myNumberOfJobs - 1
    return myStartIndex, myEndIndex
end 



println("\n -- Running main() -- \n")
@time main()