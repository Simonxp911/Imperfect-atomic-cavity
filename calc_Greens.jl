

function ana_FT_GF(lattice_type, a, e1, e2, kx=0.0, ky=0.0, L=0.0, alpha=1e-7)
    if lattice_type ∈ ["square", "triangular"]
        return ana_FT_GF_2D(lattice_type, a, e1, e2, kx, ky, L, alpha)
    elseif lattice_type ∈ ["linear"]
        return ana_FT_GF_1D(a, e1, e2, kx)
    else
        throw(ArgumentError("lattice_type = $lattice_type has not been implemented in ana_FT_GF"))
    end
end

    
function ana_FT_GF_2D(lattice_type, a, e1, e2, kx=0.0, ky=0.0, L=0.0, alpha=1e-7)
    if L == 0
        return ana_FT_GF_2D_Lzero(lattice_type, a, e1, e2, kx, ky, alpha)
    else
        return ana_FT_GF_2D_Lnonzero(lattice_type, a, e1, e2, L, kx, ky)
    end
end

    
function ana_FT_GF_2D_Lzero(lattice_type, a, e1, e2, kx=0.0, ky=0.0, alpha=1e-7)
    # Calculate the self-energy/collective energies of the momentum modes of a planar lattice
    # The real contribution from the origin in the Fourier sum (ri = rj) is excluded (inf. Lamb shift)
    qx = kx/wa  
    qy = ky/wa
    Gk = 0.0 + 0.0im
    summands = Array{ComplexF64}(undef, 0)
    i = 0
    unit_cell_area = get_unit_cell_area(lattice_type, a)
    lattice_vectors = get_rec_lattice_vectors(lattice_type, a)
    
    while i < 4 || abs(sum(summands)) > 1e-4
        ms = get_ms(lattice_type, lattice_vectors, i)
        
        # Construct kappa_perp^2, kappa_z, and kappa vector
        kappa_perp2 = (qx .- ms[:, 1]).^2 + (qy .- ms[:, 2]).^2
        kappa_z = sqrt.(complex(1 .- kappa_perp2))
        kappa = stack([qx .- ms[:, 1], qy .- ms[:, 2]], dims=1)
        
        # Calculate summands (note that cross terms xz and yz are fixed to be zero)
        summands = (e1' * e2 .- transpose(e1[1:2]'*kappa) .* kappa'*e2[1:2] .- conj(e1[3])*kappa_z.^2*e2[3])./kappa_z
        
        # Multiply by appropriate exponential
        summands .*= exp.(-alpha*wa^2*kappa_perp2)
        
        # Perform sum
        Gk += sum(summands)
        i += 1
    end
    
    # Multiply with front constant to get minus the self-energy
    Gk *= 3im*π/(wa^2*unit_cell_area)
    
    # Subtract R[G(0, wa)]
    mat = Diagonal([1/2*(1 - 1/(2*alpha*wa^2)) , 1/2*(1 - 1/(2*alpha*wa^2)) , 1 + 1/(wa*sqrt(π*alpha))])
    Gk -= 3/(4*wa)*sqrt(π/alpha)*exp(-alpha*wa^2) * (e1' * mat * e2)  
     
    return Gk
end

    
function ana_FT_GF_2D_Lnonzero(lattice_type, a, e1, e2, L, kx=0.0, ky=0.0)
    # Calculate the self-energy/collective energies of the momentum modes of a planar lattice
    # The real contribution from the origin in the Fourier sum (ri = rj) is excluded (inf. Lamb shift)
    # The calculation assumes the sign of L to be positive
    qx = kx/wa  
    qy = ky/wa
    kpL = wa*L
    Gk = 0.0 + 0.0im
    summands = Array{ComplexF64}(undef, 0)
    i = 0
    unit_cell_area = get_unit_cell_area(lattice_type, a)
    lattice_vectors = get_rec_lattice_vectors(lattice_type, a)
    
    while i < 4 || abs(sum(summands)) > 1e-4
        ms = get_ms(lattice_type, lattice_vectors, i)
        
        # Construct kappa_perp^2, kappa_z, and kappa vector
        kappa_perp2 = (qx .- ms[:, 1]).^2 + (qy .- ms[:, 2]).^2
        kappa_z = sqrt.(complex(1 .- kappa_perp2))
        kappa = stack([qx .- ms[:, 1], qy .- ms[:, 2], kappa_z], dims=1)
        
        # Calculate summands
        summands = (e1' * e2 .- transpose(e1'*kappa) .* kappa'*e2)./kappa_z
        
        # Multiply by appropriate exponential
        summands .*= exp.(1im*kappa_z*kpL)
        
        # Perform sum
        Gk += sum(summands)
        i += 1
    end
    
    # Multiply with front constant to get minus the self-energy
    Gk *= 3im*π/(wa^2*unit_cell_area)
     
    return Gk
end


function ana_FT_GF_matrix(lattice_type, a, kx=0.0, ky=0.0, L=0.0, alpha=1e-7)
    if lattice_type ∈ ["square", "triangular"]
        return ana_FT_GF_2D_matrix(lattice_type, a, kx, ky, L, alpha)
    elseif lattice_type ∈ ["linear"]
        return ana_FT_GF_1D_matrix(a, kx)
    else
        throw(ArgumentError("lattice_type = $lattice_type has not been implemented in ana_FT_GF"))
    end
end


function ana_FT_GF_2D_matrix(lattice_type, a, kx=0.0, ky=0.0, L=0.0, alpha=1e-7)
    # Calculate the self-energy/collective energies of the momentum modes of a planar lattice
    # The real contribution from the origin in the Fourier sum (ri = rj) is excluded (inf. Lamb shift)
    # The calculation assumes the sign of L to be positive
    qx = kx/wa  
    qy = ky/wa
    kpL = wa*L
    Gk = fill(0.0 + 0.0im, 3, 3)
    summands = Array{ComplexF64}(undef, 0)
    i = 0
    unit_cell_area = get_unit_cell_area(lattice_type, a)
    lattice_vectors = get_rec_lattice_vectors(lattice_type, a)
    
    while i < 4 || any(abs.(sum(summands)) .> 1e-4)
        ms = get_ms(lattice_type, lattice_vectors, i)
        
        # Construct kappa_perp^2, kappa_z, and kappa vector
        kappa_perp2 = (qx .- ms[:, 1]).^2 + (qy .- ms[:, 2]).^2
        kappa_z = sqrt.(complex(1 .- kappa_perp2))
        kappa = [[qx - ms[i, 1], qy - ms[i, 2], kappa_z[i]] for i in eachindex(kappa_z)]
        
        # Calculate summands
        summands = (Ref(I) .- kappa .* adjoint.(kappa))./kappa_z
        
        # Multiply by appropriate exponential
        if L == 0
            summands .*= exp.(-alpha*wa^2*kappa_perp2)
        else
            summands .*= exp.(1j*kappa_z*kpL)
        end
        
        # Perform sum
        Gk += sum(summands)
        i += 1
    end
    
    # Multiply with front constant to get minus the self-energy
    Gk *= 3im*π/(wa^2*unit_cell_area)
    
    if L == 0
        # Subtract R[G(0, wa)]
        mat = Diagonal([1/2*(1 - 1/(2*alpha*wa^2)) , 1/2*(1 - 1/(2*alpha*wa^2)) , 1 + 1/(wa*sqrt(π*alpha))])
        Gk -= 3/(4*wa)*sqrt(π/alpha)*exp(-alpha*wa^2) * mat
        
        # Set xy and yz entries to zero (the above approach/formulas don't find the correct result for these specific entries)
        Gk[[3, 6, 7, 8]] .= 0
    end

    return Gk
end


function get_unit_cell_area(lattice_type, a)
    if lattice_type == "square"
        return a^2
    elseif lattice_type == "triangular"
        return sqrt(3)/2*a^2
    end
end


function get_rec_lattice_vectors(lattice_type, a)
    if lattice_type == "square"
        return [1, 0]/a, [0, 1]/a
    elseif lattice_type == "triangular"
        return [2/sqrt(3), 0]/a, [1/sqrt(3), 1]/a, [-1/sqrt(3), 1]/a
    end
end


function get_ms(lattice_type, lattice_vectors, i)
    if i == 0
        return [0;;0]
    end
    
    if lattice_type == "square"
        m_1D = -i:i - 1
        return stack(vcat(Ref( i*lattice_vectors[1]) .+ m_1D.*Ref(lattice_vectors[2]),
                          Ref( i*lattice_vectors[2]) .- m_1D.*Ref(lattice_vectors[1]),
                          Ref(-i*lattice_vectors[1]) .- m_1D.*Ref(lattice_vectors[2]),
                          Ref(-i*lattice_vectors[2]) .+ m_1D.*Ref(lattice_vectors[1])), dims=1)
        
    elseif lattice_type == "triangular"
        m_1D = 0:i - 1
        return stack(vcat(Ref( i*lattice_vectors[1]) .+ m_1D.*Ref(lattice_vectors[3]),
                          Ref( i*lattice_vectors[2]) .- m_1D.*Ref(lattice_vectors[1]),
                          Ref( i*lattice_vectors[3]) .- m_1D.*Ref(lattice_vectors[2]),
                          Ref(-i*lattice_vectors[1]) .- m_1D.*Ref(lattice_vectors[3]),
                          Ref(-i*lattice_vectors[2]) .+ m_1D.*Ref(lattice_vectors[1]),
                          Ref(-i*lattice_vectors[3]) .+ m_1D.*Ref(lattice_vectors[2])), dims=1)
    end
end


function ana_FT_GF_1D(a, e1, e2, kx=0.0)
    return e1' * ana_FT_GF_1D_matrix(a, kx) * e2
end


function ana_FT_GF_1D_matrix(a, kx=0.0)
    # This assumes a FT along x, i.e. a chain lattice along the x axis    
    # Get the parallel and perpendicular components of the energy shifts and decay rates
    tildeDelta_para, tildeDelta_perp = ana_tildeDeltas_1D(a, kx)
    tildeGamma_para, tildeGamma_perp = ana_tildeGammas_1D(a, kx)
    
    # Put together the result (note the sign on the real part, due to the defintion of the energy shift)
    return Diagonal([-tildeDelta_para + 1im*tildeGamma_para, -tildeDelta_perp + 1im*tildeGamma_perp, -tildeDelta_perp + 1im*tildeGamma_perp])
end


function ana_tildeDeltas_1D(a, kx=0.0)
    # # Calculate the self-energy/collective energies of the momentum modes of a 1D (linear) lattice
    tildeDelta_para = real( 
                            -          ( polylog(3, exp(1im*(wa + kx)*a)) + polylog(3, exp(1im*(wa - kx)*a)) )
                            + 1im*wa*a*( polylog(2, exp(1im*(wa + kx)*a)) + polylog(2, exp(1im*(wa - kx)*a)) )
                          )
    
    tildeDelta_perp = real( 
                                         polylog(3, exp(1im*(wa + kx)*a)) + polylog(3, exp(1im*(wa - kx)*a)) 
                            - 1im*wa*a*( polylog(2, exp(1im*(wa + kx)*a)) + polylog(2, exp(1im*(wa - kx)*a)) )
                            - wa^2*a^2*( polylog(1, exp(1im*(wa + kx)*a)) + polylog(1, exp(1im*(wa - kx)*a)) )
                          )
    
    
    # Multiply with front constant to get minus the self-energy
    tildeDelta_para *= 3/(wa^3*a^3)
    tildeDelta_perp *= 3/(2*wa^3*a^3)
    
    return tildeDelta_para, tildeDelta_perp
end


function ana_tildeGammas_1D(a, kx=0.0)
    # Calculate the self-energy/collective energies of the momentum modes of a 1D (linear) lattice
    tildeGamma_para = 0.0
    tildeGamma_perp = 0.0
    q0 = 2*π/a
    m = 0
    
    while abs(kx + q0*m) ≤ wa        
        # Perform sum
        tildeGamma_para += 1 - (kx + q0*m)^2/wa^2
        tildeGamma_perp += 1 + (kx + q0*m)^2/wa^2
        
        m += 1
    end
    
    # Multiply with front constant to get minus the self-energy
    tildeGamma_para *= 3*π/(2*wa*a)
    tildeGamma_perp *= 3*π/(4*wa*a)
    
    return tildeGamma_para, tildeGamma_perp
end


function realspace_GF(ri, rj, e1, e2)
    # return e1' * realspace_GF_matrix(ri, rj) * e2
    
    # [optimized calculation]
    GF_matrix = realspace_GF_matrix(ri, rj)
    result = 0.0im
    for b in 1:3, a in 1:3
        result += conj(e1[a]) * GF_matrix[a, b] * e2[b]
    end
    return result
end


function realspace_GF_matrix(ri, rj)
    # # Calculate relative vector and its norm
    # r_vec = ri - rj
    # r = norm(r_vec)
    
    # # If ri == rj return zero (this removes the corresponding terms in sums of the EoMs)
    # if isapprox(r, 0)
    #     return zeros(ComplexF64, 3,3)
    # end
    
    # # Prepare rr^dagger matrix
    # rr_hat = r_vec*r_vec'/r^2
    
    # # Calculate GF tensor
    # GF = exp(1im*wa*r)/(4π*r) * ( (1 + (1im*wa*r - 1)/(wa*r)^2)*I - (1 + 3*(1im*wa*r - 1)/(wa*r)^2)*rr_hat )
    
    # # Return the GF
    # return 6*π/wa*GF
    
    
    # Calculate norm of relative vector [optimized calculation]
    r = sqrt( (ri[1] - rj[1])^2 + (ri[2] - rj[2])^2 + (ri[3] - rj[3])^2 )
    
    # If ri == rj return zero (this removes the corresponding terms in sums of the EoMs)
    if isapprox(r, 0)
        return zeros(ComplexF64, 3,3)
    end

    # Calculate and return GF tensor [optimized calculation]
    GF = zeros(ComplexF64, 3, 3)
    for b in 1:3, a in 1:3
        GF[a, b] = 6*π/wa*exp(1im*wa*r)/(4π*r) * ( (1 + (1im*wa*r - 1)/(wa*r)^2)*(a == b) - (1 + 3*(1im*wa*r - 1)/(wa*r)^2)*(ri[a] - rj[a])*(ri[b] - rj[b])/r^2)
    end
    return GF
end
