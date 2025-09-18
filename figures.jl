

function fig_array(array)
    # Extract x and y coordinates of sites
    x = [site[1] for site in array]
    y = [site[2] for site in array]
    z = [site[3] for site in array]
    
    # Start figure
    fig = plot(reuse=false)
    
    # Plot the single site cumulants
    scatter!(x, y, z)
    # scatter!(x, y, z, texts=eachindex(x))
    
    # Finish figure
    plot!(legend=false, aspect_ratio=:equal, ticks=:native)
    display(fig)
end


function fig_Delta_scan(Delta_range, scan, SP)
    
    # Start figure
    fig = plot(reuse=false, size=(900, 600), ticks=:native)
    
    # Plot the single site cumulants
    for i in axes(scan, 1)
        plot!(Delta_range, scan[i, :], label=false)
    end
    
    # Title and y-label
    titl = "lattice_type, N_sheets, radius, cc, L = $(SP.lattice_type), $(SP.N_sheets), $(SP.radius), $(SP.cut_corners), $(round(SP.L, sigdigits=4)) \n" *
           "ff, pos_unc_ratio, N_inst = $(SP.ff), $(SP.pos_unc_ratio), $(SP.N_inst) \n" *
           "drive, w0_ratio, e1, detec_mode = $(SP.drive_type), $(SP.w0_ratio), $(SP.e1_label), $(SP.detec_mode)"
    ylabl = L"Transmission coefficient, $ T=|t|^2 $"
    
    # Finish figure
    xlims!(Delta_range[1], Delta_range[end])
    ylims!(0, 1)
    xlabel!(L"$ \Delta/\gamma $")
    ylabel!(ylabl)
    title!(titl)
    display(fig)
    
end


function fig_Delta_scan_stats(Delta_range, means, stds, T_inf_k0, T_inf_k, SP)
    # Prepare title and y-label
    titl = "lattice_type, N_sheets, radius, cc, L = $(SP.lattice_type), $(SP.N_sheets), $(SP.radius), $(SP.cut_corners), $(round(SP.L, sigdigits=4)) \n" *
           "ff, pos_unc_ratio, N_inst = $(SP.ff), $(SP.pos_unc_ratio), $(SP.N_inst) \n" *
           "drive, w0_ratio, e1, detec_mode = $(SP.drive_type), $(SP.w0_ratio), $(SP.e1_label), $(SP.detec_mode)"
    ylabl = L"Transmission coefficient, $ T=|t|^2 $"
    
    if false == T_inf_k0 == T_inf_k
        # Start figure
        fig = plot(reuse=false, size=(900, 600), ticks=:native)
        
        # Plot the finite system scan statistics
        plot!(Delta_range, means, ribbon=stds, fillalpha=0.35, c=1, label=false)
        
        # Add labels, etc. (do it before adding the second axis to avoid double labels...)
        xlims!(Delta_range[1], Delta_range[end])
        xlabel!(L"$ \Delta/\gamma $")
        ylabel!(ylabl)
        title!(titl)
    else
        # Start figure
        fig = plot(reuse=false, size=(900, 600), ticks=:native)
        
        # Plot the infinite system, normal-incidence, plane-wave scan
        plot!(Delta_range, T_inf_k0, linestyle=:dash, c=:red, linewidth=0.5, label="Inf. array, k=0")
        plot!(Delta_range, T_inf_k, linestyle=:dash, c=:black, linewidth=0.5, label="Inf. array, $(SP.drive_type)")
        plot!(legend=:bottomleft)
        
        # Add labels, etc. (do it before adding the second axis to avoid double labels...)
        xlims!(Delta_range[1], Delta_range[end])
        ylims!(0, 1)
        xlabel!(L"$ \Delta/\gamma $")
        ylabel!(ylabl)
        title!(titl)
        
        # Plot the finite system scan statistics (with a separate y-axis)
        plot!(twinx(), Delta_range, means, ribbon=stds, fillalpha=0.35, c=1, xticks=:none, label="Fin. array", legend=:bottomright)
        ylims!(0, 1)
    end
    display(fig) 
end


function fig_Delta_scan_stats_comparison(Delta_range, means, stds, SP)
    
    # Prepare colors for plotting
    colors = distinguishable_colors(SP.ff_specs[3]*SP.pos_unc_ratio_specs[3], [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
    # colors = palette(:rainbow, SP.ff_specs[3])
    
    # Start figure
    fig = plot(reuse=false, size=(900, 600), ticks=:native)
    
    # Plot the finite system scan statistics
    for (j, pos_unc_ratio) in enumerate(SP.pos_unc_ratio_range)
        for (i, ff) in enumerate(SP.ff_range)
            plot!(Delta_range, means[i, j], ribbon=stds[i, j], fillalpha=0.35, c=colors[i + (j-1)*SP.ff_specs[3]], label=L"$ ff = %$(round(ff, sigdigits=2)) $, pos_unc $ = %$(round(pos_unc_ratio, sigdigits=2)) $")
        end
    end
        
    # Make title 
    titl = "lattice_type, N_sheets, radius, cc, L = $(SP.lattice_type), $(SP.N_sheets), $(SP.radius), $(SP.cut_corners), $(round(SP.L, sigdigits=4)) \n" *
           "N_inst = $(SP.N_inst) \n" *
           "drive, w0_ratio, e1, detec_mode = $(SP.drive_type), $(SP.w0_ratio), $(SP.e1_label), $(SP.detec_mode)"
    
    # Finish figure
    xlims!(Delta_range[1], Delta_range[end])
    # ylims!(0, 1)
    xlabel!(L"$ \Delta/\gamma $")
    ylabel!(L"Transmission coefficient, $ T=|t|^2 $")
    title!(titl)
    display(fig)
    
end


function fig_Delta_scan_stats_comparison_fixed_pos_unc(Delta_range, means, stds, SP)
    
    # Prepare colors for plotting
    colors = distinguishable_colors(SP.ff_specs[3], [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
    # colors = palette(:rainbow, SP.ff_specs[3])
    
    # Make a figure for each value of pos_unc_ratio
    for (j, pos_unc_ratio) in enumerate(SP.pos_unc_ratio_range)
        fig = plot(reuse=false, size=(900, 600), ticks=:native)    
        
        # Plot the finite system scan statistics, one for each value of ff
        for (i, ff) in enumerate(SP.ff_range)
            plot!(Delta_range, means[i, j], ribbon=stds[i, j], fillalpha=0.35, c=colors[i], label=L"$ ff = %$(round(ff, sigdigits=2)) $")
        end
        
        # Make title 
        titl = "lattice_type, N_sheets, radius, cc, L = $(SP.lattice_type), $(SP.N_sheets), $(SP.radius), $(SP.cut_corners), $(round(SP.L, sigdigits=4)) \n" *
               "pos_unc_ratio, N_inst = $(pos_unc_ratio), $(SP.N_inst) \n" *
               "drive, w0_ratio, e1, detec_mode = $(SP.drive_type), $(SP.w0_ratio), $(SP.e1_label), $(SP.detec_mode)"
        
        # Finish figure
        xlims!(Delta_range[1], Delta_range[end])
        # ylims!(0, 1)
        xlabel!(L"$ \Delta/\gamma $")
        ylabel!(L"Transmission coefficient, $ T=|t|^2 $")
        title!(titl)
        display(fig)
    end
end


function fig_Delta_scan_stats_comparison_fixed_ff(Delta_range, means, stds, SP)
    
    # Prepare colors for plotting
    colors = distinguishable_colors(SP.pos_unc_ratio_specs[3], [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
    # colors = palette(:rainbow, SP.ff_specs[3])
    
    # Make a figure for each value of pos_unc_ratio
    for (i, ff) in enumerate(SP.ff_range)
        fig = plot(reuse=false, size=(900, 600), ticks=:native)    
        
        # Plot the finite system scan statistics, one for each value of ff
        for (j, pos_unc_ratio) in enumerate(SP.pos_unc_ratio_range)
            plot!(Delta_range, means[i, j], ribbon=stds[i, j], fillalpha=0.35, c=colors[j], label="pos_unc_ratio = " * L"$%$(round(pos_unc_ratio, sigdigits=2)) $")
        end
        
        # Make title 
        titl = "lattice_type, N_sheets, radius, cc, L = $(SP.lattice_type), $(SP.N_sheets), $(SP.radius), $(SP.cut_corners), $(round(SP.L, sigdigits=4)) \n" *
               "ff, N_inst = $(ff), $(SP.N_inst) \n" *
               "drive, w0_ratio, e1, detec_mode = $(SP.drive_type), $(SP.w0_ratio), $(SP.e1_label), $(SP.detec_mode)"
            
        # Finish figure
        xlims!(Delta_range[1], Delta_range[end])
        # ylims!(0, 1)
        xlabel!(L"$ \Delta/\gamma $")
        ylabel!(L"Transmission coefficient, $ T=|t|^2 $")
        title!(titl)
        display(fig)
    end
end


function fig_Efield_intensity(SP)
    # Get collective energies and choose the detuning of perfect transmission
    Gk = ana_FT_GF(SP.lattice_type, SP.a, SP.e1, SP.e1)
    tildeDelta = -real(Gk)
    tildeGamma =  imag(Gk)
    Δ = tildeDelta - tildeGamma*tan(wa*SP.L)
    
    # Find the steady state coherences
    σ_ss = calc_σ_ss(Δ, SP.array[1], SP.N, SP.e1, SP.drivemode[1])
    
    # Define x, y, and z ranges for the plot
    x_range = range(-3*SP.radius*SP.a, 3*SP.radius*SP.a, 101)
    y_range = deepcopy(x_range)
    z_range = range(-3*SP.L, 3*SP.L, 101)
    
    # Calculate the E-field intensity  (in the xz and the zy planes)
    intensity_xz = zeros(length(x_range), length(z_range))
    for (i, x) in enumerate(x_range), (j, z) in enumerate(z_range)
        r = [x, 0.0, z]
        
        # We calculate E-field multiplied by d and divided by incoming amplitude
        Ed = calc_total_Efield_fin(r, SP.array[1], σ_ss, SP.N, SP.drive_type, SP.w0, SP.e1, "forward")
        
        intensity_xz[i, j] = Ed'*Ed
    end
    
    intensity_xy = zeros(length(x_range), length(y_range))
    for (i, x) in enumerate(x_range), (j, y) in enumerate(y_range)
        r = [x, y, maximum(z_range)]
        
        # We calculate E-field multiplied by d and divided by incoming amplitude
        Ed = calc_total_Efield_fin(r, SP.array[1], σ_ss, SP.N, SP.drive_type, SP.w0, SP.e1, "forward")
        
        intensity_xy[i, j] = Ed'*Ed
    end
    
    
    # Plot the intensities
    for (intensity, array_sites, ranges) in zip((intensity_xz, intensity_xy), 
                                                ([[site[1], site[3]] for site in SP.array[1]], [[site[1], site[2]] for site in SP.array[1]]),
                                                ((x_range, z_range), (x_range, y_range)))
        # Cut off intensity at some saturation for plotting purposes (the intensity at the positions of the atoms is very high)
        # sat = maximum(intensity[length(x_range)÷2, Int(round(length(z_range)*3.5/8)):Int(round(length(z_range)*4.5/8))])
        sat = 5
        intensity[intensity .> sat] .= sat
    
        fig = heatmap(ranges[1], ranges[2], transpose(intensity), 
                    c = :viridis,
                    #   levels=50,
                    colorbar=true,
                    size=(1000, 700),
                    reuse=false)
                    
        # Plot the atomic sites
        for site in array_sites
            scatter!([site[1]], [site[2]], c=:black, markershape=:circle, label=false)
        end
        
        # Finish figure
        title!(L"$ a = %$(round(SP.a, sigdigits=3)) $, $ L = %$(round(SP.L, sigdigits=3)) $, $ \Delta = %$(round(Δ, sigdigits=3)) $")
        display(fig)
    end
    
end


function fig_Efield_intensity_3D(SP)
    # Get collective energies and choose the detuning of perfect transmission
    Gk = ana_FT_GF(SP.lattice_type, SP.a, SP.e1, SP.e1)
    tildeDelta = -real(Gk)
    tildeGamma =  imag(Gk)
    Δ = tildeDelta - tildeGamma*tan(wa*SP.L)
    
    # Find the steady state coherences
    σ_ss = calc_σ_ss(Δ, SP.array[1], SP.N, SP.e1, SP.drivemode[1])
    
    # Define x, y, and z ranges for the plot
    n = 101
    x_range = range(-3*SP.radius*SP.a, 3*SP.radius*SP.a, n)
    y_range = deepcopy(x_range)
    z_range = range(-5*SP.L, 10*SP.L, n)
    
    # Calculate the E-field intensity (in the xz, yz, and xy planes)
    intensities = [zeros(n, n) for i in 1:4]
    for i in 1:n, j in 1:n
        for (r, intensity) in zip(([x_range[i], 0.0, z_range[j]],
                                   [0.0, y_range[i], z_range[j]],
                                   [x_range[i], y_range[j], SP.detec_z],
                                   [x_range[i], y_range[j], z_range[end]]),
                                   intensities)
            # We calculate E-field multiplied by d and divided by incoming amplitude
            Ed = calc_total_Efield_fin(r, SP.array[1], σ_ss, SP.N, SP.drive_type, SP.w0, SP.e1, "forward")
            
            intensity[i, j] = Ed'*Ed
        end
    end
    
    # Cut off intensities at some saturation value
    for intensity in intensities
        # sat = maximum(intensity[length(x_range)÷2, Int(round(length(z_range)*3.5/8)):Int(round(length(z_range)*4.5/8))])
        # sat = 0.1
        sat = maximum(intensities[3])
        intensity[intensity .> sat] .= sat
    end
    
    # Temporarily use the plotlyjs backend to make an interactive 3D plot
    Plots.with(:plotlyjs) do
        # Start figure
        fig = plot(reuse=false, size=(1000, 1000))
        
        # Plot the intensity surfaces
        # yz-plane (really xz)
        xx, yy, zz = x_range.*zeros(n)',  y_range.*ones(n)', z_range.*ones(n)'
        surface!(zz', xx, yy, surfcolor=intensities[1], c=:viridis)
        
        # xz-plane (really xy)
        xx, yy, zz = x_range.*ones(n)',  y_range.*zeros(n)', z_range.*ones(n)'
        surface!(zz', xx, yy, surfcolor=intensities[2], c=:viridis)
        
        # detection plane xy-plane (really yz)
        xx, yy, zz = x_range.*ones(n)',  y_range.*ones(n)', SP.detec_z*ones(n, n)
        detection_plane = xx'.^2 + yy.^2 .<= SP.detec_radius^2
        flt = ones(size(detection_plane)); flt[.!detection_plane] .= NaN        
        surface!(zz.*flt, xx'.*flt, yy.*flt, surfcolor=intensities[3], c=:viridis)
        
        # end xy-plane (really yz)
        xx, yy, zz = x_range.*ones(n)',  y_range.*ones(n)', z_range[end]*ones(n, n)
        surface!(zz, xx', yy, surfcolor=intensities[4], c=:viridis)
        
        # Plot the atomic sites
        xs = [site[1] for site in SP.array[1]]
        ys = [site[2] for site in SP.array[1]]
        zs = [site[3] for site in SP.array[1]]
        scatter!(zs, xs, ys, c=:black, markershape=:circle, label=false)
        
        # Finish figure
        plot!(camera=(-30,30), aspect_ratio=:equal, colorbar=false)
        xlabel!("z")
        ylabel!("x")
        zlabel!("y")
        xlims!(extrema(z_range))
        ylims!(extrema(x_range))
        zlims!(extrema(y_range))
        title!("a = $(round(SP.a, sigdigits=3)), L = $(round(SP.L, sigdigits=3)), Delta = $(round(Δ, sigdigits=3))")
        display(fig)
    end
end