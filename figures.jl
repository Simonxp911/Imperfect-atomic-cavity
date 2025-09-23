    

function fig_array(array)
    # Start figure 
    fig = Figure(size=(900, 600))
    
    # Make title and axis
    Label(fig[1, 1], L"The atomic array and fiber$$", tellwidth=false)
    xmax = maximum([site[1] for site in array])
    xmin = minimum([site[1] for site in array])
    ymax = maximum([site[2] for site in array])
    ymin = minimum([site[2] for site in array])
    zmax = maximum([site[3] for site in array])
    zmin = minimum([site[3] for site in array])
    Axis3(fig[2, 1], limits=(zmin, zmax, ymin, ymax, xmin, xmax), 
                     yreversed = true,
                     xlabel=L"$ z/λ $", 
                     ylabel=L"$ y/λ $", 
                     zlabel=L"$ x/λ $", 
                     aspect=(zmax - zmin, ymax - ymin, xmax - xmin)./maximum((zmax - zmin, ymax - ymin, xmax - xmin)))
    
    # Plot the atoms
    radius = 0.1
    θs = range(0, π, 20)
    φs = range(0, 2π, 20)
    xSph = radius.*[cos(φ)*sin(θ) for θ in θs, φ in φs]
    ySph = radius.*[sin(φ)*sin(θ) for θ in θs, φ in φs]
    zSph = radius.*[cos(θ) for θ in θs, φ in φs]
    for site in array
        surface!(zSph .+ site[3], ySph .+ site[2], xSph .+ site[1], colormap=:grays)
    end
    
    # Finish figure
    display(GLMakie.Screen(), fig)
end


function fig_Delta_scan(Delta_range, scan, SP)
    # Title and y-label
    titl = "lattice_type, N_sheets, radius, cc, L = $(SP.lattice_type), $(SP.N_sheets), $(SP.radius), $(SP.cut_corners), $(round(SP.L, sigdigits=4)) \n" *
           "ff, pos_unc_ratio, N_inst = $(SP.ff), $(SP.pos_unc_ratio), $(SP.N_inst) \n" *
           "drive, w0_ratio, e1, detec_mode = $(SP.drive_type), $(SP.w0_ratio), $(SP.e1_label), $(SP.detec_mode)"
    
    # Start figure
    fig = Figure(size=(900, 600))
    Label(fig[1, 1], titl, tellwidth=false)
    ax1 = Axis(fig[2, 1], limits=(extrema(Delta_range)..., 0, 1), 
               xlabel=L"$ Δ/γ $", 
               ylabel=L"Transmission coefficient, $ T=|t|^2 $")
    
    # Plot the single site cumulants
    for i in axes(scan, 1)
        plot!(ax1, Delta_range, scan[i, :], label=false)
    end
        
    # Finish figure
    display(GLMakie.Screen(), fig)
end


function fig_Delta_scan_stats(Delta_range, means, stds, T_inf_k0, T_inf_k, SP)
    # Prepare title and y-label
    titl = "lattice_type, N_sheets, radius, cc, L = $(SP.lattice_type), $(SP.N_sheets), $(SP.radius), $(SP.cut_corners), $(round(SP.L, sigdigits=4)) \n" *
           "ff, pos_unc_ratio, N_inst = $(SP.ff), $(SP.pos_unc_ratio), $(SP.N_inst) \n" *
           "drive, w0_ratio, e1, detec_mode = $(SP.drive_type), $(SP.w0_ratio), $(SP.e1_label), $(SP.detec_mode)"
    
    fig = Figure(size=(900, 600))
    Label(fig[1, 1], titl, tellwidth=false)
    ax1 = Axis(fig[2, 1], limits=(extrema(Delta_range)..., 0, 1), 
               xlabel=L"$ Δ/γ $", 
               ylabel=L"Transmission coefficient, $ T=|t|^2 $")
    
    # Plot finite system scan statistics 
    lines!(ax1, Delta_range, means, color=:blue, label="Fin. array")
    band!( ax1, Delta_range, means + stds, means - stds , color=(:blue, 0.35))
    axislegend(ax1, position=:rb)
    
    # Plot the infinite system, normal-incidence, plane-wave scan (with a separate y-axis)
    if false != T_inf_k0 && false != T_inf_k
        ax2 = Axis(f[2, 1], limits=(extrema(Delta_range)..., 0, 1), 
                   yticklabelcolor=:red, yaxisposition=:right)
        hidespines!(ax2)
        hidexdecorations!(ax2)
    
        lines!(ax2, Delta_range, T_inf_k0, linestyle=:dash, color=:red, linewidth=0.5, label="Inf. array, k=0")
        lines!(ax2, Delta_range, T_inf_k, linestyle=:dash, color=:black, linewidth=0.5, label="Inf. array, $(SP.drive_type)")
        axislegend(ax2, position=:lb)
    end
    
    display(GLMakie.Screen(), fig)
end


function fig_Delta_scan_stats_comparison(Delta_range, means, stds, SP)
    # Prepare colors and make title 
    colors = distinguishable_colors(SP.ff_specs[3]*SP.pos_unc_ratio_specs[3], [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
    titl = "lattice_type, N_sheets, radius, cc, L = $(SP.lattice_type), $(SP.N_sheets), $(SP.radius), $(SP.cut_corners), $(round(SP.L, sigdigits=4)) \n" *
           "N_inst = $(SP.N_inst) \n" *
           "drive, w0_ratio, e1, detec_mode = $(SP.drive_type), $(SP.w0_ratio), $(SP.e1_label), $(SP.detec_mode)"
           
    # Start figure
    fig = Figure(size=(900, 600))
    Label(fig[1, 1], titl, tellwidth=false)
    ax1 = Axis(fig[2, 1], limits=(extrema(Delta_range)..., 0, 1), 
               xlabel=L"$ Δ/γ $", 
               ylabel=L"Transmission coefficient, $ T=|t|^2 $")
    
    # Plot the finite system scan statistics
    for (j, pos_unc_ratio) in enumerate(SP.pos_unc_ratio_range)
        for (i, ff) in enumerate(SP.ff_range)
            lines!(ax1, Delta_range, means[i, j], color=colors[i + (j-1)*SP.ff_specs[3]], label=L"$ ff = %$(round(ff, sigdigits=2)) $, pos_unc $ = %$(round(pos_unc_ratio, sigdigits=2)) $")
            band!( ax1, Delta_range, means[i, j] + stds[i, j], means[i, j] - stds[i, j] , color=colors[i + (j-1)*SP.ff_specs[3]], alpha=0.35)
        end
    end
    
    # Finish figure
    axislegend(ax1, position=:lb)
    display(GLMakie.Screen(), fig)
end


function fig_Delta_scan_stats_comparison_fixed_pos_unc(Delta_range, means, stds, SP)
    # Prepare colors 
    colors = distinguishable_colors(SP.ff_specs[3], [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
    
    # Make a figure for each value of pos_unc_ratio
    for (j, pos_unc_ratio) in enumerate(SP.pos_unc_ratio_range)
        fig = Figure(size=(900, 600))
        titl = "lattice_type, N_sheets, radius, cc, L = $(SP.lattice_type), $(SP.N_sheets), $(SP.radius), $(SP.cut_corners), $(round(SP.L, sigdigits=4)) \n" *
               "pos_unc_ratio, N_inst = $(pos_unc_ratio), $(SP.N_inst) \n" *
               "drive, w0_ratio, e1, detec_mode = $(SP.drive_type), $(SP.w0_ratio), $(SP.e1_label), $(SP.detec_mode)"
        Label(fig[1, 1], titl, tellwidth=false)
        ax1 = Axis(fig[2, 1], limits=(extrema(Delta_range)..., 0, 1), 
                xlabel=L"$ Δ/γ $", 
                ylabel=L"Transmission coefficient, $ T=|t|^2 $")
        
        # Plot the finite system scan statistics, one for each value of ff
        for (i, ff) in enumerate(SP.ff_range)
            lines!(ax1, Delta_range, means[i, j], color=colors[i], label=L"$ ff = %$(round(ff, sigdigits=2)) $")
            band!( ax1, Delta_range, means[i, j] + stds[i, j], means[i, j] - stds[i, j] , color=colors[i], alpha=0.35)
        end
        
        # Finish figure
        axislegend(ax1, position=:lb)
        display(GLMakie.Screen(), fig)
    end
end


function fig_Delta_scan_stats_comparison_fixed_ff(Delta_range, means, stds, SP)
    # Prepare colors
    colors = distinguishable_colors(SP.pos_unc_ratio_specs[3], [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
    
    # Make a figure for each value of pos_unc_ratio
    for (i, ff) in enumerate(SP.ff_range)
        fig = Figure(size=(900, 600))
        titl = "lattice_type, N_sheets, radius, cc, L = $(SP.lattice_type), $(SP.N_sheets), $(SP.radius), $(SP.cut_corners), $(round(SP.L, sigdigits=4)) \n" *
               "ff, N_inst = $(ff), $(SP.N_inst) \n" *
               "drive, w0_ratio, e1, detec_mode = $(SP.drive_type), $(SP.w0_ratio), $(SP.e1_label), $(SP.detec_mode)"
        Label(fig[1, 1], titl, tellwidth=false)
        ax1 = Axis(fig[2, 1], limits=(extrema(Delta_range)..., 0, 1), 
                xlabel=L"$ Δ/γ $", 
                ylabel=L"Transmission coefficient, $ T=|t|^2 $")
        
        # Plot the finite system scan statistics, one for each value of ff
        for (j, pos_unc_ratio) in enumerate(SP.pos_unc_ratio_range)
            lines!(ax1, Delta_range, means[i, j], color=colors[j], label=L"pos unc ratio $ = %$(round(pos_unc_ratio, sigdigits=2)) $")
            band!( ax1, Delta_range, means[i, j] + stds[i, j], means[i, j] - stds[i, j] , color=colors[j], alpha=0.35)
        end
        
        # Finish figure
        axislegend(ax1, position=:lb)
        display(GLMakie.Screen(), fig)
    end
end


function fig_Efield_intensity(x_range, y_range, z_range, intensity_xz, intensity_xy, array)
    # Plot the intensities
    for (intensity, array_site_indices, ranges, labels) in zip((intensity_xz, intensity_xy), 
                                                               ((1, 3), (1, 2)),
                                                               ((x_range, z_range), (x_range, y_range)),
                                                               ((L"$ x/λ $", L"$ z/λ $"), (L"$ x/λ $", L"$ y/λ $")))
        
        # Start figure 
        width  = maximum(ranges[1]) - minimum(ranges[1])
        height = maximum(ranges[2]) - minimum(ranges[2])
        fig = Figure(size=(width/height*600, 600))
        
        # Make title and axis
        # titl = L"$ a = %$(round(SP.a, sigdigits=3)) $, $ L = %$(round(SP.L, sigdigits=3)) $, $ \Delta = %$(round(Δ, sigdigits=3)) $"
        Label(fig[1, 1], "title", tellwidth=false)
        Axis(fig[2, 1], limits=(extrema(ranges[1]), extrema(ranges[2])), 
                        xlabel=labels[1], 
                        ylabel=labels[2], 
                        aspect=DataAspect())
        
        # Plot the E-field intensity
        # sat = maximum(intensity[length(x_range)÷2, Int(round(length(z_range)*3.5/8)):Int(round(length(z_range)*4.5/8))])
        sat = 1
        intensity[intensity .> sat] .= sat
        hm = heatmap!(ranges[1], ranges[2], intensity, colormap =:viridis)
        Colorbar(fig[2, 2], hm)
        
        # Plot a representation of the atomic array
        scatter!([site[array_site_indices[1]] for site in array], [site[array_site_indices[2]] for site in array], color=:black, marker=:circle, markersize=10)
        
        # Finish figure
        display(GLMakie.Screen(), fig)
    end
end


function fig_Efield_intensity_3D(x_range, y_range, z_range, intensities, array, detec_z, detec_radius)
    # Cut off intensities at some saturation value
    for intensity in intensities
        # sat = maximum(intensity[length(x_range)÷2, Int(round(length(z_range)*3.5/8)):Int(round(length(z_range)*4.5/8))])
        # sat = 0.1
        sat = maximum(intensities[3])
        intensity[intensity .> sat] .= sat
    end
    
    # Start figure 
    fig = Figure(size=(900, 600))
    
    # Make title and axis
    # titl = "a = $(round(SP.a, sigdigits=3)), L = $(round(SP.L, sigdigits=3)), Delta = $(round(Δ, sigdigits=3))"
    Label(fig[1, 1], "title", tellwidth=false)
    zWidth  = maximum(z_range) - minimum(z_range)
    xHeight = maximum(x_range) - minimum(x_range)
    yDepth  = maximum(y_range) - minimum(y_range)
    Axis3(fig[2, 1], limits=(extrema(z_range), extrema(x_range), extrema(y_range)), 
                     yreversed = true,
                     xlabel=L"$ z/λ $", 
                     ylabel=L"$ x/λ $", 
                     zlabel=L"$ y/λ $", 
                     aspect=(zWidth, xHeight, yDepth)./maximum((zWidth, xHeight, yDepth)))
    
        
    # Plot the intensity surfaces
    n = length(x_range)
    # yz-plane (really xz)
    xx, yy, zz = x_range.*zeros(n)',  y_range.*ones(n)', z_range.*ones(n)'
    surface!(zz', xx, yy, color=intensities[1], colormap=:viridis)
    
    # xz-plane (really xy)
    xx, yy, zz = x_range.*ones(n)',  y_range.*zeros(n)', z_range.*ones(n)'
    surface!(zz', xx, yy, color=intensities[2], colormap=:viridis)
    
    # detection plane xy-plane (really yz)
    xx, yy, zz = x_range.*ones(n)',  y_range.*ones(n)', detec_z*ones(n, n)
    detection_plane = xx'.^2 + yy.^2 .<= detec_radius^2
    flt = ones(size(detection_plane)); flt[.!detection_plane] .= NaN        
    surface!(zz.*flt, xx'.*flt, yy.*flt, color=intensities[3], colormap=:viridis)
    
    # end xy-plane (really yz)
    xx, yy, zz = x_range.*ones(n)',  y_range.*ones(n)', z_range[end]*ones(n, n)
    surface!(zz, xx', yy, color=intensities[4], colormap=:viridis)
    
    # Plot the atoms
    radius = 0.3
    θs = range(0, π, 20)
    φs = range(0, 2π, 20)
    xSph = radius.*[cos(φ)*sin(θ) for θ in θs, φ in φs]
    ySph = radius.*[sin(φ)*sin(θ) for θ in θs, φ in φs]
    zSph = radius.*[cos(θ) for θ in θs, φ in φs]
    for site in array
        surface!(zSph .+ site[3], ySph .+ site[2], xSph .+ site[1], colormap=:grays)
    end
    
    # Finish figure
    display(GLMakie.Screen(), fig)
end 