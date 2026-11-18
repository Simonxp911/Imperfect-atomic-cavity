    

function fig_array(array)
    # Start figure 
    fig = Figure(size=(900, 600))
    
    # Make title and axis
    Label(fig[1, 1], L"The atomic array and fiber$$", tellwidth=false)
    zoom = 1.1
    shift = 0.5
    xmax = maximum([site[1] for site in array]) * zoom + shift
    xmin = minimum([site[1] for site in array]) * zoom - shift
    ymax = maximum([site[2] for site in array]) * zoom + shift
    ymin = minimum([site[2] for site in array]) * zoom - shift
    zmax = maximum([site[3] for site in array]) * zoom + shift
    zmin = minimum([site[3] for site in array]) * zoom - shift
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
    title = "lattice_type, N_sheets, radius, cc, L = $(SP.AP.lattice_type), $(SP.AP.N_sheets), $(SP.AP.radius), $(SP.AP.cut_corners), $(round(SP.AP.L, sigdigits=4)) \n" *
           "ff, pos_unc_ratio, N_inst = $(SP.AP.ff), $(SP.pos_unc_ratio), $(SP.AP.N_inst) \n" *
           "drive, w0_ratio, dipoleMoment, detec_type = $(SP.DrP.drive_type), $(SP.w0_ratio), $(SP.EP.dipoleMoment_label), $(SP.DeP.detec_type)"
    
    # Start figure
    fig = Figure(size=(900, 600))
    Label(fig[1, 1], title, tellwidth=false)
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
    title = "lattice_type, N_sheets, radius, cc, L = $(SP.AP.lattice_type), $(SP.AP.N_sheets), $(SP.AP.radius), $(SP.AP.cut_corners), $(round(SP.AP.L, sigdigits=4)) \n" *
           "ff, pos_unc_ratio, N_inst = $(SP.AP.ff), $(SP.pos_unc_ratio), $(SP.AP.N_inst) \n" *
           "drive, w0_ratio, dipoleMoment, detec_type = $(SP.DrP.drive_type), $(SP.w0_ratio), $(SP.EP.dipoleMoment_label), $(SP.DeP.detec_type)"
    
    fig = Figure(size=(900, 600))
    Label(fig[1, 1], title, tellwidth=false)
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
        lines!(ax2, Delta_range, T_inf_k, linestyle=:dash, color=:black, linewidth=0.5, label="Inf. array, $(SP.DrP.drive_type)")
        axislegend(ax2, position=:lb)
    end
    
    display(GLMakie.Screen(), fig)
end


function fig_Delta_TRscan_stats(Delta_range, T_means, T_stds, R_means, R_stds, SP)
    # Prepare title and y-label
    title = "lattice_type, N_sheets, radius, cc, L = $(SP.AP.lattice_type), $(SP.AP.N_sheets), $(SP.AP.radius), $(SP.AP.cut_corners), $(round(SP.AP.L, sigdigits=4)) \n" *
           "ff, pos_unc_ratio, N_inst = $(SP.AP.ff), $(SP.pos_unc_ratio), $(SP.AP.N_inst) \n" *
           "drive, w0_ratio, dipoleMoment, detec_type = $(SP.DrP.drive_type), $(SP.w0_ratio), $(SP.EP.dipoleMoment_label), $(SP.DeP.detec_type)"
    
    # Start figure
    fig = Figure(size=(900, 600))
    Label(fig[1, 1:2], title, tellwidth=false)
    ax1 = Axis(fig[2, 1], limits=(extrema(Delta_range)..., 0, 1), 
               xlabel=L"$ Δ/γ $", 
               ylabel=L"Transmission coefficient, $ T=|t|^2 $")
    ax2 = Axis(fig[2, 2], limits=(extrema(Delta_range)..., 0, 1), 
               xlabel=L"$ Δ/γ $", 
               ylabel=L"Reflection coefficient, $ R=|r|^2 $")
    
    # Plot 
    lines!(ax1, Delta_range, T_means, color=:blue)
    band!( ax1, Delta_range, T_means + T_stds, T_means - T_stds , color=(:blue, 0.35))
    lines!(ax2, Delta_range, R_means, color=:red)
    band!( ax2, Delta_range, R_means + R_stds, R_means - R_stds , color=(:red, 0.35))
    
    # Finish figure
    display(GLMakie.Screen(), fig)
end


function fig_scanStatsComparison(Delta_range, T_means, T_stds, R_means, R_stds, labels, title)
    # Prepare colors
    colors = distinguishable_colors(length(labels), [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
    
    #    Start figure
    fig = Figure(size=(1200, 600))
    Label(fig[1, 1:2], title, tellwidth=false)
    ax1 = Axis(fig[2, 1], limits=(extrema(Delta_range)..., 0, 1), 
               xlabel=L"$ Δ/γ $", 
               ylabel=L"Transmission coefficient, $ T=|t|^2 $")
    ax2 = Axis(fig[2, 2], limits=(extrema(Delta_range)..., 0, 1), 
               xlabel=L"$ Δ/γ $", 
               ylabel=L"Reflection coefficient, $ R=|r|^2 $")
    
    for (i, label) in enumerate(labels)
        lines!(ax1, Delta_range, T_means[i], color=colors[i], label=label)
        band!( ax1, Delta_range, T_means[i] + T_stds[i], T_means[i] - T_stds[i] , color=colors[i], alpha=0.35)
        lines!(ax2, Delta_range, R_means[i], color=colors[i], label=label)
        band!( ax2, Delta_range, R_means[i] + R_stds[i], R_means[i] - R_stds[i] , color=colors[i], alpha=0.35)
    end
    
    # Finish figure
    Legend(fig[3, 1:2], ax1, orientation=:horizontal, nbanks=2) 
    display(GLMakie.Screen(), fig)
    # save("C:\\Users\\Simon\\Downloads\\explore_a370_varyROI.png", fig)
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
        # title = L"$ a = %$(round(SP.AP.a, sigdigits=3)) $, $ L = %$(round(SP.AP.L, sigdigits=3)) $, $ \Delta = %$(round(Δ, sigdigits=3)) $"
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
    # title = "a = $(round(SP.AP.a, sigdigits=3)), L = $(round(SP.AP.L, sigdigits=3)), Delta = $(round(Δ, sigdigits=3))"
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
    
    # detection trans plane xy-plane (really yz)
    xx, yy, zz = x_range.*ones(n)',  y_range.*ones(n)', detec_z*ones(n, n)
    detection_plane = xx'.^2 + yy.^2 .<= detec_radius^2
    flt = ones(size(detection_plane)); flt[.!detection_plane] .= NaN        
    surface!(zz.*flt, xx'.*flt, yy.*flt, color=intensities[3], colormap=:viridis)
    
    # detection refl plane xy-plane (really yz)
    xx, yy, zz = x_range.*ones(n)',  y_range.*ones(n)', -detec_z*ones(n, n)
    detection_plane = xx'.^2 + yy.^2 .<= detec_radius^2
    flt = ones(size(detection_plane)); flt[.!detection_plane] .= NaN        
    surface!(zz.*flt, xx'.*flt, yy.*flt, color=intensities[4], colormap=:viridis)
    
    # end xy-plane (really yz)
    xx, yy, zz = x_range.*ones(n)',  y_range.*ones(n)', z_range[end]*ones(n, n)
    surface!(zz, xx', yy, color=intensities[5], colormap=:viridis)
    
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


function fig_Efield_intensity_3D_crossSections(x_range, y_range, z_range, intensities, array, detec_z, detec_radius)
    # Expect the same input as fig_Efield_intensity_3D, which means some input is not used
    
    # Cut off intensities at some saturation value
    for intensity in intensities
        # sat = maximum(intensity[length(x_range)÷2, Int(round(length(z_range)*3.5/8)):Int(round(length(z_range)*4.5/8))])
        # sat = 0.1
        sat = maximum(intensities[3])
        intensity[intensity .> sat] .= sat
    end
    
    # Start figure 
    fig = Figure(size=(1200, 900))
    
    # Make title and axis
    # title = "a = $(round(SP.AP.a, sigdigits=3)), L = $(round(SP.AP.L, sigdigits=3)), Delta = $(round(Δ, sigdigits=3))"
    Label(fig[1, 1:2], "title", tellwidth=false)
    zWidth  = maximum(z_range) - minimum(z_range)
    xHeight = maximum(x_range) - minimum(x_range)
    yDepth  = maximum(y_range) - minimum(y_range)
    ax1 = Axis3(fig[2, 1], limits=(extrema(z_range), extrema(x_range), extrema(y_range)), 
                    yreversed = true,
                    xlabel=L"$ z/λ $", 
                    ylabel=L"$ x/λ $", 
                    zlabel=L"$ y/λ $", 
                    aspect=(zWidth, xHeight, yDepth)./maximum((zWidth, xHeight, yDepth)),
                    azimuth=0.85*2π)
    ax2 = Axis(fig[2, 2], limits=(extrema(x_range), extrema(y_range)), 
                    xlabel=L"$ x/λ $", 
                    ylabel=L"$ y/λ $",
                    aspect=xHeight/yDepth)
    ax3 = Axis(fig[3, 1], limits=(extrema(z_range), extrema(y_range)), 
                    xlabel=L"$ z/λ $", 
                    ylabel=L"$ y/λ $",
                    aspect=zWidth/yDepth)
    ax4 = Axis(fig[3, 2], limits=(extrema(z_range), extrema(x_range)), 
                    xlabel=L"$ z/λ $", 
                    ylabel=L"$ x/λ $",
                    aspect=zWidth/xHeight)
    
        
    # Plot the intensity surfaces
    n = length(x_range)
    # yz-plane (really xz)
    xx, yy, zz = x_range.*zeros(n)',  y_range.*ones(n)', z_range.*ones(n)'
    surface!(ax1, zz', xx, yy, color=intensities[1], colormap=:viridis, shading=false)
    heatmap!(ax3, z_range, y_range, intensities[1]', colormap=:viridis)
    
    # xz-plane (really xy)
    xx, yy, zz = x_range.*ones(n)',  y_range.*zeros(n)', z_range.*ones(n)'
    surface!(ax1, zz', xx, yy, color=intensities[2], colormap=:viridis, shading=true)
    heatmap!(ax4, z_range, x_range, intensities[2]', colormap=:viridis)
    
    # end xy-plane (really yz)
    xx, yy, zz = x_range.*ones(n)',  y_range.*ones(n)', z_range[end]*ones(n, n)
    surface!(ax1, zz, xx', yy, color=intensities[5], colormap=:viridis, shading=false)
    heatmap!(ax2, x_range, y_range, intensities[5], colormap=:viridis)
    
    # Plot the atoms
    radius = 0.3
    θs = range(0, π, 20)
    φs = range(0, 2π, 20)
    xSph = radius.*[cos(φ)*sin(θ) for θ in θs, φ in φs]
    ySph = radius.*[sin(φ)*sin(θ) for θ in θs, φ in φs]
    zSph = radius.*[cos(θ) for θ in θs, φ in φs]
    for site in array
        surface!(ax1, zSph .+ site[3], ySph .+ site[2], xSph .+ site[1], colormap=:grays)
        scatter!(ax3, [site[3]], [site[2]], color=RGBA(0, 0, 0, 0), markersize=10, strokewidth=1, strokecolor=:black)
        scatter!(ax4, [site[3]], [site[1]], color=RGBA(0, 0, 0, 0), markersize=10, strokewidth=1, strokecolor=:black)
        
    end
    
    # Colorbar(fig[3, 2:3], limits=colorrange, colormap=:viridis, label=L"Intensity, any unit$$")
    
    # Finish figure
    display(GLMakie.Screen(), fig)
    save("C:\\Users\\Simon\\Downloads\\Efield_intensity_3D.png", fig)
end 


function fig_excitation_distribution(sites_x, sites_y, sheets_z, N_sheets, pops, title)
    
    # Start figure
    fig = Figure(size=(1200, 600))
    Label(fig[1, 1:N_sheets+1], title, tellwidth=false)
    axes = [
            Axis(fig[2, i][1, 1], #limits=(extrema(Delta_range)..., 0, 1), 
            xlabel=L"$ x/λ $")
           for i in 1:N_sheets]
    axes[1].ylabel = L"$ y/λ $"
    linkyaxes!(axes...)
    
    colorrange = [minimum([pop for pop in pops if !isnan(pop)]), maximum([pop for pop in pops if !isnan(pop)])]
    for (i, ax, pop, sheet_z) in zip(1:N_sheets, axes, eachrow(pops), sheets_z)
        # heatmap!(ax, sites_x, sites_y, pop)
        # tricontourf!(ax, sites_x, sites_y, pop, colormap=:viridis)
        scatter!(ax, sites_x, sites_y, color=pop, markersize=25, strokewidth=1, strokecolor=:black, colormap=:viridis, colorrange=colorrange)
        ax.title = L"$ z = %$(round(sheet_z, sigdigits=4)) $"
    end
    # Colorbar(fig[2, N_sheets+1], sc)
    Colorbar(fig[2, N_sheets+1], limits=colorrange, colormap=:viridis, label=L"Population$$")
    
    # Finish figure
    display(GLMakie.Screen(), fig)
    # save("C:\\Users\\Simon\\Downloads\\.png", fig)
end


    
    # save("C:\\Users\\Simon\\Downloads\\FILENAME.png", fig)