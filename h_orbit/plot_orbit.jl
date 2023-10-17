using Plots

include("sampler.jl")

function getSamples(n, l, m; grid_size=1024)
    sampler = newSampler(n, l, m, grid_size)
    sample_xz_plane(sampler, 0)
end

function plot_orbit(n, l, m, grid_size=1024)
    (grid, _) = getSamples(n, l, m; grid_size)
    posi_phase_points_x = []
    posi_phase_points_z = []
    nega_phase_points_x = []
    nega_phase_points_z = []
    for (i, cell) in enumerate(grid)
        row_i = i % grid_size 
        col_i = trunc(Int, i / grid_size)
        if cell > 1e-6
            append!(posi_phase_points_x, [col_i])
            append!(posi_phase_points_z, [row_i])
        elseif cell < -1e-6
            append!(nega_phase_points_x, [col_i])
            append!(nega_phase_points_z, [row_i])
        end
    end
    p = scatter(nega_phase_points_x, nega_phase_points_z, 
            markersize=0.5,markercolor=:blue, markerstrokewidth=0)

    p = scatter(p, posi_phase_points_x, posi_phase_points_z, 
        markersize=0.5,markercolor=:red, markerstrokewidth=0)

    title!(p, "lables: n="*string(n)*" l="*string(l)*" m="*string(m))

    println(length(posi_phase_points_x))
    println(length(nega_phase_points_x))
    savefig(p, "n="*string(n)*" l="*string(l)*" m="*string(m)*".png")   
end

plot_orbit(10,7,0)

