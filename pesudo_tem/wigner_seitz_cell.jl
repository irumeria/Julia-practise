
function lattice2recip(lattice_vectors::Matrix)
    reciprocal_lattice = zeros(3, 3)
    detlattice = det(lattice_vectors)
    reciprocal_lattice[:, 1] = 2 * pi *
                               cross(lattice_vectors[:, 2], lattice_vectors[:, 3]) /
                               detlattice
    reciprocal_lattice[:, 2] = 2 * pi *
                               cross(lattice_vectors[:, 3], lattice_vectors[:, 1]) /
                               detlattice
    reciprocal_lattice[:, 3] = 2 * pi *
                               cross(lattice_vectors[:, 1], lattice_vectors[:, 2]) /
                               detlattice
    reciprocal_lattice
end

function expandLattice(
    cell::Matrix,
    factor::Int;
    include_negative=true,
    exclude_origin=false,
    in_hkl=false
)
    start_factor = include_negative ? -factor : 0
    expanded_lat = in_hkl ? Dict() : []
    for efx = start_factor:factor
        for efy = start_factor:factor
            for efz = start_factor:factor
                if exclude_origin && efx == 0 && efy == 0 && efz == 0
                    continue
                end
                new_lat = efx .* cell[1, :] + efy .* cell[2, :] + efz .* cell[3, :]
                if in_hkl
                    hkl = string(efx) * "," * string(efy) * "," * string(efz)
                    expanded_lat[hkl] = new_lat
                else
                    expanded_lat = [expanded_lat; reshape(new_lat, 1, length(new_lat))]
                end

            end
        end
    end
    expanded_lat
end

using LinearAlgebra

function to_voronoi(points::Matrix, i::Int; return_lines=false)
    num_points = size(points, 1)
    voronoi_lines = []
    center_point = points[i, :]

    for j in 1:num_points
        if j == i
            continue
        end
        midpoint = (center_point + points[j, :]) / 2
        direction = points[j, :] - center_point
        perpendicular = [-direction[2], direction[1]]
        line = [midpoint - 100 * perpendicular, midpoint + 100 * perpendicular]
        theta = atan(perpendicular[2], perpendicular[1])
        if return_lines
            push!(voronoi_lines, line)
        else
            push!(voronoi_lines, (midpoint, theta))
        end
        
    end

    voronoi_lines
end


lattice = [1 0 0; 1/2 1 0; 1 0 1]
recip_lattice = lattice2recip(lattice)

lattice = expandLattice(lattice,1)
@show lattice

recip_lattice = expandLattice(recip_lattice,1)
@show recip_lattice

# @show draw_voronoi(recip_lattice[:, 1:2])
wigner_borders = to_voronoi(lattice[:, 1:2], Int((size(lattice)[1]+1)/2); return_lines=true)
wigner_borders_recip = to_voronoi(recip_lattice[:, 1:2], Int((size(recip_lattice)[1]+1)/2); return_lines=true)

using Plots

function draw_lines_with_points(lines, points, name, limit)
    for line in lines
        x = [line[1][1], line[2][1]]
        y = [line[1][2], line[2][2]]
        xlims!(-limit, limit)
        ylims!(-limit, limit)
        plot!(x, y, color=:blue, legend = false)
    end
    scatter!(points[:, 1], points[:, 2], color=:red, legend = false)
    
    savefig(name)
end

draw_lines_with_points(wigner_borders, lattice, "wigner_borders.png", 5)
p = plot()
draw_lines_with_points(wigner_borders_recip, recip_lattice, "wigner_borders_recip.png", 10)
