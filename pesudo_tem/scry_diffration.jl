using Plots
using LinearAlgebra

# Each row in the matrix stands for a cell vector
# some example lattices with a=1
const global FCC = 0.5 .* [1 1 0; 0 1 1; 1 0 1]
const global BCC = 0.5 .* [1 1 -1; -1 1 1; 1 -1 1]
const global CUBE = diagm([1, 1, 1])

const global FCC_BASIC_PRIMITIVE_COORDINATE = [0 0 0; 0 1/2 1/2; 1/2 0 1/2; 1/2 1/2 0]
const global BCC_BASIC_PRIMITIVE_COORDINATE = [0 0 0; 1/2 1/2 1/2]

function get_sqr_complex(complex_number)
    real_part = real(complex_number)
    img_part = real(im * complex_number)
    real_part^2 + img_part^2
end

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
    in_hkl=true
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

function diffaction_strength(coords::Matrix, v::Vector)
    # Calculate the 'S_G' * distantce_to_center
    # assume that the strength decay in normal distribution
    sigma = 0.4 # this parameter simulates the size of Ewald's sphere
    distance_factor = 1 / sigma * exp(-sum(v .^ 2) / 2 * sigma^2)
    get_sqr_complex(
        sum([
            exp(
                im * 2 * pi *
                sum(coords[pindex, :] .* v)
            )
            for pindex = 1:size(coords)[1]
        ])
    ) * distance_factor
end

function get_diffraction_recip_point(
        recip_points, 
        primitive_coords
    )

    diffraction_recip_point = []
    strengths = []
    hkls = []
    for (hkl, rp) in recip_points
        hkl_array = [parse(Int, s) for s in split(hkl, ",")]
        strength = diffaction_strength(primitive_coords, hkl_array) # TODO: if I need to abs.(hkl) ?
        if abs(strength) < 1e-6
            continue
        end
        diffraction_recip_point = [diffraction_recip_point; reshape(rp, 1, length(rp))]
        strengths = [strengths; strength]
        hkls = [hkls; hkl]
    end
    diffraction_recip_point, strengths, hkls
end


alpha = 0.287
cell = alpha .* BCC
primitive_coords = BCC_BASIC_PRIMITIVE_COORDINATE

# println(latpoints)
@show recip_cell = lattice2recip(cell)
@show diffaction_strength(primitive_coords, [2, 0, 2]) # for test

recip_points = expandLattice(recip_cell, 2)

diffraction_point, strengths, hkls = get_diffraction_recip_point(
    recip_points,
    primitive_coords
)

# selected the hkl with l=0
reserved_indexs = vec(mapslices(col -> abs(col[3]) < 1e-6, diffraction_point, dims=2))
diffraction_point = diffraction_point[reserved_indexs, :]
strengths = strengths[reserved_indexs, :]
hkls = hkls[reserved_indexs, :]

p = scatter(diffraction_point[:, 1], diffraction_point[:, 2],
    markersize=strengths, markercolor=:blue, markerstrokewidth=0, legend=false,
    series_annotations = text.(hkls, :bottom))

title!(p, "electron diffraction partarn, BCC alpha=" * string(alpha))

savefig(p, "BCC with alpha=" * string(alpha) * ".png")

