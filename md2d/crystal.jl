using LinearAlgebra

const global FCC_BASIC_PRIMITIVE_COORDINATE = [0 0 0; 0 1/2 1/2; 1/2 0 1/2; 1/2 1/2 0]

const global FCC_COORDINATE = [
    0 0 0; 
    1 0 0;
    0 1 0;
    0 0 1;
    1 1 0;
    1 0 1;
    0 1 1;
    1 1 1;
    0 1/2 1/2; 
    1/2 0 1/2; 
    1/2 1/2 0;
    1 1/2 1/2;
    1/2 1/2 1;
    1/2 1 1/2;
    1/2 1/2 1/2
]

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
    include_negative=false,
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

function expandCoords(
        posi_array::Matrix,
        factor::Int;
        include_negative=false
    )
    expanded_posi = copy(posi_array)
    start_factor = include_negative ? -factor : 0
    for efx = start_factor:factor
        for efy = start_factor:factor
            for efz = start_factor:factor
                new_posi = zeros(size(posi_array))
                new_posi[:,1] = posi_array[:, 1] .+ efx
                new_posi[:,2] = posi_array[:, 2] .+ efy
                new_posi[:,3] = posi_array[:, 3] .+ efz
                expanded_posi = [expanded_posi; new_posi]
            end
        end
    end
    unique(expanded_posi, dims = 1)

end

# @show size(expandedCoords(FCC_COORDINATE, 1))