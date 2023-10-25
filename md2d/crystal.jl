
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