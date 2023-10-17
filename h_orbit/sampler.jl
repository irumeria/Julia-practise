using ProgressMeter

include("orbit.jl")
include("constants.jl")

struct Sampler
    orbital
    grid_size
    sample_amount
    cell_volume
    norm_factor
end

struct Posi
    z_posi
    x_posi
end

function newSampler(n, l, m, grid_size)
    r_bound = H_ATOM_RADIUS
    sample_amount = 1e8
    orbital = newOrbital(n, l, m)
    cell_volume = (r_bound / (grid_size / 2.0 - 1.0))^3
    norm_factor = r_bound / (grid_size / 2.0 - 1.0)

    Sampler(
        orbital,
        grid_size,
        sample_amount,
        cell_volume,
        norm_factor,
    )
end

function sample_xz_plane(sampler::Sampler, y)
    grid = zeros((sampler.grid_size, sampler.grid_size))
    posiGrid = fill(Posi(0, 0), (sampler.grid_size, sampler.grid_size))

    @showprogress for i = 1:sampler.grid_size
        z = (sampler.grid_size / 2 - i) *
            sampler.norm_factor
        for j = 1:sampler.grid_size
            x = (sampler.grid_size / 2 - j) *
                sampler.norm_factor
            posiGrid[i, j] = Posi(z, x)

            coord = [x, y, z]
            (p, phase) = probability_with_phase(
                sampler.orbital,
                coord,
                sampler.cell_volume)
            # This calculates the probability at this point after sample_amount of sampling
            p = 1.0 - (1.0 - p)^(sampler.sample_amount)
            target = rand()
            if p > target
                grid[i, j] = phase
            end
        end
    end

    (grid, posiGrid)

end

