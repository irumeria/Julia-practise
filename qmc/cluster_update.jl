using Distributions
using BenchmarkTools
using ProgressMeter

mutable struct SystemState
    MC_steps::Int32         # Number of monte carlo steps
    eq_data_points::Int32   # Number of equilibrium data points in MC steps
    L::Int32                # Grid size in 1 dimension
    T::Float32              # Initial temperature
    T_steps::Int32        # Number of temperature steps
    dT::Float32             # Temperature increment
    J::Float32              # Coupling J (Keep at 1)
    kb::Float32             # Boltzman constant (Keep at 1)
    spin_init               # Initial spin grid (up, down or random)
    algorithm               # SW (Swendesen Wang) or SF (Spin Flip)
    bs_trials::Int32        # Number of boostrap trials, between 500-2000 suffieces
    spin_site_total_number  # Total number of spin sites
    MCS                     # Montecarlo step size
    time_steps              # Montecarlo time to regular time steps
    function SystemState(
        MC_steps=1000,
        L=40,
        T=3.5,
        T_steps=8,
        dT=-0.025,
        J=1,
        kb=1,
        spin_init="up",
        algorithm="SW",
        bs_trials=1000,
    )
        spin_site_total_number = L^2
        MCS = spin_site_total_number
        time_steps = MCS * MC_steps
        eq_data_points = MC_steps

        @assert T + dT * T_steps > 0
        @assert algorithm == "SW"
        @assert L < 45

        new(MC_steps, eq_data_points, L, T, T_steps, dT, J, kb, spin_init, algorithm, bs_trials, spin_site_total_number, MCS, time_steps)
    end
end

function nonzero_in_matrix(x)
    inds = Tuple.(findall(!iszero, x))
    a = first.(inds)
    b = last.(inds)
    cat(a, b; dims=2)
end

function spin_site_energy(sys, x, y, grid_spins)
    spin_site_energy = 0
    # spin_neigbour_x = (spin_site_x + np.array([1, 0, -1, 0]))%(self.L)
    # spin_neigbour_y = (spin_site_y + np.array([0, -1, 0, 1]))%(self.L)
    
    n_y_up = y + 1 > sys.L ? y + 1 - sys.L : y + 1
    n_y_down = y - 1 < 1 ? sys.L + (y - 1) : y - 1
    n_x_left = x - 1 < 1 ? sys.L + (x - 1) : x - 1
    n_x_right = x + 1 > sys.L ? x + 1 - sys.L : x + 1

    neighbours_coords = [
        x n_y_down;
        x n_y_up;
        n_x_left y;
        n_x_right y
    ]

    for i  = 1:4
        nx = neighbours_coords[i, 1]
        ny = neighbours_coords[i, 2]
        spin_value_center = grid_spins[x, y]
        spin_value_neighbour = grid_spins[nx, ny]

        spin_site_energy += -sys.J*spin_value_center*spin_value_neighbour 
    end
    spin_site_energy 
end

function system_energy(sys, grid_coordinates, grid_spins)
    
    sys_energy = 0
    for i = 1:sys.L
        for j = 1:sys.L
            spin_site_x = grid_coordinates[1][i]
            spin_site_y = grid_coordinates[2][j]
            sys_energy += spin_site_energy(sys, spin_site_x, spin_site_y, grid_spins)
        end
    end
    sys_energy = sys_energy/2 # To counter double counting of the links
    
    sys_energy
end

function swendsen_wang(
    sys::SystemState,
    grid_coordinates::Vector{<:Vector{<:Int}},
    grid_spins::Matrix
)
    islands = []
    cluster_flips = []
    bonds = bond_eval(sys, grid_spins)
    not_visited = fill(true, (sys.L, sys.L))

    for i = 1:sys.L
        for j = 1:sys.L
            flip_cluster = 2 .* rand(0:1) - 1 # -1 or 1
            # @show i j
            spin_site_x = grid_coordinates[1][i] 
            spin_site_y = grid_coordinates[2][j] 
            # println(spin_site_x, spin_site_y)
            # @show "-1" spin_site_x spin_site_y
            cluster, grid_spins, not_visited = back_track(
                sys,
                spin_site_x,
                spin_site_y,
                bonds,
                grid_spins,
                flip_cluster,
                not_visited
            )

            if length(cluster) > 0
                append!(islands, [cluster]) # TODO
                append!(cluster_flips, flip_cluster)
            end
        end
    end
    islands, grid_spins, cluster_flips
end

function bond_eval(sys::SystemState, grid_spins::Matrix)
    #=
    Goes over all the spins in the system and checks the bonds
    bonds : 3D array (2, L, L). contains al the bonds present in the system. 
    The first 2D array gives the horizontal bonds
    The second 2D array gives the vertical bonds
    Element (1,i,j) gives the relation
        between spin_site (i,j) and (i,j+1)
    =#
    bonds = zeros(Float32, (2, sys.L, sys.L))
    chance_value = minimum([1, exp(-2 * sys.J / (sys.kb * sys.T))])

    delta_spin_hor = abs.(
        grid_spins +
        circshift(grid_spins, (0, -1))
    ) / 2
    non_zero_index_hor = nonzero_in_matrix(delta_spin_hor)

    delta_spin_ver = abs.(
        grid_spins +
        circshift(grid_spins, (-1, 0))
    ) / 2
    non_zero_index_ver = nonzero_in_matrix(delta_spin_ver)

    bino = Binomial(1, chance_value)
    for i = 1:size(non_zero_index_hor)[1]
        if rand(bino) == 1
            bonds[1, non_zero_index_hor[i, 1], non_zero_index_hor[i, 2]] = 0
        else
            bonds[1, non_zero_index_hor[i, 1], non_zero_index_hor[i, 2]] = Inf
        end
    end

    for i = 1:size(non_zero_index_ver)[1]
        if rand(bino) == 1
            bonds[2, non_zero_index_ver[i, 1], non_zero_index_ver[i, 2]] = 0
        else
            bonds[2, non_zero_index_ver[i, 1], non_zero_index_ver[i, 2]] = Inf
        end
    end

    # @show sum([rand(bino) for _ = 1:3200])

    # @show sum([if abs(e) > 1e-6 1 else 0 end for e in bonds])

    bonds
end

function back_track(
    sys::SystemState,
    x::Int,
    y::Int,
    bonds::Array{<:AbstractFloat,3},
    grid_spins::Matrix,
    flip_cluster::Int,
    not_visited::Matrix{Bool};
    cluster=[]
)
    #=
    Checks the neighbours of the spin, if they are 
    equal this functions jumps over to that spin and 
    repeats itself. The spins that are already visited 
    are skipped. Everytime an equal bond is found, this
    spind is added to the cluster.
    =#
    # @show "0" x y
    if not_visited[x, y]
        # @show x y
        not_visited[x, y] = false
        append!(cluster, [[x, y]])
        grid_spins[x, y] = grid_spins[x, y] * flip_cluster
        n_y_up = y + 1 > sys.L ? y + 1 - sys.L : y + 1
        n_y_down = y - 1 < 1 ? sys.L + (y - 1) : y - 1
        n_x_left = x - 1 < 1 ? sys.L + (x - 1) : x - 1
        n_x_right = x + 1 > sys.L ? x + 1 - sys.L : x + 1
        # @show n_y_up n_y_down n_x_left n_x_right
        @assert n_y_up <= sys.L && y <= sys.L
        if bonds[1, x, y] == Inf
            n_x = x
            n_y = n_y_up
            # @show "1" n_x n_y
            cluster, grid_spins, not_visited = back_track(
                sys,
                n_x,
                n_y,
                bonds,
                grid_spins,
                flip_cluster,
                not_visited;
                cluster
            )
        end
        if bonds[1, x, n_y_down] == Inf
            n_x = x
            n_y = n_y_down
            # @show "2" n_x n_y
            cluster, grid_spins, not_visited = back_track(
                sys,
                n_x,
                n_y,
                bonds,
                grid_spins,
                flip_cluster,
                not_visited;
                cluster
            )
        end
        if bonds[2, x, y] == Inf
            n_x = n_x_right
            n_y = y
            # @show "3" n_x n_y
            cluster, grid_spins, not_visited = back_track(
                sys,
                n_x,
                n_y,
                bonds,
                grid_spins,
                flip_cluster,
                not_visited;
                cluster
            )
        end
        if bonds[2, n_x_left, y] == Inf
            n_x = n_x_left
            n_y = y
            # @show "4" n_x n_y
            cluster, grid_spins, not_visited = back_track(
                sys,
                n_x,
                n_y,
                bonds,
                grid_spins,
                flip_cluster,
                not_visited;
                cluster
            )
        end
    end
    cluster, grid_spins, not_visited
end

function assign_spin(sys::SystemState)
    if sys.spin_init == "random"
        grid_spins = rand(sys.L, sys.L)
        grid_spins[map(e -> e >= 0.5, grid_spins)] = 1
        grid_spins[map(e -> e >= 0.5, grid_spins)] = -1
    elseif sys.spin_init == "up"
        grid_spins = ones(Int32, (sys.L, sys.L))
    else
        sys.spin_init == "down"
        grid_spins = -1 .* ones(Int32, (sys.L, sys.L))
    end
    grid_spins
end

function grid_init(sys::SystemState)

end

function MC_simulation(sys::SystemState)
    grid_spins = assign_spin(sys)
    grid_coordinates = [
        range(1, sys.L) |> collect,
        range(1, sys.L) |> collect
    ]

    T_total = zeros((sys.T_steps, 1))
    h_total = zeros((sys.T_steps, 1))
    energy = zeros((sys.T_steps, 1))
    magnetisation = zeros((sys.T_steps, 2))
    chi = zeros((sys.T_steps, 2))
    c_v = zeros((sys.T_steps, 2))

    energy_i = zeros((sys.MC_steps, 1))
    magnetisation_i = zeros((sys.MC_steps, 1))
    chi_i = zeros((sys.MC_steps, 1))

    @showprogress for (j, temp) in enumerate(range(1, sys.T_steps))
        for (i, t) in enumerate(range(1, sys.MC_steps))
            islands, grid_spins, cluster_flips =
                swendsen_wang(
                    sys,
                    grid_coordinates,
                    grid_spins
                )
            energy_i[i] = system_energy(sys, grid_coordinates, grid_spins)
        end
    end
    @show energy_i
end

sys = SystemState()
sys.J

MC_simulation(sys)

