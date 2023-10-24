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
        T_steps=81,
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


function swendsen_wang(
    sys::SystemState,
    grid_coordinates::Vector{<:Vector{<:AbstractFloat}},
    grid_spins::Matrix
)
    islands = []
    cluster_flips = []
    not_visited = ones((sys.L, sys.L))
    bonds = bond_eval(sys, grid_spins)

    for i = 1:sys.spin_site_total_number
        for j = 1:sys.spin_site_total_number 
            cluster = []
            flip_cluster = 2 .* rand(0:1) - 1 # -1 or 1
            spin_site_x = grid_coordinates[1][i]
            spin_site_y = grid_coordinates[2][j]
            cluster, grid_spins = back_track(sys, spin_site_x, spin_site_y, bonds, not_visited, cluster, grid_spins, flip_cluster)

            if length(cluster) > 0
                append!(islands, [cluster]) # TODO
                append!(cluster_flips, flip_cluster)
            end
        end
    end
    islands, grid_spins, cluster_flips
end

function bond_eval(sys::SystemState, grid_spins::Matrix)
    # TODO
end

function back_track(
    sys::SystemState,
    x,
    y,
    bonds,
    not_visited,
    cluster,
    grid_spins,
    flip_cluster
)

end

function assign_spin(sys::SystemState)
    if sys.spin_init == "random"
        grid_spins = rand(sys.L, sys.L)
        grid_spins[map(e -> e>=0.5 , grid_spins)] = 1
        grid_spins[map(e -> e>=0.5 , grid_spins)] = -1
    elseif sys.spin_init == "up"
        grid_spins = ones(Int32, (sys.L, sys.L))
    else sys.spin_init == "down"
        grid_spins = -1 .* ones(Int32, (sys.L, sys.L))
    end
    grid_spins
end

function grid_init(sys::SystemState)
    
end

function MC_simulation(sys::SystemState)
    grid_spins = assign_spin(sys)
    grid_coordinates = [
        range(1, sys.L) |> collect |> float,
        range(1, sys.L) |> collect |> float
        ]
    grid_coordinates
    T_total = zeros((sys.T_steps, 1))
    h_total = zeros((sys.T_steps, 1))
    energy = zeros((sys.T_steps, 1))
    magnetisation = zeros((sys.T_steps, 2))
    chi = zeros((sys.T_steps, 2))
    c_v = zeros((sys.T_steps, 2))       
    
    energy_i = zeros((sys.MC_steps, 1))
    magnetisation_i = zeros((sys.MC_steps, 1))
    chi_i = zeros((sys.MC_steps, 1))

    for (j, temp) in enumerate(range(1, sys.T_steps))
        for (i, t) in enumerate(range(1, sys.MC_steps))
            islands, grid_spins, cluster_flips = 
                swendsen_wang(
                    sys, 
                    grid_coordinates, 
                    grid_spins
                )
            # energy_i[i] =  
        end
    end
end

@show sys = SystemState()
sys.J

MC_simulation(sys)