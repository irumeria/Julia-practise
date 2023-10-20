using Plots

function get_energy_harmonic_oscillator(path, sum_energy=true)
    path_len = length(path)
    energies = [(path[i+1] - path[i])^2 + path[i+1]^2 for i = 1:path_len-2]
    sum_energy == false ? energies : sum(energies) 
end

function test()
    # hyper params
    path_length = 100
    x_len = path_length
    half_x_len = floor(Int, x_len/2) 
    epochs = 5000
    metropolis_factor = 1

    # initial
    path = zeros(path_length)
    energy_curves = zeros(epochs)
    phi0_distribution_another = zeros(x_len)
    phi0_distribution = zeros(x_len)

    former_energy = get_energy_harmonic_oscillator(path)
    for ep = 1:epochs
        # MC
        element = rand(1:path_length-1)
        change = 2.0 * (rand() - 0.5)
        path[element] += change

        energy = get_energy_harmonic_oscillator(path)
        if energy > former_energy ||
            exp(former_energy - energy) < rand() * metropolis_factor
            path[element] = path[element] - change
        end

        x = floor(Int, path[element]*16 + half_x_len)
        x = x > x_len ? x_len : x
        x = x < 1 ? 1 : x

        phi0_distribution_another[x] += 1
        former_energy = energy
        energy_curves[ep] = energy
    end

    energies = get_energy_harmonic_oscillator(path, false)
    phi0_distribution = [
        prod(exp.(-1 * energies[1:i])) 
        for i = 1:path_length-2 
    ]

    # plot result
    plot(1:x_len, phi0_distribution_another, xlabel = "x", ylabel = "probability") 
    savefig("phi0_1.png")
    plot(1:path_length-2, phi0_distribution, xlabel = "|x|", ylabel = "probability") 
    savefig("phi0_2.png")
    plot(1:epochs, energy_curves, xlabel = "epoch", ylabel = "E_all ∝ -S")    
    savefig("energy_curves.png")
end

test()

