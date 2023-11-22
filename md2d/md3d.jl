using Distributions
using Random
using ProgressMeter
using CairoMakie
using LinearAlgebra

include("crystal.jl")

mutable struct SystemState
	posi_array::Matrix
	mass::AbstractFloat
	velocity_array::Matrix
	cell_length::AbstractFloat
end

function plot_distribution(sys::SystemState, name)
	fig = Figure()
	ax = Axis3(fig[1, 1],
		title = "distribution",
		xlabel = "x",
		ylabel = "y",
		zlabel = "z")
	scatter!(ax,
		sys.posi_array[:, 1],
		sys.posi_array[:, 2],
		sys.posi_array[:, 3])
	save(name, fig)
end

function sign(a, b)
	ret = abs(a)
	if b < 0.0
		ret = -abs(a)
	end
	ret
end

function _LJFroce(
	posi_1::Vector,
	posi_2::Vector,
	ep::AbstractFloat,
	sigma::AbstractFloat,
	cell_length::AbstractFloat;
	rcut = Inf,
)
	# default: return the force of the atom in posi_1
	posi_delta = posi_1 - posi_2
	for i = eachindex(posi_delta)
		if abs(posi_delta[i]) > 0.5 * cell_length
			posi_delta[i] = posi_delta[i] - sign(cell_length, posi_delta[i])
		end
	end
	r_delta = sum(sqrt.(posi_delta .^ 2))
	if r_delta > rcut
		return zeros(3)
	end
	LJ_partial_r = -4 * ep * (
					   12 * sigma^12 / r_delta^14 -
					   6 * sigma^6 / r_delta^8
				   )
	LJ_partial_delta_distance = LJ_partial_r .* posi_delta

	LJ_partial_delta_distance
end

function calculateLJForces(sys::SystemState, ep, sigma, rcut)
	force_array = zeros(size(sys.posi_array))
	atom_number = size(sys.posi_array)[1]
	for i ∈ 1:atom_number
		for j ∈ (i+1):atom_number
			force = _LJFroce(
				sys.posi_array[i, :],
				sys.posi_array[j, :],
				ep,
				sigma,
				sys.cell_length;
				rcut = rcut,
			)
			# @show size(force_array[i,:]) size(force)
			force_array[i, :] = force
			force_array[j, :] = -1 .* force
		end
	end
	force_array
end

function updateVelocity(
	sys::SystemState,
	force_array,
	dt,
)
	acceleration_array = force_array / sys.mass
	new_velocity_array = sys.velocity_array + acceleration_array * dt
	new_velocity_array
end

function updatePosi(sys::SystemState, dt)
	(sys.posi_array + sys.velocity_array * dt) .% sys.cell_length
end

function do_FCC_simulation(steps::Int)

	d = Normal()
	dt = 0.000005
	half_dt = dt * 0.5

	lattice_size = 2
	cell_length = lattice_size + 0.2

	fcc = expandCoords(FCC_COORDINATE, lattice_size)

	fcc = map(e -> e |> float, fcc)
	sys = SystemState(
		fcc,
		1.0,
		rand(d, size(fcc)) .* 1e-3,
		cell_length,
	)

	ep = 1e-5
	sigma = 1 / sqrt(2)
	rcut = cell_length / 2.1

	plot_distribution(sys, "figures/init.png")

	@showprogress for step in 1:steps
		force_array = calculateLJForces(sys, ep, sigma, rcut)
		sys.velocity_array = updateVelocity(sys, force_array, half_dt)
		sys.posi_array = updatePosi(sys, dt)
		if step % 1000 == 0
			plot_distribution(sys, "figures/distribution_" * string(step) * ".png")
		end
	end

	plot_distribution(sys, "figures/final.png")
end

do_FCC_simulation(10000)
