import DataFrames: DataFrame
using CairoMakie
using ProgressMeter

global atom_number = 25
global T_init = 2.0

global density = 1.0
global x = zeros(Float32, atom_number)
global y = zeros(Float32, atom_number)
global vx = zeros(Float32, atom_number)
global vy = zeros(Float32, atom_number)
global fx = zeros(Float32, (2, atom_number))
global fy = zeros(Float32, (2, atom_number))

global L = trunc(Int, 1.0 * atom_number^0.5)
global epochs = 10000
# atoms = Vector(epochs)

global T_states = zeros(epochs)
global P_states = zeros(epochs)

function twelveran()  # Average 12 rands for Gaussian
	s = 0.0
	for _ in 1:12
		s += rand()
	end

	s / 12.0 - 0.5
end

function initialposvel()          # Initialize
	i = 0
	for ix ∈ 1:L                             # x->   1  2  3  4  5
		for iy in 1:L                        # y=1   0  5  10 15 20 
			i = i + 1                        # y=2   1  6  11 16 21
			x[i] = ix                        # y=3   2  7  12 17 22
			y[i] = iy                        # y=4   3  8  13 18 23
			vx[i] = twelveran()              # y=5   4  9  14 19 24
			vy[i] = twelveran()              # numbering of 25 atoms
			vx[i] = vx[i] * sqrt(T_init)
			vy[i] = vy[i] * sqrt(T_init)
		end
	end
end


function sign(a, b)
	ret = abs(a)
	if (b < 0.0)
		ret = -abs(a)
	end
	ret
end

function Forces(t, w; returnPE = false)   # Forces
	# invr2 = 0. 
	r2cut = 9.0                              # Switch: PEorW = 1 for PE   
	PE = 0.0
	for i ∈ 1:atom_number
		fx[t, i] = fy[t, i] = 0.0
	end
	for i ∈ 1:atom_number
		for j ∈ (i+1):atom_number
			dx = x[i] - x[j]
			dy = y[i] - y[j]
			if abs(dx) > 0.50 * L
				dx = dx - sign(L, dx)  # Interact with closer image
			end
			if abs(dy) > 0.50 * L
				dy = dy - sign(L, dy)
			end
			r2 = dx * dx + dy * dy
			if r2 < r2cut
				if r2 == 0.0  # To avoid 0 denominator
					r2 = 1e-4
				end
				invr2 = 1.0 / r2
				wij = 48.0 * (invr2^3 - 0.5) * invr2^3
				fijx = wij * invr2 * dx
				fijy = wij * invr2 * dy
				fx[t, i] = fx[t, i] + fijx
				fy[t, i] = fy[t, i] + fijy
				fx[t, j] = fx[t, j] - fijx
				fy[t, j] = fy[t, j] - fijy
				PE = PE + 4.0 * (invr2^3) * ((invr2^3) - 1.0)
				w = w + wij

			end
		end
	end

	ret = w
	if returnPE
		ret = PE
	end

	ret
end

function plot_states(T_states, P_states, total_time, name)
	time_step = total_time / length(T_states)
	fig = Figure()
	ax = Axis(fig[1, 1],
		title = "MD result",
		xlabel = "time",
		ylabel = "y/ev")
	line_data = DataFrame(
		x = collect(time_step:time_step:total_time),
		T = T_states,
		P = P_states,
	)
	lines!(ax, line_data[:, :x], line_data[:, :P], label = "Potential", color = :blue)
	lines!(ax, line_data[:, :x], line_data[:, :T], label = "Travail", color = :red)
	save(name, fig)
end

function plot_distribution(name)
	global x
	global y
	fig = Figure()
	ax = Axis(fig[1, 1],
		title = "distribution",
		xlabel = "x",
		ylabel = "y")
	scatter!(ax, x[:], y[:])
	save(name, fig)
end

function timevolution()
	time_1 = 1
	PE = 0.0
	h = 0.031          # step
	h_half = h / 2.0
	# initial KE & PE via Forces
	KE = 0.0
	w = 0.0
	initialposvel()
	for i ∈ 1:atom_number
		KE = KE + (vx[i] * vx[i] + vy[i] * vy[i]) / 2.0
	end
	PE = Forces(time_1, w, returnPE = true)
	time = 1

	@showprogress for ep ∈ 1:epochs
		for i ∈ 1:atom_number
			PE = Forces(time_1, w, returnPE = true)
			x[i] = x[i] + h * (vx[i] + h_half * fx[time_1, i])
			y[i] = y[i] + h * (vy[i] + h_half * fy[time_1, i])
			if x[i] <= 0.0
				x[i] = x[i] + L            # Periodic boundary conditions
			end
			if x[i] >= L
				x[i] = x[i] - L
			end
			if y[i] <= 0.0
				y[i] = y[i] + L
			end
			if y[i] >= L
				y[i] = y[i] - L
			end
		end
		PE = 0.0
		time_2 = 2
		PE = Forces(time_2, w, returnPE = true)
		KE = 0.0
		w = 0.0
		for i ∈ 1:atom_number
			vx[i] = vx[i] + h_half * (fx[time_1, i] + fx[time_2, i])
			vy[i] = vy[i] + h_half * (fy[time_1, i] + fy[time_2, i])
			KE = KE + (vx[i] * vx[i] + vy[i] * vy[i]) / 2
		end
		w = Forces(time_2, w, returnPE = false)
		P = density * (KE + w)
		T = KE / (atom_number)
		T_states[ep] = T
		P_states[ep] = P
		time += 1
		if ep % 500 == 0 || ep == 1
			plot_distribution("distribution_" * string(ep) * ".png")
		end
	end

	# plot_states(T_states, P_states, time, "MD2D_result.png")
end

println(" === experiment start === ")
@time timevolution()
