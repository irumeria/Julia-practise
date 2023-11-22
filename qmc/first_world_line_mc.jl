using ProgressMeter

mutable struct SystemStateXY
	const size::Int32
	const T::Float32
	const algorithm::String
	const weight_matrix::Array{Float32, 2}
	spin_conf::Array{Int32, 2}
end

function weight_matrix_2D(delta_t::Float32, J::Float32) # (1,1) (1,2) (2,1) (2,2)
	weights = zeros(Float32, (4, 4))
	weights[1, 1] = weights[4, 4] = 1
	weights[2, 2] = weights[3, 3] = cosh(0.5 * delta_t * J)
	weights[2, 3] = weights[3, 2] = sinh(0.5 * delta_t * J)
	weights
end

function calRatioLocalWeight(
	expanded_spins::Array{Int32, 2},
	weights_matrix::Array{Float32, 2})

    

end

function updateCheckBoard(sys::SystemStateXY)
	odd_inds = [i for i ∈ 1:2:sys.size]
	last_odd_ind = odd_inds[end] - 2
	if last_odd_ind == sys.size
		last_odd_ind -= 2
	end
	for i ∈ 1:2:last_odd_ind # X axis
		for j ∈ 1:2:last_odd_ind # Y axis
			plaqette = [i j; i+1 j; i j+1; i+1 j+1]
			expanded_plaqette = zeros(Int32, (4, 4, 2)) 

			for i_plus ∈ -1:3
				for j_plus ∈ -1:3
					expanded_plaqette[i_plus + 2, j_plus + 2, :] = [i + i_plus, j + j_plus]
				end
			end
			plaqette_spins = sys.spin_conf[plaqette[:, 1], plaqette[:, 2]]
			expanded_spins = zeros(Int32, (4, 4))

			for i in 1:size(expanded_plaqette)[1]
				for j in 1:size(expanded_plaqette)[2]
					expanded_spins[i, j] = sys.spin_conf[expanded_plaqette[i, j, 1], expanded_plaqette[i, j, 2]]
				end
			end

			@show size(expanded_spins)
			if plaqette_spins[1] == plaqette_spins[3] &&
			   plaqette_spins[2] == plaqette_spins[4] &&
			   plaqette_spins[1] == -plaqette_spins[2]
				calRatioLocalWeight(expanded_spins, sys.weights_matrix)

			end
		end
	end
end


function run_experiment()

	steps = 10000
	delta_img_t = 1
	J = 1
	weight_matrix = weight_matrix_2D(delta_img_t, J)
	sys = SystemStateXY(
		100,
		300,
		"Simple_World_Line",
		weight_matrix,
		rand([-1, 1], (100, 100)),
	)

	@showprogress for step ∈ 1:steps
		updateCheckBoard(sys)
	end
end

run_experiment()
