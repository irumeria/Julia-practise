
global N = 100
global M = 101

function energy(path)
    sums = 0.0
    i = 1
    for _ in range(1, N - 2)
        sums += (path[i+1] - path[i]) ^ 2
        i += 1
    end
    sums += path[i+1] * path[i+1]
    sums
end

function main()
    path = zeros(Float32, M)
    oldE = energy(path)
    while true
        element = floor(Int, N * rand() + 1)
        change = 2.0*(rand() - 0.5)    
        path[element] += change    
        newE = energy(path)
        if newE > oldE && exp(- newE + oldE ) <= rand()
            path[element] += change  
        end
        oldE = newE
        println(oldE)
        # sleep(1)
    end
end

main()