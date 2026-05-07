using LinearAlgebra

function main(n::Int=10000, iters::Int=10)
    total_time = 0.0

    for iter in 1:iters
        A = rand(n, n)
        t = time()
        F = lu(A)
        total_time += time() - t
    end

    println("Average time: ", total_time / iters, " seconds")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end