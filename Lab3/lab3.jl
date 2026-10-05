using LinearAlgebra
using Random

n = 100
eps = 1e-6
max_iter = 10000

Random.seed!(42)

function generate_system(n)
    A = Float64.(rand(-10:10, n, n))
    for i in 1:n
        s = sum(abs.(A[i, :])) - abs(A[i, i])
        A[i, i] = s + rand(1:10)
    end
    b = Float64.(rand(-10:10, n))
    return A, b
end

function simple_iteration(A, b, eps, max_iter)
    n = length(b)
    P = zeros(n, n)
    g = zeros(n)
    for i in 1:n
        g[i] = b[i] / A[i, i]
        for j in 1:n
            if i != j
                P[i, j] = -A[i, j] / A[i, i]
            end
        end
    end

    # println("Матрица P:")
    # display(P)
    
    P_norm = opnorm(P, Inf)
    if P_norm < 1
        println("Норма P = ", P_norm)
    else
        println("Норма P = ", P_norm)
    end
    
    # println("Вектор g:")
    # display(g)
    
    x_prev = zeros(n)
    for k in 1:max_iter
        x = P * x_prev + g
        delta_x = x - x_prev
        error = norm(delta_x)
        if error <= eps
            return x, k, error
        end
        x_prev = x
    end
    error("Метод не сошелся за $max_iter итераций")
end

A, b = generate_system(n)
x, iterations, err = simple_iteration(A, b, eps, max_iter)

println("Количество итераций: ", iterations)
println("x = ", x)
println("delta_x^(k): ", err)

x_check = A \ b
println("Решение из Ax=b: ", x_check)
println("Разница: ", x - x_check)
println("Разница: ", norm(x - x_check))
