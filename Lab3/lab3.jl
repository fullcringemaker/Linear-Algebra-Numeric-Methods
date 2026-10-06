using LinearAlgebra
using Random

Random.seed!(42)

n = 100
k = 3
eps = 1e-6
max_iter = 10000

function generate_system(n, k)
    A = rand(1.0:100.0, n, n)
    for i in 1:n
        S = sum(abs, A[i, :]) - abs(A[i, i])
        A[i, i] = round(k * S, digits=2)
    end
    return A
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
    P_norm = opnorm(P, Inf)
    println("Норма P = ", P_norm)
    if P_norm < 1
        println("Условие сходимости ||P|| < 1 выполнено")
    else
        println("Условие сходимости ||P|| < 1 не выполнено")
    end
    x_prev = zeros(n)
    for iteration in 1:max_iter
        x = P * x_prev + g
        delta_x = x - x_prev
        err = norm(delta_x)
        if err <= eps
            return x, iteration, err
        end
        x_prev = x
    end
    error("Метод не сошелся за $max_iter итераций")
end

A = generate_system(n, k)
x_exact = ones(n)
b = A * x_exact
x, iterations, err = simple_iteration(A, b, eps, max_iter)

println("k = ", k)
println("Количество итераций: ", iterations)
println("Полученный x: ", x)
println("delta_x^(k): ", err)
println("Норма разницы: ", norm(x - x_exact))

function seidel(A, b, eps, max_iter)
    n = length(b)
    x_prev = zeros(n)
    for iteration in 1:max_iter
        x = copy(x_prev)
        for i in 1:n
            sum1 = 0.0
            sum2 = 0.0
            for j in 1:i-1
                sum1 = sum1 + A[i, j] * x[j]
            end
            for j in i+1:n
                sum2 = sum2 + A[i, j] * x_prev[j]
            end
            x[i] = (b[i] - sum1 - sum2) / A[i, i]
        end
        delta_x = x - x_prev
        err = norm(delta_x)
        if err <= eps
            return x, iteration, err
        end
        x_prev = x
    end
    error("Метод не сошелся за $max_iter итераций")
end

x_seidel, iterations_seidel, err_seidel = seidel(A, b, eps, max_iter)

println("k = ", k)
println("Количество итераций: ", iterations_seidel)
println("Полученный x: ", x_seidel)
println("delta_x^(k): ", err_seidel)
println("Норма разницы: ", norm(x_seidel - x_exact))


