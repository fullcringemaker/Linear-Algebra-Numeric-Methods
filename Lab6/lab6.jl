using LinearAlgebra, Plots

n = 100
epsilon = 1e-8
A = zeros(n, n)

for row in 1:n
    for col in row:n
        value = round(1 + 9 * rand(); digits=2)
        A[row, col] = value
        A[col, row] = value
    end
    A[row, row] = sum(A[row, :])
end

for row in 1:n
    A[row, row] = sum(A[row, :])
end

x_star = ones(n)
x0 = zeros(n)
f = A * x_star
println("Матрица A:")
display(A)

lambda = eigvals(Symmetric(A))
lambda_min, lambda_max = extrema(lambda)

tau_limit = 2 / lambda_max
tau_opt = 2 / (lambda_min + lambda_max)
q0 = (lambda_max - lambda_min) / (lambda_max + lambda_min)

println("Интервал собственных значений: [", lambda_min, ", ", lambda_max, "]")
println("Оптимальное tau = ", tau_opt)
println("q0 = ", q0)

tau_step = min(tau_opt, tau_limit - tau_opt) / 3
tau_values = [tau_opt - 2 * tau_step, tau_opt - tau_step, tau_opt, tau_opt + tau_step, tau_opt + 2 * tau_step]
q_values = [max(abs(1 - tau * lambda_min), abs(1 - tau * lambda_max)) for tau in tau_values]

iteration_counts = Int[]
error_histories = Vector{Float64}[]

for tau in tau_values
    x = copy(x0)
    k = 0
    error_norm = norm(x - x_star)
    errors = [error_norm]
    while error_norm > epsilon && k < 10000
        x = x - tau * (A * x - f)
        k += 1
        error_norm = norm(x - x_star)
        push!(errors, error_norm)
    end
    push!(iteration_counts, k)
    push!(error_histories, errors)

    println("tau = ", tau)
    println("Количество итераций: ", k)
    println("Итоговая ошибка: ", error_norm)
    println("Левое расстояние до нуля: ", abs(1 - tau * lambda_min))
    println("Правое расстояние до нуля: ", abs(1 - tau * lambda_max))
    println()
end

geometry = plot(
    [0, lambda_max, lambda_max, 0],
    [1, 1, -1, 1],
    label="Граница сектора",
    color=:black,
    xlabel="lambda",
    ylabel="mu",
    title="Геометрическая интерпретация"
)

hline!(
    geometry,
    [0],
    label=false,
    color=:gray
)

plot!(
    geometry,
    [0, lambda_max],
    [1, -q0],
    label="mu=1-tau_opt*lambda",
    color=:blue
)

plot!(
    geometry,
    [lambda_min, lambda_min],
    [0, q0],
    label="Расстояние q0",
    color=:red,
    linestyle=:dash,
    linewidth=2
)

plot!(
    geometry,
    [lambda_max, lambda_max],
    [-q0, 0],
    label=false,
    color=:red,
    linestyle=:dash,
    linewidth=2
)

scatter!(
    geometry,
    [lambda_min, lambda_max],
    [q0, -q0],
    label=false,
    color=:red
)

comparison = plot(
    tau_values,
    q_values,
    marker=:circle,
    label="Рассматриваемые значения",
    xlabel="tau",
    ylabel="q(tau)",
    title="Зависимость q от tau"
)

scatter!(
    comparison,
    [tau_opt],
    [q0],
    color=:red,
    label="Минимум q0"
)

plot(
    tau_values,
    iteration_counts,
    marker=:circle,
    legend=false,
    xlabel="tau",
    ylabel="Количество итераций",
    title="Зависимость количества итераций от tau"
)
