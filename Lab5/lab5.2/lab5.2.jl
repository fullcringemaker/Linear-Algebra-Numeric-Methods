import Pkg
Pkg.add("Polynomials")

using Polynomials

A = [1.0 2.0 3.0;
     2.0 4.0 5.0;
     3.0 5.0 6.0]

eps = 1e-6

n = size(A, 1)

# n = 3
# A = zeros(n, n)

# for i in 1:n
#     for j in i:n
#         value = round(1 + 9 * rand(), digits=2)
#         A[i,j] = value
#         A[j,i] = value
#     end
# end

println("A =")
display(A)

E = zeros(n, n)
for i in 1:n
    E[i,i] = 1.0
end

D = copy(A)
B_total = copy(E)

for k in n:-1:2
    # step = n - k + 1
    B = copy(E)
    B_inv = copy(E)
    for j in 1:n
        if j == k-1
            B[k-1,j] = 1.0 / D[k,k-1]
        else
            B[k-1,j] = -D[k,j] / D[k,k-1]
        end
        B_inv[k-1,j] = D[k,j]
    end
    D = B_inv * D * B
    B_total = B_total * B
    # println("B$step=")
    # display(B)
    # println("B$step^(-1)=")
    # display(B_inv)
    # println("D$step=")
    # display(D)
end

P = D
p = collect(P[1,:])
println("\nМатрица Фробениуса:")
display(P)

println("\nКоэффициенты:")
for i in 1:n
    println("p$i=", p[i])
end

poly = Polynomial(vcat(-reverse(p), 1.0), :lambda)
println("\nХарактеристический многочлен:")
println(poly)

f(x) = poly(x)

lambda_lib = sort(real.(roots(poly)))

println("Библиотечные собственные значения:")
for i in 1:n
    println("lambda", i, " = ", lambda_lib[i])
end

left_bounds = Float64[]
right_bounds = Float64[]

println("\nКруги Гершгорина:")
for i in 1:n
    center = A[i,i]
    radius = 0.0
    for j in 1:n
        if i != j
            radius = radius + abs(A[i,j])
        end
    end
    a = center - radius
    b = center + radius
    push!(left_bounds, a)
    push!(right_bounds, b)
    println("Круг ", i)
    println("Центр: ", center)
    println("Радиус: ", radius)
    println("a = ", a)
    println("b = ", b)
    println()
end

g_left = minimum(left_bounds)
g_right = maximum(right_bounds)
println("Общие границы: [", g_left, ", ", g_right, "]")

intervals = Tuple{Float64,Float64}[]
parts = 10

for i in 1:parts
    a = g_left + (g_right - g_left)*(i-1)/parts
    b = g_left + (g_right - g_left)*i/parts
    if f(a)*f(b) < 0
        push!(intervals, (a, b))
    end
end

println("\nИнтервалы собственных значений:")
for i in 1:length(intervals)
    a, b = intervals[i]
    println("lambda", i, ": [", a, ", ", b, "]")
end

function bisection(a, b)
    if f(a)*f(b) >= 0
        error("На концах интервала функция должна иметь разные знаки")
    end
    history = []
    n = 0
    while (b - a)/2 > eps
        c = (a + b)/2
        push!(history, (a, b, c))
        n = n + 1
        if f(c) == 0
            return c, n, history
        elseif f(a)*f(c) < 0
            b = c
        else
            a = c
        end
    end
    return (a + b)/2, n, history
end

b_values = Float64[]
b_iterations = Int[]
b_histories = []

println("\nМетод половинного деления:")
for i in 1:length(intervals)
    a, b = intervals[i]
    value, count, history = bisection(a, b)
    push!(b_values, value)
    push!(b_iterations, count)
    push!(b_histories, history)
    println("Собственное значение ", i, ": ", value)
    println("Количество итераций: ", count)
end

tau = 0.005

function relaxation(x0, tau)
    history = [x0]
    n_iter = 0
    while abs(f(x0)) > eps
        x0 = x0 - tau*sign(df(x0))*f(x0)
        push!(history, x0)
        n_iter = n_iter + 1
        if n_iter > 1000
            error("Метод релаксаций не сошелся")
        end
    end
    return x0, n_iter, history
end

r_values = Float64[]
r_iterations = Int[]
r_histories = Vector{Float64}[]

println("\nМетод релаксаций:")

for i in 1:n
    a, b = intervals[i]
    x0 = (a + b)/2
    value, count, history = relaxation(x0, tau)
    push!(r_values, value)
    push!(r_iterations, count)
    push!(r_histories, history)
    println("Собственное значение ", i, ": ", value)
    println("Количество итераций: ", count)
end

function secant(x0, x1)
    history = [x0, x1]
    n = 0
    while abs(f(x1)) > eps
        x_new = x1 - f(x1)*(x1 - x0)/(f(x1) - f(x0))
        push!(history, x_new)
        x0 = x1
        x1 = x_new
        n = n + 1
    end
    return x1, n, history
end

s_values = Float64[]
s_iterations = Int[]
s_histories = Vector{Float64}[]

println("\nМетод секущих:")
for i in 1:length(intervals)
    a, b = intervals[i]
    value, count, history = secant((a+b)/2, b)
    push!(s_values, value)
    push!(s_iterations, count)
    push!(s_histories, history)
    println("Собственное значение ", i, ": ", value)
    println("Количество итераций: ", count)
end

println("\nСравнение с библиотечным решением:")

for i in 1:n
    println("\nlambda", i, ":")
    println("Точное решение: ", lambda_lib[i])
    println("Метод половинного деления:")
    println("Корень: ", b_values[i])
    println("Сравнение с точным: ", abs(lambda_lib[i] - b_values[i]))
    println("Метод релаксаций:")
    println("Корень: ", r_values[i])
    println("Сравнение с точным: ", abs(lambda_lib[i] - r_values[i]))
    println("Метод секущих:")
    println("Корень: ", s_values[i])
    println("Сравнение с точным: ", abs(lambda_lib[i] - s_values[i]))
end
