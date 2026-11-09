import Pkg
Pkg.add("Polynomials")

using Polynomials

A = [1.0 2.0 3.0;
     2.0 4.0 5.0;
     3.0 5.0 6.0]

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

lambda = sort(real.(roots(poly)))
println("\nСобственные значения: ", round.(lambda, digits=2))

println("\nB =")
display(B_total)

x = Vector{Vector{Float64}}()

for i in 1:n
    yi = [lambda[i]^(n-j) for j in 1:n]
    xi = B_total * yi
    push!(x, xi)
    println("\nlambda$i=", lambda[i])
    println("y$i=", yi)
    println("x$i=", xi)
end

println("\nПроверка ортогональности:")

for i in 1:n-1
    for j in i+1:n
        s = sum(x[i][k] * x[j][k] for k in 1:n)
        println("(x$i,x$j)=", s)
    end
end

y0 = zeros(n)
y0[1] = 1.0
y = [y0]

for k in 1:n
    push!(y, A * y[end])
end

for k in 0:n
    println("y$k=", y[k+1])
end

M = hcat(y[n:-1:1]...)
C = hcat(M, y[n+1])

println("\nСЛАУ:")
display(C)

for k in 1:n-1
    for i in k+1:n
        factor = C[i,k] / C[k,k]
        for j in k:n+1
            C[i,j] -= factor * C[k,j]
        end
    end
end

p = zeros(n)

for i in n:-1:1
    p[i] = C[i,n+1]
    for j in i+1:n
        p[i] -= C[i,j] * p[j]
    end
    p[i] /= C[i,i]
end

println("\nКоэффициенты:")
for i in 1:n
    println("p$i=", p[i])
end

poly = Polynomial(vcat(-reverse(p), 1.0), :lambda)
println("\nХарактеристический многочлен:")
println(poly)

lambda = sort(real.(roots(poly)))
println("\nСобственные значения:", round.(lambda, digits=2))

x = Vector{Vector{Float64}}()

for i in 1:n
    q = zeros(n)
    q[1] = 1.0
    for j in 2:n
        q[j] = lambda[i] * q[j-1] - p[j-1]
    end
    xi = zeros(n)
    for j in 1:n
        xi += q[j] * y[n-j+1]
    end
    push!(x, xi)
    println("\nlambda$i=", lambda[i])
    println("q$i=", q)
    println("x$i=", xi)
end

println("\nПроверка ортогональности:")

for i in 1:n-1
    for j in i+1:n
        s = sum(x[i][k] * x[j][k] for k in 1:n)
        println("(x$i,x$j)=", s)
    end
end
