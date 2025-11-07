using Enzyme

function newton_solve(f, x0; kw...)
  grad(x) = Enzyme.autodiff(Forward, f, Duplicated(x, 1.0))[1]
  return newton_solve(f, x0, grad; kw...)
end

function newton_solve(f, x0, grad; atol=1e-10, rtol=1e-10, maxiter=20, verbose=false)
  x = float(x0)
  iter = 0
  fx = 0.0

  while iter < maxiter
    verbose && println("[iter = $iter], x = $x, \t fx = $fx")
    iter += 1
    fx = f(x)
    abs(fx) <= atol && return x

    dfx = grad(x)
    Δx = fx / dfx
    x -= Δx

    isfinite(x) || (@warn("Non-finite value at iteration $iter"); return x)
    abs(Δx) <= atol + rtol * abs(x) && return x
  end
  @warn("Newton failed to converge after $iter iterations")
  return x
end


# 测试函数：f(x) = x^2 - 4
f(x) = x^2 - 4

grad(x) = Enzyme.autodiff(Forward, f, Duplicated(x, 1.0))[1]

x0 = 1.0
x_sol = newton_solve(f, x0; verbose=true)

# 输出结果
println("Computed solution: $x_sol")
println("Exact solution: 2.0")
println("Error: $(abs(x_sol - 2.0))")
