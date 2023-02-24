push!(LOAD_PATH, "/Users/Veenty/Documents/JuliaPackages/Spectral Operators")
using SpectralOperators

#example for 1d

function recursive_coefficients(a, r, Cₙ, Legendre_coefficients, N)

    ϵ = 10^(-14)
    CN = Cₙ[N]
    if abs.(CN) <=ϵ
        CN = 0.0
    end
    a⁻¹CN = (1/a)*CN
    println(r/a)
    Cₙ[N-1] = Cₙ[N-1] + (2(N-1) + 1.0)*r*a⁻¹CN
    Cₙ[N-2] = Cₙ[N-2] -  CN

    Legendre_coefficients[N] = Legendre_coefficients[N] - a⁻¹CN
    Legendre_coefficients[N-2] = Legendre_coefficients[N-2] + a⁻¹CN

    if N-1 ==2

        a⁻¹C2 = (1/a)*Cₙ[2]
        v = (r)*a⁻¹C2 + Cₙ[1]
        a⁻¹v = (1/a)*v
        Legendre_coefficients[2] =  Legendre_coefficients[2] - a⁻¹C2
        Legendre_coefficients[1] =  Legendre_coefficients[1] - a⁻¹v

        return  -a⁻¹v + v,Legendre_coefficients

    else

        return recursive_coefficients(a, r, Cₙ, Legendre_coefficients, N-1)

    end
end



N = 10
a = -2.0
convert(Float64, a)
L = 0.1
r = 2/L


radtrans, radtrans_inv, nodes, weights = LegTransform_gen(N)
nodesL = (nodes .+ 1)/r
weightsL = weights/r


f(x) = sin(x)#r*x - 1
F = f.(nodesL)

LF = radtrans(F)

LC = zeros(N)

expcoef, Leg = recursive_coefficients(a, r, copy(LF), LC, N)

term1 = radtrans_inv(Leg)

term2 = exp.(a*nodesL)*expcoef

real_solution(t) = (exp(a*t) - a*sin(t) - cos(t))/(a^2 + 1)
real_solution_debug(t) = (- a*sin(t) - cos(t))/(a^2 + 1)

term1_debug = radtrans(real_solution_debug.(nodesL) )
real_solution.(nodesL) - (term1 + term2)
