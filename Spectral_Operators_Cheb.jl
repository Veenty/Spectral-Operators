
#Generates the function ∫ₐˣ ⋅ dx
export ChebIntegratorGen, 𝘾, 𝘾⁻¹, IntegrateCoeff, δ₁, δ₋₁,δ̂₁,δ̂₋₁,DerivativeCheb
using FFTW
# using BenchmarkTools
# using Plots





function IntegrateCoeff(f̂)

   N = size(f̂)[1]
   ∫f̂ = zeros(N)

   # ∫f̂[1] = 0.0
   # ∫f̂[2] = (f̂[1] - f̂[3]/2)
   # ∫f̂[3] = (f̂[2]/4 - f̂[4]/4)

   ∫f̂[1] = 0.0
   ∫f̂[2] = (f̂[1] - f̂[3]/2)
   ∫f̂[3] = (f̂[2]/4 - f̂[4]/4)
   for k=4:(N -1)
       # ∫f̂[k+1] = (f̂[k] - f̂[k+2])/(2*k)
       ∫f̂[k] = (f̂[k-1] - f̂[k+1])/(2*(k-1))
   end

   #∫f̂[N] = f̂[N-1]/(2*(N-1))
   ∫f̂[N] = f̂[N-1]/(2*(N-1))


   return ∫f̂

end


function der_matrix(deg::Int, a::Real=-1, b::Real=1)
    N = deg
    D = zeros(Float64, N, N+1)
    for i=1:N, j=1:N+1
        if i == 1
            if iseven(i + j)
                continue
            end
            D[i, j] = 2*(j-1)/(b-a)
        else
            if j < i || iseven(i+j)
                continue
            end
            D[i, j] = 4*(j-1)/(b-a)
        end
    end

    D
end


function DerivativeCheb(f, L, DCT, DCT⁻¹, type)

   f̂ = 𝘾(f, DCT, type)

   println("testing sizes")
   println(length(f̂))
   println(size(der_matrix(length(f) - 1, 0, L)))
   df̂ = der_matrix(length(f) - 1, 0, L) * f̂
   append!(df̂ , 0)

   return 𝘾⁻¹(df̂, DCT⁻¹, type)



end


function 𝘾(f, DCT, type)

   N = size(f)[1]
   f̂ = DCT*f

   #Now we need to renormalize depending on the type

   if type == 1

      f̂[1] = f̂[1]/2
      f̂[size(f)[1]] = f̂[size(f)[1]]/2
      f̂ = f̂/(N-1)

   elseif type == 2

      f̂[1] = f̂[1]/2
      f̂ = f̂/N

   end

   return f̂


end


function 𝘾⁻¹(f̂, DCT⁻¹, type)

   N = size(f̂)[1]
   f̂_ = copy(f̂)

   #Now we need to renormalize depending on the type

   if type == 1

      f̂_[1] = f̂[1]*2 #*(N-1)
      f̂_[end] = f̂[end]*2#*(N-1)


   elseif type == 2

      f̂_[1] = f̂_[1]*2#*(N)

   end

   return (DCT⁻¹*f̂_)/2


end


function ∫indefinite(f, L, type, DCT, DCT⁻¹)

   #type corresponds to the type of cosine transform
   #options are 1 or 2
   #We still no normalize depending on the interval length, and starting point

   f̂ = 𝘾(f, DCT, type)
   ∫f̂ = IntegrateCoeff(f̂)
   ∫f = 𝘾⁻¹(f̂, DCT⁻¹, type)

   return L/2*∫f

end


function δ₁(f, DCT ,type)

   if type ==1

      return f[1]#f[end]

   else
      return δ̂₁(𝘾(f, DCT, type))

   end
end


function δ₋₁(f, DCT ,type)

   if type ==1

      return f[end]#f[1]

   else
      return δ̂₋₁(𝘾(f, DCT, type))

   end
end




function δ̂₁(f̂)

   return sum(f̂)

end


function δ̂₋₁(f̂)

   f₋₁ = 0
   N = size(f̂)[1]
   for i=1:N
      f₋₁+= f̂[i]*(-1)^(i-1)
   end

   return f₋₁


end

function Cheb_∫_Gen(N, a, b, type)

   if type == 1

      DCT = FFTW.plan_r2r(ones(N) ,FFTW.REDFT00, 1, flags = FFTW.MEASURE)
      DCT⁻¹ = DCT

   else
      #type ==2

      DCT = FFTW.plan_r2r(ones(N) ,FFTW.REDFT10, 1, flags = FFTW.MEASURE)
      DCT⁻¹ = FFTW.plan_r2r(ones(N) ,FFTW.REDFT01, 1, flags = FFTW.MEASURE)


   end

end


function def_int_left(f, DCT, L)

   N = size(f)[1]
   f̂ = DCT * f #un normalized
   f̂[2:N] = f̂[2:N]/(N-1)
   f̂[1] = f̂[1]/(2*N-2)

   #println(f̂)

   ∫f̂ = zeros(N)

   ∫f̂[1] = 0.0
   ∫f̂[2] = (f̂[1] - f̂[3]/2)
   ∫f̂[3] = (f̂[2]/4 - f̂[4]/4)
   for k=4:(N -1)
       # ∫f̂[k+1] = (f̂[k] - f̂[k+2])/(2*k)
       ∫f̂[k] = (f̂[k-1] - f̂[k+1])/(2*(k-1))
   end

   #∫f̂[N] = f̂[N-1]/(2*(N-1))
   ∫f̂[N] = f̂[N-1]/((N-1))



   #we do stuff with the coefficients

   ∫f̂[1] = ∫f̂[1]*(2*N-2)
   ∫f̂[2:N] = ∫f̂[2:N]*(N-1)

   ∫f = DCT * ∫f̂/ (2*N-2)  #This is the primitive

   ∫f = L/2*(∫f .- ∫f[end] )#Scale to the proper length of the interval

   return ∫f

end

function def_int_right(f, DCT, L)

   N = size(f)[1]
   f̂ = DCT * f #un normalized
   f̂[2:N] = f̂[2:N]/(N-1)
   f̂[1] = f̂[1]/(2*N-2)

   #println(f̂)

   ∫f̂ = zeros(N)

   ∫f̂[1] = 0.0
   ∫f̂[2] = (f̂[1] - f̂[3]/2)
   ∫f̂[3] = (f̂[2]/4 - f̂[4]/4)
   for k=4:(N -1)
       # ∫f̂[k+1] = (f̂[k] - f̂[k+2])/(2*k)
       ∫f̂[k] = (f̂[k-1] - f̂[k+1])/(2*(k-1))
   end

   #∫f̂[N] = f̂[N-1]/(2*(N-1))  #the last coef has a special normalization as well
   ∫f̂[N] = f̂[N-1]/((N-1))


   #we do stuff with the coefficients

   ∫f̂[1] = ∫f̂[1]*(2*N-2)
   ∫f̂[2:N] = ∫f̂[2:N]*(N-1)

   ∫f = DCT * ∫f̂/ (2*N-2)  #This is the primitive

   ∫f = L/2*(∫f[1] .- ∫f  )#Scale to the proper length of the interval

   return ∫f

end


function ChebIntegratorGen2(N, a, b, cumulative = true)

   precomp_dct = FFTW.plan_r2r(ones(N) ,FFTW.REDFT00, 1, flags = FFTW.MEASURE)

   if cumulative

      ∫ₐˣ(f) = def_int_left(f, precomp_dct, b-a)
      return ∫ₐˣ

   else

      ∫ₓᵇ(f) = def_int_right(f, precomp_dct, b-a)
      return ∫ₓᵇ

   end

end
#
# N = 14
# θ = 1/N*collect(0:N)
# x = cos.(π*θ)
#
#
# n = 13
# y = x .^ n
# ∫y = 1/(n+1)*x .^ (n+1) .- 1/(n+1)*x[end]^(n+1)
#
# p = FFTW.plan_r2r(y ,FFTW.REDFT00, 1, flags = FFTW.MEASURE)
#
#
# # Fk_2 = FFTW.r2r(y, FFTW.REDFT00)
# #
# # @btime FFTW.r2r($y, FFTW.REDFT00)
# # @btime $p * $y
# ∫ = ChebIntegratorGen(N +1, -1, 1,  false)
#
# a = ∫(y)
#
# plot(x, a)
# plot!(x, ∫y)
#
# error_sup = maximum(abs.(∫y - a))
#
# println(error_sup)
