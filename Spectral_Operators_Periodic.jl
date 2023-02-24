#Periodic Operators
export periodic_∂ₓ, periodic_∂ₓ², periodic_∂ₓc∂ₓ, periodic_∂ₓⁿ, generator_prod, character_mult
export interface_prod, OverSampling2D, DownSampling2D
using FFTW
# To do list
# Precompute FFT and make function generators
# Change ∂ₓc∂ₓ so it actually is SDP (to handle the nyquist term)
# Over sampling
# Downsampling


function periodic_∂ₓ(f, L)

    f̂ = fft(f)
    N = size(f)[1]
    f̂′ = zeros(ComplexF64, N)

    for k=0:(N-1)

        if k< N/2

            f̂′[k+1] = f̂[k+1]*(2*π*1im/L) *k

        elseif k >N/2

            f̂′[k+1] = f̂[k+1]*(2*π*1im/L) * (k-N)

        else

            f̂′[k+1] = 0.0
        end


    end


    return real.(ifft(f̂′))

end

function periodic_∂ₓ²(f, L)

    f̂ = fft(f)
    N = size(f)[1]
    f̂′′ = zeros(ComplexF64, N)


    for k=0:(N-1)

        if k <= N/2

            f̂′′[k+1] = -f̂[k+1]*(2*π*k/L)^2

        elseif k >N/2

            f̂′′[k+1] = -f̂[k+1]*(2*π*(k-N)/L)^2


        end


    end

    return real.(ifft(f̂′′))

end


function periodic_∂ₓⁿ(f,L,n)

    N = size(f)[1]
    derivatives_2 = n ÷ 2

    f⁽²ⁿ⁾ = zeros(N)

    for j =1:n

        f⁽²ⁿ⁾ = periodic_∂ₓ²(f⁽²ⁿ⁾, L)

    end

    if N%2 ==1

        return periodic_∂ₓ(f⁽²ⁿ⁾, L)

    else

        return f⁽²ⁿ⁾

    end


end

function periodic_∂ₓc∂ₓ(f, c, L)

    N = size(f)[1]


    if N%2 == 1
        #We dont need to consider alianising effects
        return periodic_∂ₓ(c.*periodic_∂ₓ(f, L), L)

    else

        f̂ = fft(f)
        f̂′ = zeros(ComplexF64, N)
        f̂_alianising =  f̂[N÷2 + 1]


        for k=0:(N-1)



            if k< N/2

                f̂′[k+1] = f̂[k+1]*(2*π*1im/L) *k

            elseif k >N/2

                f̂′[k+1] = f̂[k+1]*(2*π*1im/L) * (k-N)

            else

                f̂′[k+1] = 0.0
            end



        end

        f′ = real.(ifft(f̂′))

        v = c.*f′
        v̂ = fft(v)
        v̂′ = zeros(ComplexF64, N)

        for k=0:(N-1)

            if k< N/2

                v̂′[k+1] = v̂[k+1]*(2*π*1im/L) *k

            elseif k >N/2

                v̂′[k+1] = v̂[k+1]*(2*π*1im/L) * (k-N)

            else

                v̂′[k+1] = 0.0#-sum(c)/N *(π/L * N^2)*f̂_alianising
            end



        end


        return real.(ifft(v̂′))

    end

end

function OverSampling1D( f , N, fftn, fft⁻¹N)

    n = size(f)[1]
    f̂ = fftn*f

    if n%2 == 1

        f̂ = vcat(f̂[1:div(n,2)+1],zeros(N-n), f̂[div(n,2)+2:end])

    else

        f̂ₙ½ = f̂[div(n,2)+1]
        println(f̂ₙ½)
        f̂ = vcat(f̂[1:div(n,2)],[f̂ₙ½/2],zeros(N-n-1),[conj(f̂ₙ½)/2], f̂[div(n,2)+2:end])

    end

    return real.(fft⁻¹N*f̂)*(N/n)

end

function DownSampling1D( f , n, fftN, fft⁻¹n)

    N = size(f)[1]
    f̂ = fftN*f

    if n%2 == 1

        f̂ = vcat(f̂[1:div(n,2)+1],f̂[N - div(n,2) + 1:end])

    else

        f̂ₙ½ = (f̂[div(n,2)+1] + f̂[N - (div(n,2))])/2
        f̂ = vcat(f̂[1:div(n,2)],[f̂ₙ½], f̂[end - div(n,2) + 2:end])

    end

    return real.(fft⁻¹n*f̂)*(n/N)

end


function OverSampling2D(f, N, fftn, fft⁻¹N)

    #N is an array of the size to over sample to each dim
    #fftn is a precomputed d-dimensional dft
    #fft⁻¹N is a precomputed d-dimensional inv_dft


    f̂ = fftn*f
    #d = size(N)[1]
    n₁ = size(f)[1]
    n₂ = size(f)[2]


    if n₁%2 == 1

        f̂ =  vcat(f̂[1:(div(n₁,2) +1),:], zeros(N[1] - n₁, n₂), f̂[div(n₁,2)+2:end,:] )

    else

        f̂ₙ½ = f̂[div(n₁,2)+1,:]

        #aliasing term

        f̂ = vcat(f̂[1:div(n₁,2),:],transpose(f̂ₙ½)/2,zeros(N[1]-n₁-1, n₂),conj.(transpose(f̂ₙ½)/2), f̂[div(n₁,2)+2:end,:])



    end


    if n₂%2 == 1

        f̂ =  hcat(f̂[:,1:div(n₂,2)+1], zeros(N[1], N[2]- n₂),  f̂[:,div(n₂,2)+2:end] )

    else

        #aliasing term
        f̂ₙ½ = f̂[:,div(n₂,2)+1]
        f̂ =  hcat(f̂[:,1:div(n₂,2)], f̂ₙ½/2, zeros(N[1], N[2]- n₂-1), f̂ₙ½/2, f̂[:,div(n₂,2)+2:end] )

    end


    return real.((fft⁻¹N*f̂))*(N[1]*N[2])/(n₁*n₂)


end


function DownSampling2D(f, n, fftN, fftn⁻¹)


    N₁ = size(f)[1]
    N₂ = size(f)[2]
    f̂ = fftN*f

    if n[1]%2 == 1

        f̂ = vcat(f̂[1:div(n[1],2)+1,:],f̂[N₁ - div(n[1],2) + 1:end,:])

    else


        f̂ₙ½ = (f̂[div(n[1],2)+1, :] + conj(f̂[N₁ - div(n[1],2),:]) )/2
        f̂ = vcat(f̂[1:div(n[1],2),:],transpose(f̂ₙ½), f̂[end - div(n[1],2)+2:end,:])

    end


    if n[2]%2 == 1

        f̂ = hcat(f̂[:,1:div(n[2],2)+1],f̂[:,N₂ - div(n[2],2)+1:end])

    else

        f̂ₙ½ = (f̂[:,div(n[2],2)+1, :] + conj(f̂[:,N₂ - div(n[2],2),:]) )/2
        f̂ = hcat(f̂[:,1:div(n[2],2)],f̂ₙ½, f̂[:,end - div(n[2],2)+2:end])

    end

    return real.(fftn⁻¹*f̂)*(n[1]*n[2])/(N₁*N₂)


end

function interface_prod(w, v, fftN, fftn,fftN⁻¹, fftn⁻¹)

    n = size(w)
    N = 2 .*n
    w_over = OverSampling2D(w, N, fftn, fftN⁻¹)
    v_over = OverSampling2D(v, N, fftn, fftN⁻¹)

    return DownSampling2D(w_over .* v_over, n, fftN, fftn⁻¹)

end



function interface_prod_Fourier(w, v, fftN, fftn,fftN⁻¹, fftn⁻¹)

    n = size(w)
    N = 2 .*n
    w_over = OverSampling2D(w, N, fftn, fftN⁻¹)
    v_over = OverSampling2D(v, N, fftn, fftN⁻¹)

    return DownSampling2D(w_over .* v_over, n, fftN, fftn⁻¹)

end

function  generator_prod(n)
    #n are dimensions
    N = 2 .*n

    fftn = plan_fft(zeros(n),flags = FFTW.MEASURE)
    fftN⁻¹ = plan_ifft(zeros(N),flags = FFTW.MEASURE)

    fftN = plan_fft(zeros(N),flags = FFTW.MEASURE)
    fftn⁻¹ = plan_ifft(zeros(n),flags = FFTW.MEASURE)

    prod(w,v) = interface_prod(w,v, fftN, fftn,fftN⁻¹, fftn⁻¹)

    return prod

end


function character_mult(f, ĉ, fftn, fftn⁻¹)


    f̂ = fftn*f

    return fftn⁻¹*(ĉ .* f̂)



end









# N = 10
# L = 1.0
#
# f = cos.((0:N-1)*2*π/N)
# c = sin.((0:N-1)*2*π/N)
#
# f′ = periodic_∂ₓ(f, L )
#
# result = - 8*π^2 *f.*c
#
# calc = periodic_∂ₓc∂ₓ(f, c, 1.0)
#
# error = abs.(result - calc)
