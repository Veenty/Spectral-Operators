push!(LOAD_PATH, "/Users/Veenty/Documents/JuliaPackages/Spectral Operators")
using SpectralOperators
using FFTW
using Plots




n = 20
L = 1.0
N = 2*n

fftn = plan_fft(zeros(n,n),flags = FFTW.MEASURE)
fft⁻¹N = plan_ifft(zeros(N,N),flags = FFTW.MEASURE)
fftN = plan_fft(zeros(N,N),flags = FFTW.MEASURE)
fft⁻¹n = plan_ifft(zeros(n,n),flags = FFTW.MEASURE)

h = L/n
H = 1/N
x = h*collect(0:n-1)
X = H*collect(0:N-1)



f(x,y)  = exp(cos(2*π*x/L)* sin(2*π*y/L) )

fx = [f(x1, x2) for x1=x, x2=x]
fX = [f(x1, x2) for x1=X, x2=X]

Mult = SpectralOperators.interface_prod_Fourier(fx, fx, fftN, fftn,fft⁻¹N, fft⁻¹n)

maximum(abs.(Mult - fx.*fx))
