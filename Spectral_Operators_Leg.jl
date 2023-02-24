export LegTransform_gen, RadauTransform_gen, RadauBackwardTransform_gen, ∫φ_Leg

using LegendrePolynomials
using FastGaussQuadrature

function generate_nodes_weights_leg(n)


    nodes, weights = gausslegendre( n )

    return nodes, weights

end

function generate_Λ_matrix(n, nodes)

    Λ = zeros((n,n))

    for i=1:n
        for j=1:n
            Λ[i,j] = Plm(i-1,0,nodes[j])
        end
    end

    return Λ

end



# function leg_trans(φ, Λ, weights, normLeg_poly)
#
#     return normLeg_poly.*( LT*(weights.*φ) )
#
# end


function leg_inv(φ̂, Λ)

    return transpose(Λ)*φ̂
end


function LegTransform_gen(n)

    normLeg_poly = 0.5:1.0:(n - 0.5)
    nodes, weights = gaussradau(n)
    Λ = generate_Λ_matrix(n, nodes)
    legtrans(φ) = normLeg_poly.*( Λ*(weights.*φ) )
    legtrans_inv(φ̂) = transpose(Λ)*φ̂
    return legtrans, legtrans_inv, nodes, weights

end

function RadauTransform_gen(n)
    #the grid includes the initial point

    normLeg_poly = 0.5:1.0:(n - 0.5)
    nodes, weights = gaussradau(n)
    Λ = generate_Λ_matrix(n, nodes)
    radtrans(φ) = normLeg_poly.*( Λ*(weights.*φ) )
    radtrans_inv(φ̂) = transpose(Λ)*φ̂
    return radtrans, radtrans_inv, nodes, weights

end

function RadauBackwardTransform_gen(n)
    #the grid includes the end point

    normLeg_poly = 0.5:1.0:(n - 0.5)
    nodes, weights = gaussradau(n)
    nodes = -reverse(nodes)
    weights = reverse(weights)
    Λ = generate_Λ_matrix(n, nodes)
    radtrans(φ) = normLeg_poly.*( Λ*(weights.*φ) )
    radtrans_inv(φ̂) = transpose(Λ)*φ̂
    return radtrans, radtrans_inv, nodes, weights

end


function integrator_legspace(φ̂)

    n = size(φ̂)[1]
    ∫φ̂ = zeros((n,1))


    ∫φ̂[1] = φ̂[1] - φ̂[2]/3

    ∫φ̂[2] =  φ̂[1]-φ̂[3]/5

    ∫φ̂[n] = φ̂[n-1]/(2*(n-1) -1)

    for j=2:(n-2)

        ∫φ̂[j+1] =  φ̂[j]/(2*j - 1) - φ̂[j+2]/(2*j+3)

    end

    return ∫φ̂

end


function δ₁_leg(φ̂)

    return sum(φ̂)

end



function ∫φ_Leg(φ, legtrans, legtrans_inv, L)

    #L is the length of the domain
    φ̂ = legtrans(φ)
    ∫φ̂ = integrator_legspace(φ̂)*(L/2)

    return legtrans_inv(∫φ̂)

end

function ∫ₐᵇφ(φ, weights, L)

    return sum(φ .* weights)*(L/2)
end




# R, R⁻¹, nodesR, weightsR = RadauBackwardTransform_gen(n)
