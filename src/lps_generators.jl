function quadruples_4k_plus1(p::Integer)
    # @assert p % 4 == 1
    # @assert isprime(p)
    N = floor(Int, sqrt(p))
    N = iseven(N) ? N : N + 1
    quads = NTuple{4, Int}[]
    for a = 1:2:N
        s1 = a^2
        for b = -N:2:N
            s2 = s1 + b^2
            for c = -N:2:N
                s3 = s2 + c^2
                for d = -N:2:N
                    s4 = s3 + d^2
                    s4 == p && push!(quads, (a, b, c, d))
                end
            end
        end
    end
    return quads
end

function quadruples_4k_plus3(p::Integer)
    # @assert p % 4 == 3
    # @assert isprime(p)
    N = floor(Int, sqrt(p))
    N = iseven(N) ? N + 1 : N
    quads = NTuple{4, Int}[]
    for a = 0:2:N
        s1 = a^2
        b_range = (a == 0 ? StepRange(1, 2, N) : StepRange(-N, 2, N))
        for b in b_range
            s2 = s1 + b^2
            for c = -N:2:N
                s3 = s2 + c^2
                for d = -N:2:N
                    s4 = s3 + d^2
                    s4 == p && push!(quads, (a, b, c, d))
                end
            end
        end
    end
    return quads
end

function quadruples(p::Integer)
    @assert p > 0
    @assert isprime(p)

    if p % 4 == 1
        return quadruples_4k_plus1(p)
    elseif p % 4 == 3
        return quadruples_4k_plus3(p)
    end
end

# Proposition 2.5.2 of Davidoff et al.
# do not use commas or the PGL type cast won't work for a vector
generator(a₀, a₁, a₂, a₃, x, y) = [
    a₀ + x * a₁ + a₃ * y a₂ + x * a₃ - a₁ * y
    -a₂ + x * a₃ - a₁ * y a₀ - x * a₁ - y * a₃
]

function get_roots(p::Integer, q::Integer)
    if p % 4 == 1 && q % 4 == 1
        x = sqrt(GF{q}(q - 1))
        y = GF{q}(0)
    elseif p % 4 == 3 && q % 4 == 1
        # storing a in (i, i^2) pairs and doing argmax will also give x, y but this isn't a bottleneck
        a = maximum([i^2 % q for i in 1:q - 1])
        x = sqrt(GF{q}(a))
        y = sqrt(GF{q}(-a - 1))
        @assert x^2 + y^1 + 1 % q == 0
    else
        throw(DomainError("Haven't yet figured out this case."))
    end
    return x, y
end

function PGLtype(p::Integer, q::Integer)
    legendre = legendresymbol(p, q)

    if legendre == -1
        return PGL₂{q}
    elseif legendre == 1
        return PSL₂{q}
    else
        throw("legendresymbol(p, q) = $legendre")
    end
end

function lps_generators(p::Integer, q::Integer)
    # q must be an odd prime power in order for the roots to exist
    @assert p > 2
    @assert q > 2
    @assert p ≠ q
    @assert isprime(p)
    @assert isprime(q)

    x, y = get_roots(p, q)
    mats = [generator(a₀, a₁, a₂, a₃, x, y) for (a₀, a₁, a₂, a₃) in quadruples(p)]
    GL_t = PGLtype(p, q)
    S = GL_t.(mats)
    S = unique([S; inv.(S)])
    @assert all(inv(s) in S for s in S)
    @assert all((!isone).(S))
    return S
end

function lps(p::Integer, q::Integer)
    S = lps_generators(p, q)
    G, verts, vlabels, elabels = cayley_graph(order(eltype(S)), S)
    @assert order(eltype(S)) == length(verts)
    @assert all(isequal(p + 1), Graphs.degree(G))
    return G, verts, vlabels, elabels
end

function lps(p::Integer, q::Integer, radius::Integer)
    S = lps_generators(p, q)
    G, verts, vlabels, elabels = cayley_graph(S, radius = radius)
    radius < diameter_ub(p, q) &&
        @warn "Radius given is smaller than its upper bound, Cayley graph might be not complete!"
    return G, verts, vlabels, elabels
end

# from Lubotzky-Phillips-Sarnak
function diameter_ub(p::Integer, q::Integer)
    n = order(PGLtype(p, q))
    return floor(Int, 2 * log(p, n) + 2 * log(p, 2) + 1)
end
