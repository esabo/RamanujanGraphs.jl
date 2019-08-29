abstract type GL₂{q} <: AbstractMatrix{IntMod} end

Base.size(::GL₂) = (2,2)
Base.length(::GL₂) = 4

Base.one(::Type{T}) where T <: GL₂ = T(1, 0, 0, 1)
Base.one(::T) where T <: GL₂ = one(T)

Base.IndexStyle(::GL₂) = IndexLinear()

function Base.getindex(m::GL₂, i::Integer)
    i == 1 && return m.a
    i == 2 && return m.c # julia assumes column wise storage for IndexLinear
    i == 3 && return m.b
    i == 4 && return m.d
end

function Base.setindex!(m::GL₂{q}, v, i::Integer) where q
    if i == 1
        m.a = v
    elseif i == 2
        m.c = v
    elseif i == 3
        m.b = v
    elseif i == 4
        m.d = v
    end
    return v
end

ints(m::GL₂) = Int(m[1]), Int(m[2]), Int(m[3]), Int(m[4])

function Base.:(==)(m::T, n::T) where T <: GL₂
    m = normalform!(m)
    n = normalform!(n)
    return ints(m) == ints(n)
end

function Base.hash(m::T, h::UInt) where {q, T<:GL₂{q}}
    m = normalform!(m)
    a,c,b,d = ints(m)
    val = d + q*(b + q*(c + q^3*a)) # q-adic expression of m
    return hash(T, hash(val, h))
end

function Base.:(*)(m::T, n::T) where T <: GL₂
    a,c,b,d = ints(m)
    A,C,B,D = ints(n)
    new_a = a*A + b*C
    new_b = a*B + b*D
    new_c = c*A + d*C
    new_d = c*B + d*D
    return T(new_a, new_c, new_b, new_d)
end

function mul!(x::Number, m::T) where T<:GL₂
    m[1] = x*Int(m[1])
    m[2] = x*Int(m[2])
    m[3] = x*Int(m[3])
    m[4] = x*Int(m[4])
    return m
end

LinearAlgebra.det(m::GL₂{q}) where q = ((a,c,b,d) = ints(m); IntMod{q}(a*d-c*b))

function Base.inv(m::T) where {q, T <: GL₂{q}}
    a,c,b,d = ints(m)
    p = a*d - b*c
    p¯¹ = invmod(p, q)

    return T(p¯¹*d, -p¯¹*c+q, -p¯¹*b+q, p¯¹*a)
end



############################################
#
# PGL₂{q}
#
############################################

mutable struct PGL₂{q} <: GL₂{q}
    a::IntMod{q}
    c::IntMod{q}
    b::IntMod{q}
    d::IntMod{q}

    function PGL₂{q}(a,c,b,d) where q
        @assert q isa Integer
        q > 1 || error(ArgumentError("$q (the modulus) must be > 1"))
        m = new{q}(a,c,b,d)
        m = normalform!(m)
        det(m) ≠ 0 || throw(ArgumentError("Singular Matrix in PGL₂{$q}: $m"))
        return m
    end

    PGL₂{q}(m::AbstractMatrix) where q = PGL₂{q}(m[1,1], m[2,1], m[1,2], m[2,2])
end

isnormal(m::PGL₂) = isone(m[1]) || (iszero(m[1]) && isone(m[2]))

function normalform!(m::PGL₂)
    isnormal(m) && return m
    if !iszero(m[1])
        a = inv(m[1])
    else
        a = inv(m[2])
    end
    m = mul!(Int(a), m)
    return m
end

order(::Type{PGL₂{q}}) where q = q^3 - q

############################################
#
# PSL₂{q}
#
############################################
mutable struct PSL₂{q} <: GL₂{q}
    a::IntMod{q}
    c::IntMod{q}
    b::IntMod{q}
    d::IntMod{q}

    function PSL₂{q}(a,c,b,d) where q
        @assert q isa Integer
        q > 1 || error(ArgumentError("$q (the modulus) must be > 1"))
        m = new{q}(a,c,b,d)
        m = normalform!(m)
        # det(m) == 1 || throw(ArgumentError("Matrix of determinant ≠ 1 in PGL₂{$q}: $m"))
        return m
    end

    PSL₂{q}(m::AbstractMatrix) where q = PSL₂{q}(m[1,1], m[2,1], m[1,2], m[2,2])
end

function isnormal(m::PSL₂{q}) where q
    isone(det(m)) || return false
    if 0 < m[1] <= div(q-1,2)
        return true
    elseif iszero(m[1]) && 0 < m[2] <= div(q-1,2)
        return true
    else
        return false
    end
end

function normalform!(m::PSL₂{q}) where q
    isnormal(m) && return m
    p = det(m)
    if isone(p)
        xinv = 1
    else
        x = sqrt(p)
        xinv = Int(inv(x))
    end

    elt = ifelse(iszero(m[1]), m[2], m[1])

    if xinv*elt > div(q-1, 2)
        xinv = q - xinv
    end

    m = mul!(xinv, m)

    return m
end

order(::Type{PSL₂{q}}) where q = div(q^3 - q, 2)
