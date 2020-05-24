struct GF{q} <: Number
    value::Int16

    function GF{q}(n) where q
        @assert q > 1
        k = (-q<n<q ? n : n%q)
        k = ifelse(k >= 0, k, k + q)
        return new{q}(k)
    end
end

int(n::GF) = n.value
characteristic(::Type{GF{q}}) where q = q
characteristic(::GF{q}) where q = q

Base.:(==)(n::GF{q}, m::GF{q}) where q = int(n) == int(m)
# hash(RamanujanGraphs.GF) == 0x04fd9e474909f8bf
Base.hash(n::GF{q}, h::UInt) where q = xor(0x04fd9e474909f8bf, hash(q, hash(int(n), h)))

Base.:+(n::GF{q}, m::GF{q}) where q = GF{q}(int(n) + int(m))
Base.:-(n::GF{q}, m::GF{q}) where q = GF{q}(int(n) - int(m))
Base.:*(n::GF{q}, m::GF{q}) where q = GF{q}(int(n) * int(m))
Base.:/(n::GF{q}, m::GF{q}) where q = n*inv(m)

Base.:-(n::GF{q}) where q = GF{q}(q - int(n))
Base.inv(n::GF{q}) where q = GF{q}(invmod(int(n), q))

function Base.:^(n::GF{q}, i::Integer) where q
   i < 0 && return inv(n)^-i
   return GF{q}(powermod(int(n), i, q))
end

Base.zero(::Type{GF{q}}) where q = GF{q}(0)
Base.one(::Type{GF{q}}) where q = GF{q}(1)
Base.iszero(n::GF) = int(n) == 0
Base.isone(n::GF) = int(n) == 1

Base.promote_rule(::Type{GF{q}}, ::Type{I}) where {q, I<:Integer} = GF{q}

# taken from ValidatedNumerics, under under the MIT "Expat" License:
# https://github.com/JuliaIntervals/ValidatedNumerics.jl/blob/master/LICENSE.md
function subscriptify(n::Integer)
    subscript_0 = Int(0x2080) # Char(0x2080) -> subscript 0
    return join((Char(subscript_0 + d) for d in reverse(digits(n))), "")
end

Base.show(io::IO, n::GF{q}) where q = print(io, "$(int(n))"*subscriptify(q))

import Base: <, <=
for ord in [:<, :(<=)]
    @eval begin
        $ord(n::GF{q}, m::GF{q}) where q = $ord(int(n),int(m))
        $ord(n::GF, y::Number) = $ord(int(n), y)
        $ord(y::Number, n::GF) = $ord(y, int(n))
    end
end

function legendresymbol(n, q)
    iszero(mod(n, q)) && return zero(n)
    isone(powermod(n, (q-1)÷2, q)) && return one(n)
    return -one(n)
end

function generator(n::GF{q}) where q
    for i in 2:q-1
        isone(-legendresymbol(i, q)) && return GF{q}(i)
    end
    return zero(n) # never hit, to keep compiler happy
end

Base.sqrt(n::GF{q}) where q = GF{q}(sqrtmod(int(n), q))
issqrt(n::GF{q}) where q = legendresymbol(int(n), q) >= 0

function sqrtmod(n::Integer, q::Integer)
    l = legendresymbol(n, q)
    l == 0 && return zero(n)
    l == -1 && throw(DomainError(n, "$n is not a square modulo $q"))
    for i in 1:q # bruteforce loop
        y = powermod(i, 2, q)
        y == n && return oftype(n, i)
    end
    return zero(n) # never hit, to keep compiler happy
end

order(::Type{GF{q}}) where q = q
Base.iterate(::Type{GF{q}}, s=1) where q = s > q ? nothing : (GF{q}(s), s+1)
Base.eltype(::Type{GF{q}}) where q = GF{q}
Base.size(gf::Type{<:GF}) = (order(gf),)
