include("edits/Substitution.jl")
include("edits/Deletion.jl")
include("edits/Insertion.jl")

"""
    Edit{S <: BioSequence, T <: BioSymbol}

An edit of either `Substitution{T}`, `Insertion{S}` or `Deletion` at a position.
If deletion: Deletion of length L at ref pos `pos:pos+L-1`
If insertion: Insertion of length L b/w ref pos `pos:pos+1`
"""
abstract type Edit{S<:BioSequence,T<:BioSymbol} end

Base.length(e::Edit) = error("length not implemented for type $(typeof(e))")
Base.:(==)(e1::Edit, e2::Edit) = typeof(e1) == typeof(e2) && e1 == e2
Base.hash(x::Edit, h::UInt) = error("hash not implemented for type $(typeof(x))")

function Base.isless(x::Edit, y::Edit)
    if leftposition(x) == leftposition(y)
        return length(x) < length(y)
    end

    return leftposition(x) < leftposition(y)
end

function BioGenerics.leftposition(e::Edit)
    return error("leftposition not implemented for type $(typeof(e))")
end
function BioGenerics.rightposition(e::Edit)
    return error("rightposition not implemented for type $(typeof(e))")
end

"""
    _lendiff(edit::Edit)

Gets the number of bases that `edit` adds to the reference sequence
"""
_lendiff(edit::Edit) = error("_lendiff not implemented for type $(typeof(edit))")

struct DeletionEdit{S,T} <: Edit{S,T}
    position::UInt
    length::UInt

    function DeletionEdit{S,T}(
        position::UInt, length::UInt
    ) where {S<:BioSequence,T<:BioSymbol}
        iszero(position) && error("Deletion cannot be at a position outside the sequence")
        iszero(length) && error("Deletion must be at least 1 symbol")
        return new(position, length)
    end
end

DeletionEdit(x::Integer, y::Integer) = DeletionEdit(convert(UInt, x), convert(UInt, y))

Base.length(d::DeletionEdit) = Int(d.length)
function Base.:(==)(d1::DeletionEdit, d2::DeletionEdit)
    return d1.position == d2.position && d1.length == d2.length
end
Base.hash(x::DeletionEdit, h::UInt) = hash(DeletionEdit, hash((x.position, x.length), h))
BioGenerics.leftposition(d::DeletionEdit) = d.position
BioGenerics.rightposition(d::DeletionEdit) = leftposition(d) + length(d) - 1
_lendiff(d::DeletionEdit) = -1 * length(d)

struct InsertionEdit{S,T} <: Edit{S,T}
    position::UInt
    seq::S

    function InsertionEdit{S,T}(position::UInt, seq::S) where {S<:BioSequence,T<:BioSymbol}
        iszero(position) && error("Insertion cannot be at a position outside the sequence")
        isempty(seq) && error("Insertion must be at least 1 symbol")
        return new(position, seq)
    end
end
Base.length(i::InsertionEdit) = length(i.seq)
Base.:(==)(x::InsertionEdit, y::InsertionEdit) = x.position == y.position && x.seq == y.seq
Base.hash(x::InsertionEdit, h::UInt) = hash(InsertionEdit, hash((x.position, x.seq), h))
BioGenerics.leftposition(i::InsertionEdit) = i.position
BioGenerics.rightposition(i::InsertionEdit) = leftposition(i) + 1
_lendiff(i::InsertionEdit) = length(i)

struct SubstitutionEdit{S,T} <: Edit{S,T}
    position::UInt
    base::T

    function SubstitutionEdit{S,T}(
        position::UInt, base::T
    ) where {S<:BioSequence,T<:BioSymbol}
        iszero(position) &&
            error("Substitution cannot be at a position outside the sequence")
        return new(position, base)
    end
end
Base.length(s::SubstitutionEdit) = 1
function Base.:(==)(x::SubstitutionEdit, y::SubstitutionEdit)
    return x.position == y.position && x.base == y.base
end
function Base.hash(x::SubstitutionEdit, h::UInt)
    return hash(SubstitutionEdit, hash((x.position, x.base), h))
end
BioGenerics.leftposition(s::SubstitutionEdit) = s.position
BioGenerics.rightposition(s::SubstitutionEdit) = leftposition(s)
_lendiff(s::SubstitutionEdit) = 0

function Base.parse(::Type{T}, s::AbstractString) where {T<:Edit{Se,Sy}} where {Se,Sy}
    return parse(T, String(s))
end

function Base.parse(::Type{<:Edit{Se,Sy}}, s::Union{String,SubString{String}}) where {Se,Sy}
    # Either "Δ1-2", "11T" or "G16C"
    if (m = match(r"^Δ(\d+)-(\d+)$", s); m) !== nothing
        pos = parse(UInt, m[1])
        stop = parse(UInt, m[2])
        stop ≥ pos || throw(ArgumentError("Non-positive deletion length: \"" * s * "\""))
        return DeletionEdit(pos, stop - pos + 1)
    elseif (m = match(r"^(\d+)([A-Za-z]+)$", s); m) !== nothing
        pos = parse(UInt, m[1])
        seq = Se(m[2])
        return InsertionEdit{Se}(pos, seq)
    elseif (m = match(r"^[A-Za-z](\d+)([A-Za-z])$", s); m) !== nothing
        pos = parse(UInt, m[1])
        sym = Sy(first(m[2]))
        return SubstitutionEdit{Sy}(pos, sym)
    else
        throw(ArgumentError("Failed to parse edit \"" * s * '"'))
    end
end
