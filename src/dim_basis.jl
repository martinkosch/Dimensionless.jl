# Helper Unions
QuantityOrUnits = Union{Unitful.AbstractQuantity, Unitful.Units}
QuantityOrUnitlike = Union{Unitful.AbstractQuantity, Unitful.Unitlike}

"""
    DimBasis(basis_vectors...) -> DimBasis

Create a dimensional basis for a number of `basis_vectors` (quantities, units or dimensions).
A string identifier can optionally be added to each basis vector. 
"""
struct DimBasis{T,N}
    basis_vectors::T
    basis_vector_names::N
    basis_dims
    dim_mat

    function DimBasis(basis_vectors::T, basis_vector_names::N=nothing) where
        T<:AbstractVector{<:QuantityOrUnitlike} where
        N<:Union{Nothing,AbstractVector{<:AbstractString}}
        basis_dims = unique_dims(basis_vectors...)
        dim_mat = dim_matrix(basis_dims, basis_vectors...)
        check_basis(dim_mat)
        return new{T,N}(basis_vectors, basis_vector_names, basis_dims, dim_mat)
    end
end

DimBasis(basis_vectors::Vararg{QuantityOrUnitlike}) =
DimBasis([basis_vectors...])

function DimBasis(basis_vectors::Vararg{Pair{<:AbstractString,<:QuantityOrUnitlike}})
    to_quantity = any([isa(bv.second, Unitful.AbstractQuantity) for bv in basis_vectors])
    new_basis_vectors = [to_quantity ? 1 * bv.second : bv.second for bv in basis_vectors]
    new_basis_vector_names = [bv.first for bv in basis_vectors]
    return DimBasis(new_basis_vectors, new_basis_vector_names)
end

# Do not broadcast DimBasis
Base.Broadcast.broadcastable(basis::DimBasis) = Ref(basis)

# Helper Types
QuantityDimBasis = DimBasis{<:AbstractVector{<:Unitful.AbstractQuantity}}

NamedDimBase = DimBasis{T,<:AbstractVector{<:AbstractString}} where T

"""
    unique_dims(all_values...)

Return a vector of unique dimensions for `all_values`, a set of quantities, units or dimensions.
"""
function unique_dims(all_values::Vararg{Unitful.Dimensions})
    basis_dims = Vector{Type{<:Unitful.Dimension}}()
    for dims in all_values
        union!(basis_dims, typeof.(typeof(dims).parameters[1]))
    end
    return basis_dims
end

unique_dims(all_values::Vararg{QuantityOrUnits}) =
    unique_dims(dimension.(all_values)...)

"""
    dim_matrix(basis_dims, all_values...)

Return the dimensional matrix for a set of basis dimensions `basis_dims` and `all_values`, a set of quantities, units or dimensions.
"""
function dim_matrix(basis_dims::Array{Type{<:Unitful.Dimension}}, all_values::Vararg{Unitful.Dimensions})
    dim_mat = zeros(Rational, length(basis_dims), length(all_values))
    for (dimension_ind,dims) in enumerate(all_values)
        for dimension in typeof(dims).parameters[1]
            basis_dim_ind = findfirst(x -> isa(dimension, x), basis_dims)
            if isnothing(basis_dim_ind) 
                error("$(typeof(dimension)) is no basis dimension.")
            end
            dim_mat[basis_dim_ind,dimension_ind] = dimension.power
        end
    end
    return dim_mat
end

dim_matrix(basis_dims::Array{Type{<:Unitful.Dimension}}, all_values::Vararg{QuantityOrUnits}) =
dim_matrix(basis_dims, dimension.(all_values)...)

"""
    check_basis(dim_mat)

Use a dimensional matrix to check if a collection of dimensional vectors is a valid basis.
Throw errors if there are to few basis vectors or if the matrix does not have full rank.
"""
function check_basis(dim_mat)
    if size(dim_mat, 2) < size(dim_mat, 1)
        plr_sgl = (size(dim_mat, 2) == 1) ? "vector" : "vectors"
        error("Invalid basis: There are $(size(dim_mat, 1)) dimensions but only $(size(dim_mat, 2)) basis $(plr_sgl).")
    end

    mat_rank = LinearAlgebra.rank(dim_mat)
    if mat_rank < size(dim_mat, 2)
        if mat_rank == 1
            error("Invalid basis: There are $(size(dim_mat, 2)) basis vectors that are all linearly dependent.")
        else
            error("Invalid basis: There are $(size(dim_mat, 2)) basis vectors of which only $(mat_rank) are linearly independent.")
        end
    end
    nothing
end
