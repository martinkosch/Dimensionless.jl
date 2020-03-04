module Dimensionless

using Unitful, LinearAlgebra

export dim_basis, dim_matrix, dim_matrix_spy, dimensionless, dimensionful, change_basis

abstract type AbstractDimBasis end

Base.Broadcast.broadcastable(basis::AbstractDimBasis) = Ref(basis)

struct UnnamedDimBasis <: AbstractDimBasis
    basis_vectors
    basis_dims
    dim_mat

    function UnnamedDimBasis(basis_vectors::Vararg{Unitful.AbstractQuantity})
        basis_dims = unique_dims(basis_vectors...)
        dim_mat = dim_matrix(basis_dims, basis_vectors...)
        check_basis(dim_mat)
        return new(collect(basis_vectors), basis_dims, dim_mat)
    end
end

struct NamedDimBasis <: AbstractDimBasis
    basis_vectors
    basis_vector_names
    basis_dims
    dim_mat

    function NamedDimBasis(named_basis_vectors::Vararg{Pair{String,<:Unitful.AbstractQuantity}})
        basis_vector_names = [named_basis_vector.first for named_basis_vector in named_basis_vectors]
        basis_vectors = [named_basis_vector.second for named_basis_vector in named_basis_vectors]
        basis_dims = unique_dims(basis_vectors...)
        dim_mat = dim_matrix(basis_dims, basis_vectors...)
        check_basis(dim_mat)
        return new(basis_vectors, basis_vector_names, basis_dims, dim_mat)
    end
end

dim_basis(basis_vectors::Vararg{Unitful.AbstractQuantity}) = UnnamedDimBasis(basis_vectors...)

dim_basis(named_basis_vectors::Vararg{Pair{String,<:Unitful.AbstractQuantity}}) = NamedDimBasis(named_basis_vectors...)


function unique_dims(all_dims::Vararg{Unitful.Dimensions})
    basis_dims = Vector{Type{<:Unitful.Dimension}}()
    for dims in all_dims
        union!(basis_dims, typeof.(typeof(dims).parameters[1]))
    end
    return basis_dims
end

unique_dims(all_quantities::Vararg{Unitful.AbstractQuantity}) =
    unique_dims(dimension.(all_quantities)...)

function dim_matrix(basis_dims::Array{Type{<:Unitful.Dimension}}, all_dims::Vararg{Unitful.Dimensions})
    dim_mat = zeros(Rational, length(basis_dims), length(all_dims))
    for (dimension_ind,dims) in enumerate(all_dims)
        for dimension in typeof(dims).parameters[1]
            basis_dim_ind = findfirst(x -> isa(dimension, x), basis_dims)
            if basis_dim_ind == nothing
                error("$(typeof(dimension)) is no basis dimension!")
            end
            dim_mat[basis_dim_ind,dimension_ind] = dimension.power
        end
    end
    return dim_mat
end

dim_matrix(basis_dims::Array{Type{<:Unitful.Dimension}}, all_quantities::Vararg{Union{Unitful.AbstractQuantity, Unitful.Units}}) =
    dim_matrix(basis_dims, dimension.(all_quantities)...)

dim_matrix_spy(all_dims::Vararg{Unitful.Dimensions}) = dim_matrix(unique_dims(all_dims...), all_dims...)

dim_matrix_spy(all_quantities::Vararg{Union{Unitful.AbstractQuantity, Unitful.Units}}) = dim_matrix(unique_dims(dimension.(all_quantities)...), dimension.(all_quantities)...)

function check_basis(dim_mat)
    if size(dim_mat, 2) < size(dim_mat, 1)
        error("Invalid basis! There are $(size(dim_mat, 1)) dimensions but only $(size(dim_mat, 2)) basis vectors.")
    end
    mat_rank = LinearAlgebra.rank(dim_mat)
    if mat_rank < size(dim_mat, 2)
        if mat_rank == 1
            error("Invalid basis! There are $(size(dim_mat, 2)) basis vectors that are all linear dependent.")
        else
            error("Invalid basis! There are $(size(dim_mat, 2)) basis vectors of which only $(mat_rank) are linear independent.")
        end
    end
end

function dimensionless(quantity::Unitful.AbstractQuantity, basis::AbstractDimBasis)
    dim_mat = dim_matrix(basis.basis_dims, quantity)
    basis_multipliers = prod(basis.basis_vectors .^ (basis.dim_mat \ dim_mat))
    return ustrip(uconvert(Unitful.NoUnits, quantity / basis_multipliers))
end

function dimensionful(value::Number, unit::Unitful.Units, basis::AbstractDimBasis)
    dim_vec = dim_matrix(basis.basis_dims, dimension(unit))
    return uconvert(unit, value * prod(basis.basis_vectors .^ (basis.dim_mat \ dim_vec)))
end

function old_to_new_multiplier(dims::Unitful.Dimensions, old_basis::AbstractDimBasis, new_basis::AbstractDimBasis)
    old_dim_mat = dim_matrix(old_basis.basis_dims, dims)
    new_dim_mat = dim_matrix(new_basis.basis_dims, dims)
    multiplier = prod(new_basis.basis_vectors .^ (new_basis.dim_mat \ new_dim_mat)) / prod(old_basis.basis_vectors .^ (old_basis.dim_mat \ old_dim_mat))
    return ustrip(uconvert(Unitful.NoUnits, multiplier))
end

change_basis(quantity::Unitful.AbstractQuantity, old_basis::AbstractDimBasis, new_basis::AbstractDimBasis) = quantity * old_to_new_multiplier(dimension(quantity), old_basis, new_basis)

change_basis(unit::Unitful.Units, old_basis::AbstractDimBasis, new_basis::AbstractDimBasis) = old_to_new_multiplier(dimension(unit), old_basis, new_basis)

end # module
