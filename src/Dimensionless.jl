module Dimensionless

using Unitful, LinearAlgebra

export dimension_basis, dimensional_matrix, dimensional_matrix_spy, dimensionless, dimensionful, old_to_new_multiplier, change_basis

abstract type AbstractDimensionBasis end

Base.Broadcast.broadcastable(basis::AbstractDimensionBasis) = Ref(basis)

struct UnnamedDimensionBasis <: AbstractDimensionBasis
    basis_vectors
    basis_dims
    dim_mat

    function UnnamedDimensionBasis(basis_vectors::Vararg{Unitful.AbstractQuantity})
        basis_dims = unique_dimensions(basis_vectors...)
        dim_mat = dimensional_matrix(basis_dims, basis_vectors...)
        check_basis(dim_mat)
        return new(collect(basis_vectors), basis_dims, dim_mat)
    end
end

struct NamedDimensionBasis <: AbstractDimensionBasis
    basis_vectors
    basis_vector_names
    basis_dims
    dim_mat

    function NamedDimensionBasis(named_basis_vectors::Vararg{Pair{String,<:Unitful.AbstractQuantity}})
        basis_vector_names = [named_basis_vector.first for named_basis_vector in named_basis_vectors]
        basis_vectors = [named_basis_vector.second for named_basis_vector in named_basis_vectors]
        basis_dims = unique_dimensions(basis_vectors...)
        dim_mat = dimensional_matrix(basis_dims, basis_vectors...)
        check_basis(dim_mat)
        return new(basis_vectors, basis_vector_names, basis_dims, dim_mat)
    end
end

dimension_basis(basis_vectors::Vararg{Unitful.AbstractQuantity}) = UnnamedDimensionBasis(basis_vectors...)

dimension_basis(named_basis_vectors::Vararg{Pair{String,<:Unitful.AbstractQuantity}}) = NamedDimensionBasis(named_basis_vectors...)


function unique_dimensions(all_dimensions::Vararg{Unitful.Dimensions})
    basis_dims = Vector{Type{<:Unitful.Dimension}}()
    for dimensions in all_dimensions
        union!(basis_dims, typeof.(typeof(dimensions).parameters[1]))
    end
    return basis_dims
end

unique_dimensions(all_quantities::Vararg{Unitful.AbstractQuantity}) =
    unique_dimensions(dimension.(all_quantities)...)

function dimensional_matrix(basis_dims::Array{Type{<:Unitful.Dimension}}, all_dimensions::Vararg{Unitful.Dimensions})
    dim_mat = zeros(Rational, length(basis_dims), length(all_dimensions))
    for (dimension_ind,dimensions) in enumerate(all_dimensions)
        for dimension in typeof(dimensions).parameters[1]
            basis_dim_ind = findfirst(x -> isa(dimension, x), basis_dims)
            if basis_dim_ind == nothing
                error("$(typeof(dimension)) is no basis dimension!")
            end
            dim_mat[basis_dim_ind,dimension_ind] = dimension.power
        end
    end
    return dim_mat
end

dimensional_matrix(basis_dims::Array{Type{<:Unitful.Dimension}}, all_quantities::Vararg{Union{Unitful.AbstractQuantity, Unitful.Units}}) =
    dimensional_matrix(basis_dims, dimension.(all_quantities)...)

dimensional_matrix_spy(all_dimensions::Vararg{Unitful.Dimensions}) = dimensional_matrix(unique_dimensions(all_dimensions...), all_dimensions...)

dimensional_matrix_spy(all_quantities::Vararg{Union{Unitful.AbstractQuantity, Unitful.Units}}) = dimensional_matrix(unique_dimensions(dimension.(all_quantities)...), dimension.(all_quantities)...)

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

function dimensionless(quantity::Unitful.AbstractQuantity, basis::AbstractDimensionBasis)
    dim_mat = dimensional_matrix(basis.basis_dims, quantity)
    basis_multipliers = prod(basis.basis_vectors .^ (basis.dim_mat \ dim_mat))
    return ustrip(uconvert(Unitful.NoUnits, quantity / basis_multipliers))
end

function dimensionful(value::Number, unit::Unitful.Units, basis::AbstractDimensionBasis)
    dim_vec = dimensional_matrix(basis.basis_dims, dimension(unit))
    return uconvert(unit, value * prod(basis.basis_vectors .^ (basis.dim_mat \ dim_vec)))
end

function old_to_new_multiplier(dimensions::Unitful.Dimensions, old_basis::AbstractDimensionBasis, new_basis::AbstractDimensionBasis)
    old_dim_mat = dimensional_matrix(old_basis.basis_dims, dimensions)
    new_dim_mat = dimensional_matrix(new_basis.basis_dims, dimensions)
    multiplier = prod(new_basis.basis_vectors .^ (new_basis.dim_mat \ new_dim_mat)) / prod(old_basis.basis_vectors .^ (old_basis.dim_mat \ old_dim_mat))
    return ustrip(uconvert(Unitful.NoUnits, multiplier))
end

old_to_new_multiplier(quantity::Union{Unitful.AbstractQuantity, Unitful.Units}, old_basis::AbstractDimensionBasis, new_basis::AbstractDimensionBasis) = old_to_new_multiplier(dimension(quantity), old_basis, new_basis)

change_basis(quantity::Unitful.AbstractQuantity, old_basis::AbstractDimensionBasis, new_basis::AbstractDimensionBasis) = quantity * old_to_new_multiplier(quantity, old_basis, new_basis)

end # module
