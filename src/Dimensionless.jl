module Dimensionless

using Unitful, LinearAlgebra

export DimensionBasis, dimensional_matrix, dimensionless, dimensionful

struct DimensionBasis
    basis_vectors
    dim_set
    dim_mat

    function DimensionBasis(basis_vectors::Vararg{Unitful.AbstractQuantity})
        dim_set = unique_dimensions(basis_vectors...)
        dim_mat = dimensional_matrix(dim_set, basis_vectors...)

        if size(dim_mat, 2) < size(dim_mat, 1)
            error("Too few basis vectors!")
        elseif LinearAlgebra.rank(dim_mat) != length(basis_vectors)
            error("Basis vectors are not linear independent!")
        end
        return new(basis_vectors, dim_set, dim_mat)
    end
end

function unique_dimensions(quantities::Vararg{Unitful.AbstractQuantity})
    dim_set = Set{Type{<:Unitful.Dimension}}()
    for quantity in quantities
        union!(dim_set, typeof.(typeof(dimension(quantity)).parameters[1]))
    end
    return dim_set
end

function dimensional_matrix(dim_set::Set{Type{<:Unitful.Dimension}}, quantities::Vararg{Unitful.AbstractQuantity})
    dim_mat = zeros(Rational, length(dim_set), length(quantities))
    for (quantity_ind,quantity) in enumerate(quantities)
        quantity_dims = typeof(typeof(quantity).parameters[2]).parameters[1]
        for quantity_dim in quantity_dims
            for (dim_ind,dim) in enumerate(dim_set)
                if isa(quantity_dim, dim)
                    dim_mat[dim_ind,quantity_ind] = quantity_dim.power
                end
            end
        end
    end
    return dim_mat
end

function dimensionless(basis::DimensionBasis, quantity::Unitful.AbstractQuantity)
    dim_mat = dimensional_matrix(basis.dim_set, quantity...)
    basis_multipliers = prod(basis.basis_vectors.^(basis.dim_mat\dim_mat))
    return ustrip(uconvert(Unitful.NoUnits, quantity ./ basis_multipliers))
end


end # module
