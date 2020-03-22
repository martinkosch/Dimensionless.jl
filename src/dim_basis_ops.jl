function number_of_dims(all_vars::Vararg{Unitful.Dimensions})
    return length(unique_dims(all_vars...))
end

number_of_dims(all_vars::Vararg{Pair{<:AbstractString,<:QuantityOrUnitlike}}) =
number_of_dims([var.second for var in all_vars]...)

number_of_dims(all_vars::Vararg{<:QuantityOrUnits}) =
number_of_dims(dimension.(all_vars)...)

function number_of_dimensionless(all_vars::Vararg{<:QuantityOrUnitlike})
    return length(all_vars) - number_of_dims(all_vars...)
end

number_of_dimensionless(all_vars::Vararg{Pair{<:AbstractString,<:QuantityOrUnitlike}}) =
number_of_dimensionless([var.second for var in all_vars]...)

function dimensionless(quantity::Unitful.AbstractQuantity, basis::QuantityDimBasis)
    dim_mat = dim_matrix(basis.basis_dims, quantity)
    basis_multipliers = prod(basis.basis_vectors .^ (basis.dim_mat \ dim_mat))
    return ustrip(uconvert(Unitful.NoUnits, quantity / basis_multipliers))
end

function dimensionful(value::Number, unit::Unitful.Units, basis::QuantityDimBasis)
    dim_vec = dim_matrix(basis.basis_dims, dimension(unit))
    return uconvert(unit, value * prod(basis.basis_vectors .^ (basis.dim_mat \ dim_vec)))
end

function old_to_new_multiplier(dims::Unitful.Dimensions, old_basis::QuantityDimBasis, new_basis::QuantityDimBasis)
    old_dim_mat = dim_matrix(old_basis.basis_dims, dims)
    new_dim_mat = dim_matrix(new_basis.basis_dims, dims)
    new_multipler = prod(new_basis.basis_vectors .^ (new_basis.dim_mat \ new_dim_mat))
    old_multipler = prod(old_basis.basis_vectors .^ (old_basis.dim_mat \ old_dim_mat))
    multiplier = new_multipler / old_multipler
    return ustrip(uconvert(Unitful.NoUnits, multiplier))
end

change_basis(quantity::Unitful.AbstractQuantity, old_basis::QuantityDimBasis, new_basis::QuantityDimBasis) =
quantity * old_to_new_multiplier(dimension(quantity), old_basis, new_basis)

change_basis(unit::Unitful.Units, old_basis::QuantityDimBasis, new_basis::QuantityDimBasis) =
old_to_new_multiplier(dimension(unit), old_basis, new_basis)
