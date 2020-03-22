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
    basis_facs = prod(basis.basis_vectors .^ (basis.dim_mat \ dim_mat))
    return ustrip(uconvert(Unitful.NoUnits, quantity / basis_facs))
end

function dimensionful(value::Number, unit::Unitful.Units, basis::QuantityDimBasis)
    dim_vec = dim_matrix(basis.basis_dims, dimension(unit))
    return uconvert(unit, value * prod(basis.basis_vectors .^ (basis.dim_mat \ dim_vec)))
end

function current_to_new_fac(dims::Unitful.Dimensions, basis::QuantityDimBasis, new_basis::QuantityDimBasis)
    dim_mat = dim_matrix(basis.basis_dims, dims)
    new_dim_mat = dim_matrix(new_basis.basis_dims, dims)
    basis_fac = prod(basis.basis_vectors .^ (basis.dim_mat \ dim_mat))
    new_basis_multipler = prod(new_basis.basis_vectors .^ (new_basis.dim_mat \ new_dim_mat))
    current_to_new_fac = new_basis_multipler / basis_fac
    return ustrip(uconvert(Unitful.NoUnits, current_to_new_fac))
end

change_basis(quantity::Unitful.AbstractQuantity, basis::QuantityDimBasis, new_basis::QuantityDimBasis) =
quantity * current_to_new_fac(dimension(quantity), basis, new_basis)

change_basis(unit::Unitful.Units, basis::QuantityDimBasis, new_basis::QuantityDimBasis) =
current_to_new_fac(dimension(unit), basis, new_basis)
