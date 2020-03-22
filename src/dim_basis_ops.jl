"""
     number_of_dims(all_vars)
Return the number of unique dimensions for the given variables `all_vars.
"""
function number_of_dims(all_vars::Vararg{Unitful.Dimensions})
    return length(unique_dims(all_vars...))
end

number_of_dims(all_vars::Vararg{Pair{<:AbstractString,<:QuantityOrUnitlike}}) =
number_of_dims([var.second for var in all_vars]...)

number_of_dims(all_vars::Vararg{<:QuantityOrUnits}) =
number_of_dims(dimension.(all_vars)...)

"""
     number_of_dimensionless(all_vars)
Return the number of dimensionless numbers that can be constructed for a problem characterized by the variables `all_vars`.
"""
function number_of_dimensionless(all_vars::Vararg{<:QuantityOrUnitlike})
    return length(all_vars) - number_of_dims(all_vars...)
end

number_of_dimensionless(all_vars::Vararg{Pair{<:AbstractString,<:QuantityOrUnitlike}}) =
number_of_dimensionless([var.second for var in all_vars]...)

"""
     dimensionless(quantity, basis)
Make a `quantity` dimensionless using a dimensional `basis`.
"""
function dimensionless(quantity::Unitful.AbstractQuantity, basis::QuantityDimBasis)
    dim_mat = dim_matrix(basis.basis_dims, quantity)
    basis_facs = prod(basis.basis_vectors .^ (basis.dim_mat \ dim_mat))
    return ustrip(uconvert(Unitful.NoUnits, quantity / basis_facs))
end

"""
     dimensionful(value, unit, basis)
Restore the `unit`s of a dimensionless `value` using a dimensional `basis`.
"""
function dimensionful(value::Number, unit::Unitful.Units, basis::QuantityDimBasis)
    dim_vec = dim_matrix(basis.basis_dims, dimension(unit))
    return uconvert(unit, value * prod(basis.basis_vectors .^ (basis.dim_mat \ dim_vec)))
end

"""
     current_to_new_fac(dims, basis, new_basis)
Return the factor that is needed to transform specified dimensions `dims` from a current `basis` to a `new_basis`.
"""
function current_to_new_fac(dims::Unitful.Dimensions, basis::QuantityDimBasis, new_basis::QuantityDimBasis)
    dim_mat = dim_matrix(basis.basis_dims, dims)
    new_dim_mat = dim_matrix(new_basis.basis_dims, dims)
    basis_fac = prod(basis.basis_vectors .^ (basis.dim_mat \ dim_mat))
    new_basis_multipler = prod(new_basis.basis_vectors .^ (new_basis.dim_mat \ new_dim_mat))
    current_to_new_fac = new_basis_multipler / basis_fac
    return ustrip(uconvert(Unitful.NoUnits, current_to_new_fac))
end

"""
     change_basis(var, basis, new_basis)
Transform the specified quantity or unit `var` from a current `basis` to a `new_basis`.
"""
change_basis(var::Unitful.AbstractQuantity, basis::QuantityDimBasis, new_basis::QuantityDimBasis) =
var * current_to_new_fac(dimension(var), basis, new_basis)

change_basis(var::Unitful.Units, basis::QuantityDimBasis, new_basis::QuantityDimBasis) =
current_to_new_fac(dimension(var), basis, new_basis)
