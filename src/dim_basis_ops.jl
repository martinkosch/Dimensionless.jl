function count_dims(all_vars::Vararg{Unitful.Dimensions})
    return length(unique_dims(all_vars...))
end

count_dims(all_vars::Vararg{Pair{String,<:Union{Unitful.Dimensions, Unitful.AbstractQuantity, Unitful.Units}}}) =
count_dims([var.second for var in all_vars]...)

count_dims(all_vars::Vararg{Union{Unitful.AbstractQuantity, Unitful.Units}}) =
count_dims(dimension.(all_vars)...)

function count_dimensionless(all_vars::Vararg{<:Union{Unitful.AbstractQuantity,Unitful.Unitlike}})
    return length(all_vars) - count_dims(all_vars...)
end

count_dimensionless(all_vars::Vararg{Pair{String,<:Union{Unitful.Dimensions, Unitful.AbstractQuantity, Unitful.Units}}}) =
count_dimensionless([var.second for var in all_vars]...)

function dimensionless(quantity::Unitful.AbstractQuantity, basis::T where
     T<:DimBasis{<:AbstractVector{<:Unitful.AbstractQuantity}})
    dim_mat = dim_matrix(basis.basis_dims, quantity)
    basis_multipliers = prod(basis.basis_vectors .^ (basis.dim_mat \ dim_mat))
    return ustrip(uconvert(Unitful.NoUnits, quantity / basis_multipliers))
end

function dimensionful(value::Number, unit::Unitful.Units, basis::T where
     T<:DimBasis{<:AbstractVector{<:Unitful.AbstractQuantity}})
    dim_vec = dim_matrix(basis.basis_dims, dimension(unit))
    return uconvert(unit, value * prod(basis.basis_vectors .^ (basis.dim_mat \ dim_vec)))
end

function old_to_new_multiplier(dims::Unitful.Dimensions, old_basis::T, new_basis::T) where
     T<:DimBasis{<:AbstractVector{<:Unitful.AbstractQuantity}}
    old_dim_mat = dim_matrix(old_basis.basis_dims, dims)
    new_dim_mat = dim_matrix(new_basis.basis_dims, dims)
    new = prod(new_basis.basis_vectors .^ (new_basis.dim_mat \ new_dim_mat))
    old = prod(old_basis.basis_vectors .^ (old_basis.dim_mat \ old_dim_mat))
    multiplier = new / old
    return ustrip(uconvert(Unitful.NoUnits, multiplier))
end

change_basis(quantity::Unitful.AbstractQuantity, old_basis::T, new_basis::T) where
 T<:DimBasis{<:AbstractVector{<:Unitful.AbstractQuantity}} =
quantity * old_to_new_multiplier(dimension(quantity), old_basis, new_basis)

change_basis(unit::Unitful.Units, old_basis::T, new_basis::T) where
 T<:DimBasis{<:AbstractVector{<:Unitful.AbstractQuantity}} =
old_to_new_multiplier(dimension(unit), old_basis, new_basis)
