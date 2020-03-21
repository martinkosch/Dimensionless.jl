function dimensionless_string(named_quantity::Pair{String,<:Union{Unitful.AbstractQuantity,Unitful.Unitlike}}, basis::DimBasis{T,AbstractVector{String}} where T, seperator::String=" ")
    vars = [basis.basis_vector_names[i] * (abs(powers[i]) != 1 ? Unitful.superscript(abs(powers[i])) : "") for i in eachindex(powers)]
    (positive_powers, negative_powers) = pos_neg_powers(named_quantity.second, basis)
    out = IOBuffer()
    write(out, named_quantity.first)
    (!isempty(positive_powers)) && write(out, " ")
    write(out, join(vars[positive_powers], seperator))
    write(out, " / ")
    (length(negative_powers) > 1) && write(out, "(")
    write(out, join(vars[negative_powers], seperator))
    (length(negative_powers) > 1) && write(out, ")")

    return String(take!(out))
end

function pos_neg_powers(quantity::Union{Unitful.AbstractQuantity,Unitful.Unitlike}, basis::DimBasis{T,AbstractVector{String}} where T)
    powers = basis.dim_mat \ dim_matrix(basis.basis_dims, named_quantity.second)
    positive_powers = findall(x -> x<0, powers)
    negative_powers = findall(x -> x>0, powers)
    return (positive_powers, negative_powers)
end
