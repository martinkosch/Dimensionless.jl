NamedDimVector = Pair{<:AbstractString,<:QuantityOrUnitlike}

# function dimensionless_string(named_quantity::NamedQuantity, basis::NamedDimBase)
#
#     out = IOBuffer()
#     write(out, named_quantity.first)
#     (!isempty(positive_powers)) && write(out, " ")
#     write(out, join(vars[positive_powers], seperator))
#     write(out, " / ")
#     (length(negative_powers) > 1) && write(out, "(")
#     write(out, join(vars[negative_powers], seperator))
#     (length(negative_powers) > 1) && write(out, ")")
#
#     return String(take!(out))
# end

function print_dimensionless(io::IO, named_dim_vector::NamedDimVector, basis::NamedDimBase)
    powers = -1 * (basis.dim_mat \ dim_matrix(basis.basis_dims, named_dim_vector.second))[:]
    perm = exps_perm(powers)

    showoperators = get(io, :showoperators, false)
    sep = showoperators ? "*" : " "
    showrep(io, named_dim_vector.first, 1//1)
    for i in perm
        if powers[i] != 0
            print(io, sep)
            showrep(io, basis.basis_vector_names[i], powers[i])
        end
    end
    nothing
end

print_dimensionless(named_dim_vector::NamedDimVector, basis::NamedDimBase) =
print_dimensionless(stdout, named_dim_vector, basis)

function showrep(io::IO, identifier::String, exp::Rational)
        print(io, identifier)
        if exp != 1
            print(io, Unitful.superscript(exp))
        end
    nothing
end

"""
    exps_perm(powers)
Calculate the permutation of exponents to show large positive values first followed by large negative values.
"""
function exps_perm(powers)
    perm = sortperm(powers; rev=true)
    first_negative = count(powers .> 0) + 1
    perm[first_negative:end] .= reverse(perm[first_negative:end])
    return perm
end
