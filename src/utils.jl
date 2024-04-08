NamedDimVector = Pair{<:AbstractString,<:QuantityOrUnitlike}

"""
    print_dimless([io, ]named_dim_vector, basis)

Print the dimensionless number that can be constructed using `named_dim_vector` in the specified dimensional `basis`.
"""
function print_dimless(io::IO, named_dim_vector::NamedDimVector, basis::NamedDimBase)
    powers = -1 * (basis.dim_mat\dim_matrix(basis.basis_dims, named_dim_vector.second))[:]
    perm = exps_perm(powers)

    showoperators = get(io, :showoperators, false)
    sep = showoperators ? "*" : " "
    showrep(io, named_dim_vector.first, 1 // 1)
    for i in perm
        if powers[i] != 0
            print(io, sep)
            showrep(io, basis.basis_vector_names[i], powers[i])
        end
    end
    nothing
end

print_dimless(named_dim_vector::NamedDimVector, basis::NamedDimBase) =
    print_dimless(stdout, named_dim_vector, basis)

"""
    showrep(io, identifier, exponent)

Print the `identifier` followed by the corresponding `exponent` to the specified stream `io`.
"""
function showrep(io::IO, identifier::String, exp::Rational)
    print(io, identifier)
    if exp != 1
        exp_str = exp.den == 1 ? "^" * string(exp.num) : "^" * replace(string(exp), "//" => "/")
        print(io, exp_str)
    end
    nothing
end

"""
    exps_perm(powers)
Return the permutation of exponents so that large positive values are showed first followed by large negative values.
"""
function exps_perm(powers)
    perm = sortperm(powers; rev=true)
    first_negative = count(powers .> 0) + 1
    perm[first_negative:end] .= reverse(perm[first_negative:end])
    return perm
end
