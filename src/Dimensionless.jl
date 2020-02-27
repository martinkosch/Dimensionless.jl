module Dimensionless

using Unitful, LinearAlgebra

export BasisQuantities, dimensionless, dimensionful

struct BasisQuantities
    quantities::Tuple
    dimensions::Set
    dimensional_matrix::Array{Rational}
    function BasisQuantities(quantities::Vararg{Unitful.Quantity})
        dimensions = Set{Type{<:Unitful.Dimension}}()
        for quantity in quantities
            union!(dimensions, typeof.(typeof(dimension(quantity)).parameters[1]))
        end

        dimensional_matrix = zeros(Rational, length(dimensions), length(quantities))
        for (quantity_ind,quantity) in enumerate(quantities)
            for quantity_dim in typeof(typeof(quantity).parameters[2]).parameters[1]
                for (dim_ind,dim) in enumerate(dimensions)
                    if isa(quantity_dim, dim)
                        dimensional_matrix[dim_ind,quantity_ind] = quantity_dim.power
                    end
                end
            end
        end

        if LinearAlgebra.rank(dimensional_matrix) != length(quantities)
            error("Base dimensions are not linear independent!")
        end
        return new(quantities, dimensions, dimensional_matrix)
    end
end

function dimensionless(quantity::Unitful.Quantity, basis_quantities::BasisQuantities)

end

function dimensionful(quantity::Union{Unitful.NoUnits, Number}, unit::Unitful.Units, basis_quantities::BasisQuantities)

end

end # module
