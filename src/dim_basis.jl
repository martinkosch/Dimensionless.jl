QuantityOrUnits = Union{Unitful.AbstractQuantity,Unitful.Units}
QuantityOrUnitlike = Union{Unitful.AbstractQuantity, Unitful.Unitlike}

struct DimBasis{T,N}
    basis_vectors::T
    basis_vector_names::N
    basis_dims
    dim_mat

    function DimBasis(basis_vectors::T, basis_vector_names::N=nothing) where
        T<:AbstractVector{<:QuantityOrUnitlike} where
        N<:Union{Nothing,AbstractVector{<:AbstractString}}
        basis_dims = unique_dims(basis_vectors...)
        dim_mat = dim_matrix(basis_dims, basis_vectors...)
        check_basis(dim_mat)
        return new{T,N}(basis_vectors, basis_vector_names, basis_dims, dim_mat)
    end
end

Base.Broadcast.broadcastable(basis::DimBasis) = Ref(basis)

QuantityDimBasis = DimBasis{<:AbstractVector{<:Unitful.AbstractQuantity}}

NamedDimBase = DimBasis{T,<:AbstractVector{<:AbstractString}} where T

DimBasis(basis_vectors::Vararg{<:QuantityOrUnitlike}) =
DimBasis([basis_vectors...])

function DimBasis(named_basis_vectors::Vararg{Pair{<:AbstractString,<:QuantityOrUnitlike}})
    basis_vectors = [named_basis_vector.second for named_basis_vector in named_basis_vectors]
    basis_vector_names = [named_basis_vector.first for named_basis_vector in named_basis_vectors]
    return DimBasis(basis_vectors, basis_vector_names)
end

function unique_dims(all_dims::Vararg{Unitful.Dimensions})
    basis_dims = Vector{Type{<:Unitful.Dimension}}()
    for dims in all_dims
        union!(basis_dims, typeof.(typeof(dims).parameters[1]))
    end
    return basis_dims
end

unique_dims(all_quantities::Vararg{QuantityOrUnits}) =
    unique_dims(dimension.(all_quantities)...)

function dim_matrix(basis_dims::Array{Type{<:Unitful.Dimension}}, all_dims::Vararg{Unitful.Dimensions})
    dim_mat = zeros(Rational, length(basis_dims), length(all_dims))
    for (dimension_ind,dims) in enumerate(all_dims)
        for dimension in typeof(dims).parameters[1]
            basis_dim_ind = findfirst(x -> isa(dimension, x), basis_dims)
            if basis_dim_ind == nothing
                error("$(typeof(dimension)) is no basis dimension!")
            end
            dim_mat[basis_dim_ind,dimension_ind] = dimension.power
        end
    end
    return dim_mat
end

dim_matrix(basis_dims::Array{Type{<:Unitful.Dimension}}, all_quantities::Vararg{<:QuantityOrUnits}) =
    dim_matrix(basis_dims, dimension.(all_quantities)...)

function check_basis(dim_mat)
    if size(dim_mat, 2) < size(dim_mat, 1)
        plr_sgl = (size(dim_mat, 2) == 1) ? "vector" : "vectors"
        error("Invalid basis! There are $(size(dim_mat, 1)) dimensions but only $(size(dim_mat, 2)) basis $(plr_sgl).")
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
