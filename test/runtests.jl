using Dimensionless, Unitful
using Test

@testset "DimBasis construction" begin
    # Check invalid basis exceptions
    @test_throws ErrorException DimBasis(9.81u"m/s^2", 1u"mm", 1u"s")
    @test_throws ErrorException DimBasis(9.81u"m/s^2")

    # Check dim_basis exception in case of unvalid dimension
    basis = DimBasis(9.81u"m/s^2", 6371u"km", 1420788u"kg", 1u"mol/g")
    @test_throws ErrorException dim_matrix(basis.basis_dims, 1u"A")

    @test basis.dim_mat == Rational.([1 1 0 0; -2 0 0 0; 0 0 1 -1; 0 0 0 1])
    for dim in [:Mass, :Length, :Time, :Amount]
        @test Unitful.Dimension{dim} in basis.basis_dims
    end
    for dim in [:Temperature, :Current, :Luminosity]
        @test !(Unitful.Dimension{dim} in basis.basis_dims)
    end

    # Check all constructor combinations for DimBasis
    basis_vector_names_template = ["α", "β", "γ"]
    basis_vectors_template = [1u"s^2", 2u"mm^-1", 3u"bar^3"]
    for (broadcast_fcn, isnamed) in Iterators.product((identity, unit, dimension), (false, true))
        # Construct bases
        input_vectors = broadcast(broadcast_fcn, basis_vectors_template)
        basis = isnamed ?
        DimBasis([Pair(basis_vector_names_template[i], input_vectors[i]) for i in eachindex(basis_vectors_template)]...) :
        DimBasis(input_vectors...)

        # Check validity
        isnamed && @test basis.basis_vector_names == basis_vector_names_template
        @test basis.dim_mat == Rational.([2 0 -6; 0 -1 -3; 0 0 3])
        @test basis.basis_dims == [Unitful.Dimension{:Time}, Unitful.Dimension{:Length}, Unitful.Dimension{:Mass}]
        @test basis.basis_vectors == input_vectors
    end
end

@testset "Change of basis" begin
    change_basis = DimBasis(1u"mA", 2u"mol", 3u"Pa*s", 4u"h", 5u"kg")
    dim_less = dimensionless(1u"A/m^2", change_basis)
    @test dimensionful(dim_less, u"A/m^2", change_basis) == 1u"A/m^2"
end
