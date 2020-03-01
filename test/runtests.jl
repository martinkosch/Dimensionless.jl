using Dimensionless, Unitful
using Test

@testset "UnnamedDimensionBasis tests" begin
    @test_throws ErrorException dimension_basis(9.81u"m/s^2", 1u"mm", 1u"s")
    @test_throws ErrorException dimension_basis(9.81u"m/s^2")

    basis = dimension_basis(9.81u"m/s^2", 6371u"km", 1420788u"kg", 1u"mol/g")
    @test_throws ErrorException dimensional_matrix(basis.basis_dims, 1u"A")

    @test basis.dim_mat == Rational.([1 1 0 0; -2 0 0 0; 0 0 1 -1; 0 0 0 1])
    for dim in [:Mass, :Length, :Time, :Amount]
        @test Unitful.Dimension{dim} in basis.basis_dims
    end
    for dim in [:Temperature, :Current, :Luminosity]
        @test !(Unitful.Dimension{dim} in basis.basis_dims)
    end
end

@testset "NamedDimensionBasis tests" begin
    named_basis = dimension_basis("a"=>1u"s", "b"=>2u"mm", "c"=>3u"bar")
    unnamed_basis = dimension_basis(1u"s", 2u"mm", 3u"bar")
    @test named_basis.dim_mat == unnamed_basis.dim_mat
    @test named_basis.basis_dims == unnamed_basis.basis_dims
    @test named_basis.basis_vectors == unnamed_basis.basis_vectors
    @test named_basis.basis_vector_names == ["a", "b", "c"]
end

@testset "Change of basis tests" begin
    change_basis = dimension_basis(1u"mA", 2u"mol", 3u"Pa*s", 4u"h", 5u"kg")
    dim_less = dimensionless(1u"A/m^2", change_basis)
    @test dimensionful(dim_less, u"A/m^2", change_basis) == 1u"A/m^2"
end
