using Dimensionless, Unitful
using Test

@testset "DimBasis construction" begin
    # Check invalid basis construction exceptions
    @test_throws ErrorException DimBasis(9.81u"m/s^2", 1u"mm", 1u"s")
    @test_throws ErrorException DimBasis(9.81u"m", 9.81u"m")
    @test_throws ErrorException DimBasis(9.81u"m/s^2")

    # Check exception in dim_matrix() in case of unvalid dimension
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
    # Test change of basis and broadcast
    dimless_dimful_basis = DimBasis(1u"mA", 2u"K", 3u"Pa*s", 4u"h", 5u"kg")
    vars = [0u"A/m^2", 1u"nm", 2u"K"]
    dim_less = dimensionless.(vars, dimless_dimful_basis)
    @test dimensionful.(dim_less, [u"A/m^2", u"nm", u"K"], dimless_dimful_basis) == vars

    # Test counting of dimensions and dimensionless variables for named basis
    quantities_named = ["a"=>100.0u"m", "b"=>100u"kg", "c"=>0u"s"]
    units_named = [Pair(var.first, unit(var.second)) for var in quantities_named]
    dims_named = [Pair(var.first, dimension(var.second)) for var in quantities_named]
    @test number_of_dims(quantities_named...) == 3
    @test number_of_dims(units_named...) == 3
    @test number_of_dims(dims_named...) == 3
    @test number_of_dimensionless(quantities_named...) == 0
    @test number_of_dimensionless(units_named...) == 0
    @test number_of_dimensionless(dims_named...) == 0

    # Test counting of dimensions and dimensionless variables for unnamed basis
    quantities_unnamed = [var.second for var in quantities_named]
    units_unnamed = [var.second for var in units_named]
    dims_unnamed = [var.second for var in dims_named]
    @test number_of_dims(quantities_unnamed...) == 3
    @test number_of_dims(units_unnamed...) == 3
    @test number_of_dims(dims_unnamed...) == 3
    @test number_of_dimensionless(quantities_unnamed...) == 0
    @test number_of_dimensionless(units_unnamed...) == 0
    @test number_of_dimensionless(dims_unnamed...) == 0

    # change_basis(quantity/unit, old_b, new_b) for named and unnamed quantity bases
    basis_named = DimBasis("x"=>1u"m", "y"=>1u"kg", "z"=>1u"s")
    new_basis_named = DimBasis("x_new"=>2u"m", "y_new"=>2u"kg", "z_new"=>2u"s")
    @test Dimensionless.current_to_new_fac(Unitful.𝐋*Unitful.𝐓*Unitful.𝐌, basis_named, new_basis_named) == 2^3
    @test change_basis(100u"cm*kg*s", basis_named, new_basis_named) == 8u"m*kg*s"
    @test change_basis(u"m*kg*s", basis_named, new_basis_named) == 8

    basis_unnamed = DimBasis(1u"m", 1u"kg", 1u"s")
    new_basis_unnamed = DimBasis(2u"m", 2u"kg", 2u"s")
    @test Dimensionless.current_to_new_fac(Unitful.𝐋*Unitful.𝐓*Unitful.𝐌, basis_unnamed, new_basis_unnamed) == 2^3
    @test change_basis(100u"cm*kg*s", basis_unnamed, new_basis_unnamed) == 8u"m*kg*s"
    @test change_basis(u"m*kg*s", basis_unnamed, new_basis_unnamed) == 8
end

@testset "Utils" begin
    var = "L"=>0.2u"m"
    quantities_named = ["ρ"=>1400.0u"kg/m^3", "v"=>0.5u"m/s", "η"=>(10^4)u"Pa*s"]
    units_named = [Pair(var.first, unit(var.second)) for var in quantities_named]
    dims_named = [Pair(var.first, dimension(var.second)) for var in quantities_named]
    res_desired = "L ρ v η^-1"
    res_buf = IOBuffer()

    # Test dimensionless string for named basis
    print_dimensionless(res_buf, var, DimBasis(quantities_named...))
    @test String(take!(res_buf)) == res_desired
    print_dimensionless(res_buf, var, DimBasis(units_named...))
    @test String(take!(res_buf)) == res_desired
    print_dimensionless(res_buf, var, DimBasis(dims_named...))
    @test String(take!(res_buf)) == res_desired
end
