# Example
The following example is taken from the [Similitude Wikipedia article](https://en.wikipedia.org/wiki/Similitude_(model)):

>Consider a submarine modeled at 1/40th scale. The application operates in sea water at 0.5 °C, moving at 5 m/s. The model will be tested in fresh water at 20 °C. Find the power required for the submarine to operate at the stated speed.
>
>A free body diagram is constructed and the relevant relationships of force and velocity are formulated using techniques from continuum mechanics. The variables which describe the system are:
>
> Variable                   | Application   | Scaled model      |Units
>:---------------------------|:--------------|:------------------|:--------------
> L (Diameter of submarine)  | 1             | 1/40              | (m)
> V (Speed)                  | 5             | Calculate         | (m/s)
> ρ (Density)                | 1028          | 998               | (kg/m^3)
> μ (Dynamic viscosity)      | 1.88x10^−3    | 1.00x10^−3        | Pa·s (N s/m^2)
> F (Force)                  | Calculate     | To be measured    |	N (kg m/s^2)

Let's analyze this exemplary problem using Dimensionless.jl. We start by calculating the number of dimensions that characterize the problem:

```julia
julia> number_of_dimensions(u"m", u"m/s", u"kg/m^3", u"Pa*s", u"N")
3
```

According to [Buckinghams π theorem](https://en.wikipedia.org/wiki/Buckingham_%CF%80_theorem), the problem is governed by two dimensionless numbers (5 variables - 3 dimensions = 2 dimensional numbers):
```julia
julia> number_of_dimensionless(u"m", u"m/s", u"kg/m^3", u"Pa*s", u"N")
2
```

In order to find the required power of the submarine, we need to find a dimensional basis for the problem. In Dimensionless.jl, this is done by creating an instance of type `DimBasis`. As calculated before, there need to be three linear independent basis vectors to span the dimensional space. While creating such a basis using the `DimBasis` function, it is helpful to name the used variables to allow pretty printing of the results later on:
```julia
julia> basis = DimBasis("L"=>u"m", "V"=>u"m/s", "ρ"=>u"kg/m^3", "μ"=>u"Pa*s") # Invalid basis
ERROR: Invalid basis! There are 4 basis vectors of which only 3 are linear independent.

julia> basis = DimBasis("L"=>u"m", "ρ"=>u"kg/m^3", "μ"=>u"Pa*s") # Valid basis
DimBasis{...}
```

A valid basis! Together with the two remaining variables ``v`` and ``F``, the two dimensionless numbers can be determined:
```julia
julia> print_dimensionless("V"=>u"m/s", basis)
V L ρ μ^-1

julia> print_dimensionless("F"=>u"N", basis)
F ρ μ^-2
```
The first found number is the Reynolds number ``\mathit{Re}`` and the second number can be identified as the product of the drag coefficient ``c_d`` and ``\mathit{Re}^2``. All equations describing the problem can now be scaled using these numbers.

Let's construct two bases that describe the submarine model and the real submarine, respectively:
```julia
julia> applications_basis = DimBasis("L"=>1u"m", "ρ"=>1028u"kg/m^3", "μ"=>1.88e-3u"Pa*s")
DimBasis{...}

julia> model_basis = DimBasis("L"=>0.025u"m", "ρ"=>998u"kg/m^3", "μ"=>1e-3u"Pa*s")
DimBasis{...}
```

It is easy to transform values from one dimensional bases to the other in order to fill the missing values in the table:
```julia
julia> v_model = change_basis(5u"m/s", applications_basis, model_basis)
109.58 m s^-1

julia> change_basis(u"N", model_basis, applications_basis)
3.4313

julia> P_application_fac = v_model * change_basis(u"W", model_basis, applications_basis)
17.156 m s^-1
```

Similarity between the model and the application is achieved with a model velocity of 109.6 m/s, which is
21.9 times higher than the application velocity. The measured drag forces acting on the model under these circumstances need to be multiplied by 3.4 to determine the drag forces acting on the full scale submarine. The necessary propulsion power equates to ``P_{application} = F_{model} \, \cdot \, v_{model} \, \cdot \, \frac{P_{application}}{P_{model}} = F_{model} \, \cdot \, 17.2 \, \frac{m}{s}``.

As a last point, it is often helpful to be able to remove and restore the dimensions of variables in a given dimensional basis:
```julia
julia> dim_vars = (1u"cm/g", 1u"mm/s")
(1 cm g^-1, 1 mm s^-1)

julia> dimless_vars = dimensionless.(dim_vars, model_basis)
(6.2375, 24.95)

julia> dimensionful.(dimless_vars, (u"cm/g", u"mm/s"), model_basis)
(1.0 cm g^-1, 1.0 mm s^-1)
```
