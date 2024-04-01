# Example
The following example is taken from the [Similitude Wikipedia article](https://en.wikipedia.org/wiki/Similitude_(model)):

>Consider a submarine modeled at 1/40th scale. The application operates in sea water at 0.5 °C, moving at 5 m/s. The model will be tested in fresh water at 20 °C. Find the power required for the submarine to operate at the stated speed.
>
>A free body diagram is constructed and the relevant relationships of force and velocity are formulated using techniques from continuum mechanics. The variables which describe the system are:
>
> Variable                   | Application   | Scaled model      |Units
>:---------------------------|:--------------|:------------------|:--------------
> ``L`` (Diameter of submarine)  | 1             | 1/40              | (m)
> ``V`` (Speed)                  | 5             | Calculate         | (m/s)
> ``ρ`` (Density)                | 1028          | 998               | (kg/m^3)
> ``μ`` (Dynamic viscosity)      | 1.88x10^−3    | 1.00x10^−3        | Pa·s (N s/m^2)
> ``F`` (Force)                  | Calculate     | To be measured    |	N (kg m/s^2)


Let's analyze this exemplary problem using Dimensionless.jl. The idea of similitude is that many physical problems, despite appearing different in scale, geometry, or conditions, can exhibit similar behavior when analyzed through the lens of dimensionless parameters and scaling laws. 

## Analyzing the number of independent variables 
According to [Buckingham's Π theorem](https://en.wikipedia.org/wiki/Buckingham_%CF%80_theorem), a problem defined by ``n`` variables with corresponding units can be transformed into an equation with ``n − m`` dimensionless quantities, where ``m`` represents the number of independent base units used. 

Translated to the exemplary submarine similitude problem stated above with ``n = 5`` variables, Dimensionless.jl can be used to calculate that there are ``m = 3`` independent base units in the present case:
```julia
julia> m = number_of_dimensions(u"m", u"m/s", u"kg/m^3", u"Pa*s", u"N")
3
```

There are ``n - m = 5 - 3 = 2`` dimensionless variables, or Pi terms, that should be selected to effectively capture the essence of the physical phenomena under investigation: 
```julia
julia> number_of_dimensionless(u"m", u"m/s", u"kg/m^3", u"Pa*s", u"N")
2
```

## Finding a dimensional basis
The first step to set up the similitude problem is to find a suitable set of ``m`` independent variables that serve as a dimensional basis. These variables must span the dimensional space of the studied problem. 
In Dimensionless.jl, constructing such a dimensional basis is achieved by creating an instance of type `DimBasis`. While creating such a basis using the `DimBasis` function, it is helpful to name the used variables to allow pretty printing of the results later on:
```julia
julia> basis = DimBasis("L"=>u"m", "ρ"=>u"kg/m^3", "μ"=>u"Pa*s") # Valid basis
DimBasis{...}
```
In this case, a suitable basis could be found by selecting the variables ``L``, ``ρ``, and ``μ`` as basis vectors. In general, this set of variables must be chosen carefully. For exmaple, selecting more variables than needed leads to an error:
```julia
julia> basis = DimBasis("L"=>u"m", "V"=>u"m/s", "ρ"=>u"kg/m^3", "μ"=>u"Pa*s") # Invalid basis
ERROR: Invalid basis! There are 4 basis vectors of which only 3 are linear independent.
```


Together with the two remaining variables ``v`` and ``F``, the two dimensionless numbers can be determined:
```julia
julia> print_dimensionless("V"=>u"m/s", basis)
V L ρ μ^-1

julia> print_dimensionless("F"=>u"N", basis)
F ρ μ^-2
```
The first found number is the Reynolds number ``\mathit{Re}`` and the second number can be identified as the product of the drag coefficient ``c_\mathrm{d}`` and ``\mathit{Re}^2``. All equations describing the problem can now be scaled using these numbers.

## Calculating the missing entries in the table
Let's construct two bases that describe the submarine model and the submarine in the real application, respectively:
```julia
julia> applications_basis = DimBasis("L"=>1u"m", "ρ"=>1028u"kg/m^3", "μ"=>1.88e-3u"Pa*s")
DimBasis{...}

julia> model_basis = DimBasis("L"=>0.025u"m", "ρ"=>998u"kg/m^3", "μ"=>1e-3u"Pa*s")
DimBasis{...}
```

It is easy to transform values from one dimensional bases to the other in order to fill the missing values in the table. Starting with the submarine's velocity of 5 m/s in the application case, the scaled velocity of the model should be selected to: 
```julia
julia> v_model = change_basis(5u"m/s", applications_basis, model_basis)
109.58 m s^-1
```
This model velocity of 109.6 m/s is 21.9 times higher than the velocity of the submarine in the real application.

The factor to convert a given force in Newton from the submarine model to the real application can also be calculated: 
```
julia> change_basis(u"N", model_basis, applications_basis)
3.4313
```

Finally, the power of the real submarine can be derived with ``P_\mathrm{application} = F_\mathrm{model} \, \cdot \, v_\mathrm{model} \, \cdot \, \frac{P_\mathrm{application}}{P_\mathrm{model}}``. The following conversion factor for both power values can be calculated using the package:
```
julia> P_application_fac = v_model * change_basis(u"W", model_basis, applications_basis)
17.156 m s^-1
```
The final power of the real submarine can thus be calculated with the formula ``P_\mathrm{application} = F_\mathrm{model} \, \cdot \, 17.2 \, \frac{\mathrm{m}}{\mathrm{s}}``.

## Removing and restoring dimensions
Arguably the two most important functions of this package are `dimensionless()` and `dimensionful()`. For a given dimensional basis, these two functions can be used to convert variables from a formulation with dimensions to a corresponding dimensionless formulation and reverse: 
```julia
julia> dim_vars = (1u"cm/g", 1u"mm/s")
(1 cm g^-1, 1 mm s^-1)

julia> dimless_vars = dimensionless.(dim_vars, model_basis)
(6.2375, 24.95)

julia> dimensionful.(dimless_vars, (u"cm/g", u"mm/s"), model_basis)
(1.0 cm g^-1, 1.0 mm s^-1)
```
