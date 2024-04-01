# Example
The following example is taken from the [Similitude Wikipedia article](https://en.wikipedia.org/wiki/Similitude_(model)):

>Consider a submarine modeled at 1/40th scale. The application operates in sea water at 0.5 °C, moving at 5 m/s. The model will be tested in fresh water at 20 °C. 
>
> The variables which describe the system are:
>
> Variable                       | Application          | Scaled model          |Units
>:-------------------------------|:---------------------|:----------------------|:-------------
> ``L`` (Diameter of submarine)  | 1.0                  | 1/40                  | (m)
> ``V`` (Speed)                  | 5.0                  | **To be calculated**  | (m/s)
> ``ρ`` (Density)                | 1028.0               | 998.0                 | (kg/m^3)
> ``μ`` (Dynamic viscosity)      | 1.88e−3              | 1.0e−3                | Pa·s (N s/m^2)
> ``F`` (Force)                  | **To be calculated** | **To be measured**    | N (kg m/s^2)
>
> **Find the power required for the submarine to operate at the stated speed.**

Let's analyze this exemplary problem using Dimensionless.jl. The idea of similitude is that physical problems, despite appearing different in scale, geometry, or conditions, can exhibit similar behavior when analyzed using a dimensionless formulation of the underlying equations.   

For the example at hand, the mechanical power `P` of the submarine can generally be expressed as the product of the propulsion force `F` and the velocity `V`. 
This relationship applies regardless of the specific problem parameters and the units used. However, the numerical values of the underlying variables change depending on the scaling of the problem and the employed units. For example, the numerical values and the used units in this equation formulated for the submarine in the real application (``P_\mathrm{app} = F_\mathrm{app} \cdot V_\mathrm{app}``) will often be different from the values and units for the scaled submarine model (``P_\mathrm{model} = F_\mathrm{model} \cdot V_\mathrm{model}``). Similitude facilitates to systematically transform all variables in this physically meaningful equation to a dimensionless formulation without the constraints of specific units or scales: ``\tilde{P} = \tilde{F} \cdot \tilde{V}``. 

In the following, it will be explained how such a transformation can be achieved and how it can be used to transfer knowledge between differently scaled versions of the same physical problem.   

## Analyzing the number of independent variables 
According to [Buckingham's Π theorem](https://en.wikipedia.org/wiki/Buckingham_%CF%80_theorem), a problem defined by `n` variables with corresponding units can be transformed into an equation with `n − m` dimensionless quantities, where `m` represents the number of independent base units used. 

Translated to the exemplary submarine similitude problem stated above with `n = 5` variables, Dimensionless.jl can be used to calculate that there are `m = 3` independent base units in the present case:
```julia
julia> m = number_of_dimensions(u"m", u"m/s", u"kg/m^3", u"Pa*s", u"N")
3
```

Correspondingly, there are `n - m = 5 - 3 = 2` dimensionless variables, or Pi terms, that should be selected to effectively capture the essence of the physical phenomena under investigation: 
```julia
julia> number_of_dimensionless(u"m", u"m/s", u"kg/m^3", u"Pa*s", u"N")
2
```

## Finding a dimensional basis
The first step to set up the similitude problem is to find a suitable set of `m` independent variables that serve as a dimensional basis. These variables must span the dimensional space of the studied problem. 
In Dimensionless.jl, constructing such a dimensional basis is achieved by creating an instance of type `DimBasis`. While creating such a basis using the `DimBasis` function, it is helpful to name the used variables to allow pretty printing of the results later on:
```julia
julia> basis = DimBasis("L"=>u"m", "ρ"=>u"kg/m^3", "μ"=>u"Pa*s") 
DimBasis{...}
```
In this case, a suitable basis could be found by selecting the `m = 3` variables `L`, `ρ`, and `μ` as basis variables. Using such a basis, every variable with a unit that lies inside the space spanned by this basis can be expressed as a linear combination of the basis variables. At the same time, each dimensioned variable can be transformed into a dimensionless representation using these basic variables. 

For example, any velocity `V` can be made dimensionless by multiplying it with a factor `L ρ / μ` as can be determined using Dimensionless.jl:  
```julia
julia> print_dimensionless("V"=>u"m/s", basis)
V L ρ μ^-1
```

The same is true for every force `F`:
```julia
julia> print_dimensionless("F"=>u"N", basis)
F ρ μ^-2
```
Note that the first found number is the Reynolds number ``\mathit{Re}`` and the second number can be identified as the product of the drag coefficient ``c_\mathrm{d}`` and ``\mathit{Re}^2``. In Dimensionless.jl, these numbers are mostly used under the hood to transform between formulations with and without dimensions.  

## Calculating the missing entries in the table
Let's use the values in the data table to construct two bases describing the submarine model and the submarine in the real application, respectively:
```julia
julia> app_basis = DimBasis("L"=>1u"m", "ρ"=>1028u"kg/m^3", "μ"=>1.88e-3u"Pa*s")
DimBasis{...}

julia> model_basis = DimBasis("L"=>0.025u"m", "ρ"=>998u"kg/m^3", "μ"=>1e-3u"Pa*s")
DimBasis{...}
```

It is easy to transform values from one dimensional base to the other in order to fill the missing values in the table row by row. Starting with the submarine's velocity of `5 m/s` in the application case, the velocity of the scaled model should be selected to `109.6 m/s` in order to behave physically similar: 
```julia
julia> V_model = change_basis(5u"m/s", app_basis, model_basis)
109.58 m s^-1
```
This model velocity is 21.9 times higher than the velocity of the submarine in the real application:
```julia
julia> change_basis(u"m/s", app_basis, model_basis)
21.916172771074063
```

Similarly, the factor needed to convert any given force in Newton from the case of the scaled submarine model to the real application can also be calculated: 
```julia
julia> change_basis(u"N", model_basis, app_basis)
3.4313
```

A similar transformation can also be used to determine the power of the real submarine from experiments with the scaled submarine model. As stated above, the mechanical power of the submarine can generally be expressed as ``P = F \cdot V``. If an exemplary force of `1 N` is measured in the scaled experiment with a suitably scaled velocity of `109.6 m/s` as calculated above, this equates to a mechanical power of `1 N * 109.6 m/s = 109.6 W`. This value can be transformed back to obtain the mechanical power of the real submarine while traveling with `V = 5 m/s`:  
```julia
julia> P_app_fac = change_basis(109.58u"W", model_basis, app_basis)
17.156144908079394 W
```
Of course, at a velocity of ``V = 109.6 \mathrm{m}/\mathrm{s}`` during the scaled experiment, the obtained force can easily be higher than `1 N`. Generally, the mechanical power of the real submarine scales proportional to the measured force with a factor ``17.156144908079394 \mathrm{W}/\mathrm{N}``. 

## Removing and restoring dimensions
Arguably the two most important functions of this package are `dimensionless()` and `dimensionful()`. For a given dimensional basis, these two functions can be used to convert single or multiple variables from a formulation with dimensions to a corresponding dimensionless formulation and reverse: 
```julia
julia> dim_vars = (1u"cm/g", 1u"mm/s")
(1 cm g^-1, 1 mm s^-1)

julia> dimless_vars = dimensionless.(dim_vars, model_basis)
(6.2375, 24.95)

julia> dimensionful.(dimless_vars, (u"cm/g", u"mm/s"), model_basis)
(1.0 cm g^-1, 1.0 mm s^-1)
```
