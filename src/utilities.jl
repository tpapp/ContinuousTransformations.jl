"""
From 

```julia
begin
    symbol1 => value1
    symbol2 => value2
    ...
end
```

return the `symbol, value` pairs as a `Dict`.

 Check for duplicate symbols and non-symbol keys.
"""
function block_key_value_dict(form)
    @capture form begin expressions__ end
    dict = Dict{Symbol,Any}()
    for expression in expressions
        @capture expression key_ = value_
        @argcheck isa(key, Symbol) "Key is not a symbol in $(expression)."
        @argcheck !haskey(dict, key) "Duplicate key in $(expression)."
        dict[key] = value
    end
    dict
end

"""
Define an `isapprox` method, comparing the given fields in type `T`.
"""
macro define_isapprox(T, fields...)
    body = foldl((a, b) -> :($(a) && $(b)),
                 :(isapprox(x.$(field), y.$(field); rtol=rtol, atol=atol))
                 for field in fields)
    quote
        function Base.isapprox(x::$T, y::$T; rtol::Real=sqrt(eps), atol::Real=0)
            $(body)
        end
    end
end
