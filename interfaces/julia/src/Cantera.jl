# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

"""
    Cantera

Julia interface to [Cantera](https://cantera.org), built on top of
Cantera's generated CLib API.

```julia
using Cantera
gas = Solution("gri30.yaml")
set_TPY!(gas, 1200.0, one_atm, "CH4:1.0, O2:2.0, N2:7.52")
temperature(gas), pressure(gas), net_production_rates(gas)
```
"""
module Cantera

include("LibCantera.jl")
using .LibCantera: LibCantera, libcantera

include("errors.jl")
include("handles.jl")
include("solution.jl")
include("thermo.jl")
include("kinetics.jl")
include("transport.jl")
include("reaction.jl")
include("func1.jl")
include("reactor.jl")
include("reactornet.jl")
include("connectors.jl")
include("onedim.jl")
include("multiphase.jl")
include("rdiag.jl")
include("solutionarray.jl")
include("utils.jl")

# Register Cantera's bundled data directory so mechanisms can be loaded by name.
function __init__()
    d = LibCantera.default_data_directory()
    if d !== nothing
        try
            LibCantera.ct_addDataDirectory(d)
        catch
        end
    end
    return nothing
end

end # module Cantera
