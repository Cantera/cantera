# ReactionPathDiagram: a façade over the CLib `rdiag_*` reaction-path
# diagram API (ctrdiag.h).  A diagram holds a handle into the CLib rdiag cabinet
# and a reference to the Solution it was built from (to keep the underlying
# Kinetics manager alive for the diagram's lifetime).

"""
    ReactionPathDiagram(gas::Solution, element::AbstractString)

Reaction-path diagram tracing the flow of a chemical `element` (e.g. `"C"`,
`"H"`, `"O"`, `"N"`) through the reaction network of `gas` at its current state.

Set options (see [`set_threshold!`](@ref), [`set_flow_type!`](@ref), ...), then
call [`build!`](@ref) and retrieve the result with [`get_dot`](@ref) (Graphviz),
[`get_data`](@ref) or [`get_log`](@ref).

```julia
gas = Solution("gri30.yaml")
set_TPX!(gas, 2000.0, one_atm, "CH4:1, O2:2, N2:7.52")
d = ReactionPathDiagram(gas, "C")
set_threshold!(d, 0.01)
set_flow_type!(d, "NetFlow")
build!(d)
print(get_dot(d))
```
"""
mutable struct ReactionPathDiagram <: CanteraObject
    handle::Int32
    solution::Solution   # keep the Kinetics manager alive for the diagram
    closed::Bool
end

function ReactionPathDiagram(gas::Solution, element::AbstractString)
    kin = _kinetics_handle(gas)
    h = check(LibCantera.rdiag_newReactionPathDiagram(kin, element))
    d = ReactionPathDiagram(h, gas, false)
    finalizer(close!, d)
    return d
end

"""
    close!(d::ReactionPathDiagram)

Release the Cantera-side object backing `d`.  Safe to call multiple times; the
finalizer calls this automatically.
"""
function close!(d::ReactionPathDiagram)
    d.closed && return d
    d.closed = true
    try
        LibCantera.rdiag_del(d.handle)
    catch
    end
    return d
end

function Base.show(io::IO, d::ReactionPathDiagram)
    if d.closed
        print(io, "ReactionPathDiagram(<closed>)")
    else
        print(io, "ReactionPathDiagram(handle=", d.handle, ")")
    end
end

# ---- build / output ---------------------------------------------------------

"""
    build!(d::ReactionPathDiagram)

Build (or rebuild) the diagram from the current state of the owning Solution.
Call after setting options and before retrieving output.
"""
function build!(d::ReactionPathDiagram)
    check(LibCantera.rdiag_build(d.handle))
    return d
end

"Graphviz `dot` description of the diagram (call [`build!`](@ref) first)."
get_dot(d::ReactionPathDiagram) =
    get_string((n, b) -> LibCantera.rdiag_getDot(d.handle, n, b))

"Raw flow data underlying the diagram (call [`build!`](@ref) first)."
get_data(d::ReactionPathDiagram) =
    get_string((n, b) -> LibCantera.rdiag_getData(d.handle, n, b))

"Diagnostic log accumulated while building the diagram."
get_log(d::ReactionPathDiagram) =
    get_string((n, b) -> LibCantera.rdiag_getLog(d.handle, n, b))

# ---- scalar (double) properties ---------------------------------------------

for (getter, setter, cget, cset) in (
        (:threshold,        :set_threshold!,        :rdiag_threshold,       :rdiag_setThreshold),
        (:bold_threshold,   :set_bold_threshold!,   :rdiag_boldThreshold,   :rdiag_setBoldThreshold),
        (:normal_threshold, :set_normal_threshold!, :rdiag_normalThreshold, :rdiag_setNormalThreshold),
        (:label_threshold,  :set_label_threshold!,  :rdiag_labelThreshold,  :rdiag_setLabelThreshold),
        (:scale,            :set_scale!,            :rdiag_scale,           :rdiag_setScale),
        (:arrow_width,      :set_arrow_width!,      :rdiag_arrowWidth,      :rdiag_setArrowWidth),
    )
    @eval begin
        $getter(d::ReactionPathDiagram) = checkd(LibCantera.$cget(d.handle))
        function $setter(d::ReactionPathDiagram, v::Real)
            check(LibCantera.$cset(d.handle, Float64(v)))
            return d
        end
    end
end

# ---- string properties ------------------------------------------------------

for (getter, setter, cget, cset) in (
        (:flow_type,   :set_flow_type!,   :rdiag_flowType,   :rdiag_setFlowType),
        (:title,       :set_title!,       :rdiag_title,      :rdiag_setTitle),
        (:font,        :set_font!,        :rdiag_font,       :rdiag_setFont),
        (:bold_color,  :set_bold_color!,  :rdiag_boldColor,  :rdiag_setBoldColor),
        (:normal_color,:set_normal_color!,:rdiag_normalColor,:rdiag_setNormalColor),
        (:dashed_color,:set_dashed_color!,:rdiag_dashedColor,:rdiag_setDashedColor),
        (:dot_options, :set_dot_options!, :rdiag_dotOptions, :rdiag_setDotOptions),
    )
    @eval begin
        $getter(d::ReactionPathDiagram) =
            get_string((n, b) -> LibCantera.$cget(d.handle, n, b))
        function $setter(d::ReactionPathDiagram, v::AbstractString)
            check(LibCantera.$cset(d.handle, v))
            return d
        end
    end
end

# ---- boolean property -------------------------------------------------------

"Whether reaction details are shown as edge labels."
show_details(d::ReactionPathDiagram) =
    check(LibCantera.rdiag_showDetails(d.handle)) != 0

"Enable or disable showing reaction details as edge labels."
function set_show_details!(d::ReactionPathDiagram, v::Bool)
    check(LibCantera.rdiag_setShowDetails(d.handle, Int32(v)))
    return d
end

# ---- other operations -------------------------------------------------------

"""
    display_only!(d::ReactionPathDiagram, species_index)

Restrict the diagram to flows into/out of a single species (1-based index), or
pass `:all` (or a negative index) to show all species again.
"""
function display_only!(d::ReactionPathDiagram, k::Integer)
    # CLib expects a 0-based index; -1 means "all".
    idx = k < 0 ? Int32(-1) : Int32(k - 1)
    check(LibCantera.rdiag_displayOnly(d.handle, idx))
    return d
end
function display_only!(d::ReactionPathDiagram, s::Symbol)
    s === :all || throw(ArgumentError("expected a species index or :all, got :$s"))
    return display_only!(d, -1)
end
