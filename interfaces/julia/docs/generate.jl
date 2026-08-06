# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

# Generate the Sphinx API reference for the Cantera Julia interface.
#
#   julia --project=interfaces/julia interfaces/julia/docs/generate.jl <output-dir>
#
# This script reads the docstrings from the loaded `Cantera` module,
# and writes one MyST page per topic, describing each object with the
# standard `describe` directive so that no Sphinx extension is needed.
# The pages are written into the Sphinx source tree at build time (see doc/SConscript)
# and are not committed.
#
# Loading `Cantera` does not require a working libcantera -- every `ccall` sits inside
# a function body -- but it does require the generated CLib bindings, which `scons
# build` scaffolds.

using Cantera

const OUTDIR = length(ARGS) >= 1 ? ARGS[1] : "build/doc/sphinx/julia"

# Topic pages, in the order they appear in the toctree. Each entry maps a page to the
# interface source files whose objects it documents; this mirrors the per-topic
# organisation of the Python and MATLAB reference pages.
const PAGES = [
    ("phases", "Phases and Solutions", ["solution.jl", "handles.jl"]),
    ("thermo", "Thermodynamic Properties", ["thermo.jl"]),
    ("kinetics", "Chemical Kinetics", ["kinetics.jl", "reaction.jl"]),
    ("transport", "Transport Properties", ["transport.jl"]),
    ("zerodim", "Zero-Dimensional Reactor Networks",
     ["reactor.jl", "reactornet.jl", "connectors.jl", "func1.jl"]),
    ("onedim", "One-Dimensional Reacting Flows", ["onedim.jl"]),
    ("multiphase", "Multiphase Mixtures", ["multiphase.jl"]),
    ("rdiag", "Reaction Path Diagrams", ["rdiag.jl"]),
    ("utilities", "Utilities and Constants",
     ["solutionarray.jl", "utils.jl", "errors.jl"]),
]

"""
An object to be documented: its name, signature, docstring body, and where it was
defined (used for ordering).
"""
struct Entry
    name::String
    label::String
    signature::String
    body::String
    file::String
    line::Int
end

"What sort of object a binding refers to; part of its cross-reference label."
function kind_of(value)
    value isa Type && return "type"
    value isa Function && return "function"
    return "constant"
end

"""
MyST target label for an object, used as the cross-reference anchor.

Sphinx normalises labels by lowercasing them, replacing `_` with `-`, and discarding
characters such as `!`, so the label has to be built to survive that: the `!` of an
in-place variant is spelled out to keep `foo!` distinct from `foo`, and the kind of
object is included because Julia distinguishes a type from a function by case alone, as
with `Reaction` and `reaction`.

The normalisation is applied here rather than left to Sphinx, so that the label is the
identifier that actually appears in the rendered page and the collision check below
compares the labels that Sphinx will really emit.
"""
label_for(name, kind) =
    "jl-" * kind * "-" * lowercase(replace(name, "!" => "-bang", "_" => "-"))

"""
Split a docstring into its leading signature block and the remaining prose.

By convention a Julia docstring begins with one or more indented signature lines,
followed by a blank line. Returns `(signatures, body)`; `signatures` is empty when the
docstring has no such block.
"""
function split_signatures(text)
    lines = split(text, '\n')
    sigs = String[]
    i = 1
    while i <= length(lines) && startswith(lines[i], "    ") &&
            !isempty(strip(lines[i]))
        push!(sigs, strip(lines[i]))
        i += 1
    end
    isempty(sigs) && return (String[], text)
    # Drop the blank line that separates the block from the prose.
    while i <= length(lines) && isempty(strip(lines[i]))
        i += 1
    end
    return (sigs, join(lines[i:end], '\n'))
end

"Signature of a single method, using its declared argument names and types."
function method_signature(name, m)
    _, decls, _, _ = Base.arg_decl_parts(m)
    args = [isempty(t) ? n : "$n::$t" for (n, t) in decls[2:end]]
    return string(name, "(", join(args, ", "), ")")
end

"Argument tuple type of a method, in the form the docsystem keys docstrings by."
function argument_types(m)
    sig = try
        Base.unwrap_unionall(m.sig)
    catch
        return nothing
    end
    length(sig.parameters) < 1 && return nothing
    return try
        Tuple{sig.parameters[2:end]...}
    catch
        nothing
    end
end

"""
Best-effort signature for a docstring that has no signature block of its own.

`key` is the method signature the docsystem attached the docstring to. The matching
method is preferred, so that the rendered signature describes the same method as the
prose beneath it — for a name such as `temperature`, documented once per argument type,
picking any other method would label the text wrongly. Failing that the method with the
most arguments is used, which for an accessor with an in-place variant is the more
informative one.
"""
function derived_signature(name, value, key=Union{})
    if key isa Type && key !== Union{} && value isa Union{Function,Type}
        for m in methods(value)
            argument_types(m) === key && return method_signature(name, m)
        end
        # No method matches (a docstring on an inner constructor, say): fall back to the
        # argument types recorded in the key itself.
        from_key = signature_from_type(name, key)
        from_key === nothing || return from_key
    end
    # A type or a constant is named on its own; only a function gets a call signature.
    value isa Function || return string(name)
    ms = collect(methods(value))
    isempty(ms) && return string(name)
    return method_signature(name, ms[argmax([m.nargs for m in ms])])
end

"""
Signature rendered from the method-signature type a docstring is keyed by.

Julia keys each of a binding's docstrings by the signature of the method it documents,
so a docstring that carries no signature block of its own can still be labelled with the
argument types it applies to. Returns `nothing` when the key is not a method signature,
as for the docstring attached to a type or a constant.
"""
function signature_from_type(name, key)
    key isa Type || return nothing
    key === Union{} && return nothing
    params = try
        Base.unwrap_unionall(key).parameters
    catch
        return nothing
    end
    # The key is the tuple of *argument* types; it does not include the function itself.
    isempty(params) && return nothing
    args = ["::" * string(p) for p in params]
    return string(name, "(", join(args, ", "), ")")
end

"""
Every docstring attached to a binding, as `(key, signatures, body)` triples in source
order.

Julia allows one docstring per method signature, and the interface uses that for names
whose meaning depends on the argument type: `temperature` is documented separately for
phases, mixtures and solution arrays. All of them are collected so that none is dropped
from the reference. `key` is the method signature the docstring is attached to, used to
label a docstring that carries no signature block of its own. Duplicates are left in
place; they can only be recognised once the signatures have been resolved.
"""
function collect_docs(multidoc)
    parts = Tuple{Any,Vector{String},String}[]
    for key in multidoc.order
        sigs, body = split_signatures(join(multidoc.docs[key].text))
        push!(parts, (key, sigs, String(strip(body))))
    end
    return parts
end

"""
Prepare a Markdown docstring body for inclusion in a MyST page.

Julia docstrings are already Markdown, so the body passes through unchanged apart from
Documenter's `[`name`](@ref)` cross-references, which are rewritten as MyST links to
the target emitted for that object. A name that is not documented here -- the `Cantera`
module itself, for instance -- becomes literal text rather than a dangling link.
"""
function prepare_body(body, labels)
    return strip(replace(body, r"\[`([^`]+)`\]\(@ref\)" => function (match)
        name = match[3:findfirst("`](@ref)", match).start - 1]
        haskey(labels, name) || return "`$name`"
        return "[`$name`](#$(labels[name]))"
    end))
end

"Source location of a documented binding, for ordering entries within a page."
function location(docstr)
    path = get(docstr.data, :path, "")
    line = get(docstr.data, :linenumber, 0)
    return (isempty(path) ? "" : basename(String(path)), Int(line))
end

"Collect every documented, exported object of the `Cantera` module."
function collect_entries()
    exported = Set(names(Cantera))
    meta = Base.Docs.meta(Cantera)

    # Labels of the objects that end up on a reference page, and so can be linked to.
    # These are resolved up front because a docstring may refer to an object that is
    # documented later, or on another page.
    labels = Dict{String,String}()
    for binding in keys(meta)
        sym = binding.var
        (sym in exported && sym !== :Cantera && isdefined(Cantera, sym)) || continue
        name = string(sym)
        label = label_for(name, kind_of(getfield(Cantera, sym)))
        if haskey(labels, label)
            # A collision would silently point every reference at one of the objects.
            error("Julia objects '$(labels[label])' and '$name' share the " *
                  "cross-reference label '$label'")
        end
        labels[label] = name
    end
    labels = Dict(name => label for (label, name) in labels)

    entries = Entry[]
    for (binding, multidoc) in meta
        sym = binding.var
        sym in exported || continue
        sym === :Cantera && continue
        isdefined(Cantera, sym) || continue
        value = getfield(Cantera, sym)

        name = string(sym)
        parts = collect_docs(multidoc)

        # The first signature heads the `describe` directive; every other call form,
        # whether it came from a further line of the same docstring or from a separate
        # docstring on another method, is reintroduced in the body as a `julia` block
        # ahead of the prose that documents it. That way a name documented once per
        # argument type keeps all of its meanings.
        #
        # Deduplication happens only once the signatures are resolved: two methods may
        # share the same one-line description while still being distinct call forms
        # worth listing, as `temperature` is for a phase and for a mixture.
        resolved = Tuple{Vector{String},String}[]
        for (key, sigs, body) in parts
            isempty(sigs) && (sigs = [derived_signature(name, value, key)])
            any(r -> r[1] == sigs && r[2] == body, resolved) && continue
            push!(resolved, (sigs, body))
        end

        code(lines) = string("```julia\n", join(lines, '\n'), "\n```")
        sigs, body = resolved[1]
        signature = first(sigs)

        sections = String[]
        length(sigs) > 1 && push!(sections, code(sigs[2:end]))
        isempty(body) || push!(sections, body)
        for (more_sigs, more_body) in resolved[2:end]
            push!(sections, code(more_sigs))
            isempty(more_body) || push!(sections, more_body)
        end
        body = join(sections, "\n\n")

        file, line = location(multidoc.docs[first(multidoc.order)])
        push!(entries, Entry(name, labels[name], signature,
                             prepare_body(body, labels), file, line))
    end
    return entries
end

"Write one topic page."
function write_page(io, title, entries)
    println(io, "# ", title)
    println(io)
    for e in entries
        println(io, "(", e.label, ")=")
        # Four colons so that a fenced code block in the body cannot close the
        # directive early.
        println(io, "::::{describe} ", e.signature)
        isempty(e.body) || println(io, e.body)
        println(io, "::::")
        println(io)
    end
end

function main()
    mkpath(OUTDIR)
    entries = collect_entries()
    documented = 0

    for (page, title, files) in PAGES
        selected = filter(e -> e.file in files, entries)
        sort!(selected, by = e -> (findfirst(==(e.file), files), e.line))
        open(joinpath(OUTDIR, page * ".md"), "w") do io
            write_page(io, title, selected)
        end
        documented += length(selected)
        println("  $(rpad(page * ".md", 20)) $(length(selected)) objects")
    end

    # The list of pages and the toctree live here rather than in the hand-written
    # index, so that a documentation build without Julia can substitute a placeholder
    # (see doc/SConscript) without referring to pages that were never written.
    open(joinpath(OUTDIR, "api-reference.md"), "w") do io
        for (page, _, _) in PAGES
            println(io, "- [](", page, ".md)")
        end
        println(io)
        println(io, "```{toctree}")
        println(io, ":hidden:")
        println(io, ":maxdepth: 1")
        println(io)
        for (page, _, _) in PAGES
            println(io, page)
        end
        println(io, "```")
    end

    # Fail loudly rather than silently dropping objects from the reference if a new
    # source file is added to the interface without being assigned to a page.
    known = Set(Iterators.flatten(files for (_, _, files) in PAGES))
    orphans = unique([e.file for e in entries if !(e.file in known)])
    if !isempty(orphans)
        error("Julia source files not assigned to a reference page: " *
              join(sort(orphans), ", "))
    end
    println("  generated $documented of $(length(entries)) documented objects")
end

main()
