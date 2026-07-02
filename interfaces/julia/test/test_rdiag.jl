using Cantera, Test
# These names are exported once src/Cantera.jl is updated; import explicitly so
# this file also runs standalone before the export list is edited.
import Cantera: ReactionPathDiagram, build!, get_dot, get_data, get_log,
    threshold, set_threshold!, bold_threshold, set_bold_threshold!,
    normal_threshold, set_normal_threshold!, label_threshold, set_label_threshold!,
    scale, set_scale!, arrow_width, set_arrow_width!, flow_type, set_flow_type!,
    title, set_title!, font, set_font!, show_details, set_show_details!,
    bold_color, set_bold_color!, normal_color, set_normal_color!,
    dashed_color, set_dashed_color!, dot_options, set_dot_options!, display_only!

@testset "ReactionPathDiagram" begin
    gas = Solution("gri30.yaml")
    set_TPX!(gas, 2000.0, one_atm, "CH4:1, O2:2, N2:7.52")

    d = ReactionPathDiagram(gas, "C")
    set_threshold!(d, 0.01)
    set_flow_type!(d, "NetFlow")
    build!(d)

    dot = get_dot(d)
    @test occursin("digraph", dot)          # valid graphviz
    @test !isempty(get_data(d))
    @test threshold(d) ≈ 0.01
    @test flow_type(d) == "NetFlow"

    # a few more property round-trips
    set_bold_threshold!(d, 0.5)
    @test bold_threshold(d) ≈ 0.5
    set_scale!(d, 2.0)
    @test scale(d) ≈ 2.0
    set_arrow_width!(d, -1.0)
    @test arrow_width(d) ≈ -1.0
    set_title!(d, "carbon flow")
    @test title(d) == "carbon flow"
    set_font!(d, "Courier")
    @test font(d) == "Courier"
    set_show_details!(d, true)
    @test show_details(d) == true
    set_bold_color!(d, "red")
    @test bold_color(d) == "red"

    # display_only!/all should not error
    display_only!(d, 1)
    display_only!(d, :all)
    build!(d)
    @test occursin("digraph", get_dot(d))

    close!(d)
    @test occursin("closed", sprint(show, d))
end
