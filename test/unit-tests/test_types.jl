using BondGraph
using BondGraph: t, D, bg, j0, j1, op, bgeffort, bgflow
using BondGraph: set_bg_metadata, get_bg_junction, update_mtk_con
import ModelingToolkit: VariableConnectType
import Symbolics: unwrap, wrap

using Test
using Symbolics

@testset "set_bg_metadata" begin
    s = @variables x(t)
    modified_s = set_bg_metadata(s[1], op)
    @test getmetadata(unwrap(modified_s), bg) == op
end

@testset "get_bg_junction" begin
    s = @variables x(t)
    s_with_metadata = set_bg_metadata(s[1], j0)
    junction = get_bg_junction(s_with_metadata)
    @test junction == j0
end

@testset "update_mtk_con" begin
    s = @variables x(t)
    updated_s = update_mtk_con(s[1], op)
    @test getmetadata(unwrap(updated_s), VariableConnectType) == op
end