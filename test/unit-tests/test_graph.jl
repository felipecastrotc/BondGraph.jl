using BondGraph
using BondGraph: t, D
import BondGraph: csets2dict, csets2adjmtx

using ModelingToolkit
import ModelingToolkit: generate_connection_set
import ModelingToolkit: ConnectionElement, namespaced_var

using SparseArrays
using Test
using Symbolics

name = :test
split_str = "₊"

# Create a very basic system to generate a sample_csets
npin = :pin
npout = :pout
@named pin = Power(name=npin)
@named pout = Power(name=npout)

# Create the connection equations
eqs = [connect(pin, pout)]
true_adjmtx = sparse([0 0 0 1; 0 0 1 0; 0 0 0 0; 0 0 0 0])

# Build the dummy system
@named mdl = ODESystem(eqs, t)
sys = compose(mdl, pin, pout)

# Generate the sample_csets
_, (_, sample_csets) = generate_connection_set(sys, nothing, nothing)


@testset "csets2dict" begin
    @testset "Basic Functionality" begin
        result = csets2dict(sample_csets)
        @test isa(result, Dict{String,ConnectionElement})
    end

    @testset "Conversion and size" begin
        result = csets2dict(sample_csets)
        # Each variable (pin and pout) has two power varaibles e and f
        @test length(result) == length(sample_csets)*2
        # Check if dictionary is built correctly
        for cset in sample_csets
            for ele in cset.set
                k = String(Symbol(namespaced_var(ele)))
                @test result[k] == ele
            end
        end
    end

    @testset "Empty Connection Sets" begin
        result = csets2dict([])
        @test isempty(result)
    end

end

@testset "csets2adjmtx" begin
    sample_str2con = csets2dict(sample_csets)

    @testset "Basic Functionality" begin
        adjmtx = csets2adjmtx(sample_csets, sample_str2con)
        @test isa(adjmtx, SparseMatrixCSC)
    end

    @testset "Filter Flow" begin
        str2con = csets2dict(sample_csets)
        # When using filter the str2con will change
        adjmtx = csets2adjmtx(sample_csets, str2con, filterstr="f", filterflow=true)
        # It should only contain the flow variable
        @test size(adjmtx) == (2, 2)
        @test length(keys(str2con)) == 2
        # Check if the keys are correct
        @test string(pin.f) in keys(str2con)
        @test string(pout.f) in keys(str2con)
    end

    @testset "Filter Effort" begin
        str2con = csets2dict(sample_csets)
        # When using filter the str2con will change
        adjmtx = csets2adjmtx(sample_csets, str2con, filterstr="e", filterflow=true)
        # It should only contain the flow variable
        @test size(adjmtx) == (2, 2)
        @test length(keys(str2con)) == 2
        # Check if the keys are correct
        @test string(pin.e) in keys(str2con)
        @test string(pout.e) in keys(str2con)
    end

    @testset "Matrix Entries" begin
        adjmtx = csets2adjmtx(sample_csets, sample_str2con)
        @test true_adjmtx == adjmtx
    end

end


# TODO: adjmtx2eqs
# TODO: retrieve_node_details
# TODO: gen_con_name

# NOTE: I don't think it is necessary to todo for the generate_graph. If other 
# functions pass it will pass. The main issue is the graphplot function which 
# depends on GraphRecipes.