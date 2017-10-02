using Fluxion
using Base.Test


@testset "Helpers" begin
    @testset "Sequence" begin
        @testset "correctness" begin

            seq = Fluxion.event_sequence(Dict(
                :a => [x / 2 for x in 0:10],
                :b => 0:4,
                :c => [1.5, 2.5, 4]
                ))

            res = []

            t = 0
            while true
                events, seq = Fluxion.event_sequence_next_until(t, seq)
                push!(res, (t, events))
                t += 1
                if Fluxion.event_sequence_done(seq)
                    break
                end
            end

            res_ref = [
                (0, [(0.0, [:a, :b])]),
                (1, [(0.5, [:a]), (1.0, [:a, :b])]),
                (2, [(1.5, [:a, :c]), (2.0, [:a, :b])]),
                (3, [(2.5, [:a, :c]), (3.0, [:a, :b])]),
                (4, [(3.5, [:a]), (4.0, [:a, :b, :c])]),
                (5, [(4.5, [:a]), (5.0, [:a])])]

            @test res == res_ref

        end
    end
end
