using DataStructures: SortedDict


function last_element(iterable)
    # the built-in `last()` is limited to O(1) indexed access collections
    try
        return last(iterable)
    catch e
        if isa(e, MethodError)
            for elem in iterable
            end
            return elem
        end
    end
end


immutable Lookahead
    iterable
    state
    peek :: Nullable{Any}
end


function lookahead(iterable)
    state = start(iterable)
    if done(iterable, state)
        peek = nothing
    else
        peek, state = next(iterable, state)
        peek = Nullable{Any}(peek)
    end
    Lookahead(iterable, state, peek)
end


function lookahead_peek(la::Lookahead)
    get(la.peek)
end


function lookahead_done(la::Lookahead)
    isnull(la.peek)
end


function lookahead_next(la::Lookahead)
    if done(la.iterable, la.state)
        Lookahead(la.iterable, la.state, nothing)
    else
        item, state = next(la.iterable, la.state)
        Lookahead(la.iterable, state, Nullable{Any}(item))
    end
end


immutable EventSequence
    iters :: SortedDict{Symbol, Lookahead}
    last_event
end


function event_sequence_done(eseq::EventSequence)
    length(eseq.iters) == 0
end


function event_sequence(seqs)
    # Assuming:
    # - all iterables produce monotonously increasing values
    iters = SortedDict(((key => lookahead(seq)) for (key, seq) in seqs)...)
    last_event = max((last_element(seq) for (key, seq) in seqs)...)
    EventSequence(iters, last_event)
end


function event_sequence_next_until(max_val, eseq::EventSequence)
    events = [] # Array{Pair{Any, Symbol}, 1}()
    new_lookaheads = Array{Pair{Symbol, Lookahead}, 1}()

    for (key, lookahead) in eseq.iters
        while !lookahead_done(lookahead)

            val = lookahead_peek(lookahead)

            # handling this case separately to ensure values are snapped to max_val
            # when possible (this will help steppers to determine whether or not
            # to interpolate)
            if isapprox(val, max_val)
                val = max_val
            elseif val > max_val
                # leave it for the next iteration
                break
            end

            push!(events, val => key)

            lookahead = lookahead_next(lookahead)
        end

        if !lookahead_done(lookahead)
            push!(new_lookaheads, key => lookahead)
        end
    end

    grouped_events = [] # Array{Pair{Any, Array{Symbol, 1}}, 1}()

    # Group events with close enough times
    for (val, key) in sort(events)
        if length(grouped_events) == 0
            push!(grouped_events, (val, [key]))
        else
            last_val = grouped_events[end][1]
            if isapprox(val, last_val)
                push!(grouped_events[end][2], key)
            else
                push!(grouped_events, (val, [key]))
            end
        end
    end

    grouped_events, EventSequence(SortedDict(new_lookaheads), eseq.last_event)
end
