"""
Kneading.jl's optional KernelAbstractions extension. It loads only when
KernelAbstractions is available and runs batched flow-kneading word integration
on any KernelAbstractions backend, such as `CPU()` or `CUDABackend()`.
"""
module KneadingKernelAbstractionsExt

using KernelAbstractions
using StaticArrays: SVector
using Kneading.FlowKneading: _batch_word

import Kneading.FlowKneading: _launch_flow_words

@kernel function _flow_word_kernel!(bits, lengths, statuses, accepted, terminal_times, events,
    rule, options, @Const(parameters), @Const(states), @Const(tangents), @Const(directions), @Const(flags))
    i = @index(Global)
    word = _batch_word(rule, options, parameters[i], states[i], tangents[i], directions[i],
        flags[i], eltype(bits), events, i)
    bits[i] = word[1]
    lengths[i] = word[2] % Int32
    statuses[i] = word[3]
    accepted[i] = word[4] % Int32
    terminal_times[i] = word[5]
end

function _device(backend, values)
    device = KernelAbstractions.allocate(backend, eltype(values), size(values))
    copyto!(device, values)
    return device
end

_device_zeros(backend, T, dims...) = KernelAbstractions.zeros(backend, T, dims...)

function _launch_flow_words(backend::Backend, rule, options, parameters,
    states::AbstractVector{SVector{N,Float64}}, tangents, directions, flags, ::Type{B}, slots;
    chunk_size = 65536, workgroup_size = nothing) where {N,B}
    count = length(states)
    outputs = (; bits = zeros(B, count), lengths = zeros(Int32, count), statuses = zeros(Int8, count),
        accepted = zeros(Int32, count), terminal_times = zeros(Float64, count))
    if slots > 0
        outputs = merge(outputs, (; times = fill(NaN, slots, count),
            states = fill(zero(SVector{N,Float64}), slots, count),
            tangents = fill(zero(SVector{N,Float64}), slots, count),
            rates = fill(NaN, slots, count), components = fill(NaN, slots, count)))
    end
    groups = something(workgroup_size, backend isa CPU ? 1 : 64)
    kernel! = _flow_word_kernel!(backend, groups)
    for start in 1:chunk_size:count
        range = start:min(start + chunk_size - 1, count)
        m = length(range)
        bits = _device_zeros(backend, B, m)
        lengths = _device_zeros(backend, Int32, m)
        statuses = _device_zeros(backend, Int8, m)
        accepted = _device_zeros(backend, Int32, m)
        terminal_times = _device_zeros(backend, Float64, m)
        events = slots > 0 ? (;
            times = _device_zeros(backend, Float64, slots, m),
            states = _device(backend, fill(zero(SVector{N,Float64}), slots, m)),
            tangents = _device(backend, fill(zero(SVector{N,Float64}), slots, m)),
            rates = _device_zeros(backend, Float64, slots, m),
            components = _device_zeros(backend, Float64, slots, m)) : nothing
        kernel!(bits, lengths, statuses, accepted, terminal_times, events, rule, options,
            _device(backend, parameters[range]), _device(backend, states[range]),
            _device(backend, tangents[range]), _device(backend, directions[range]),
            _device(backend, flags[range]); ndrange = m)
        KernelAbstractions.synchronize(backend)
        outputs.bits[range] .= Array(bits)
        outputs.lengths[range] .= Array(lengths)
        outputs.statuses[range] .= Array(statuses)
        outputs.accepted[range] .= Array(accepted)
        outputs.terminal_times[range] .= Array(terminal_times)
        if slots > 0
            for name in (:times, :states, :tangents, :rates, :components)
                getfield(outputs, name)[:, range] .= Array(getfield(events, name))
            end
        end
    end
    return outputs
end

end
