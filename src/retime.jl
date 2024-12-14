# Abstract types for interpolation, aggregation, and extrapolation methods
abstract type InterpolationMethod end
abstract type AggregationMethod end
abstract type ExtrapolationMethod end

# Interpolation methods
struct Linear <: InterpolationMethod end
struct Previous <: InterpolationMethod end
struct Next <: InterpolationMethod end
struct Nearest <: InterpolationMethod end

# Aggregation methods
struct Mean <: AggregationMethod end
struct Min <: AggregationMethod end
struct Max <: AggregationMethod end
struct Count <: AggregationMethod end
struct Sum <: AggregationMethod end
struct Median <: AggregationMethod end
struct First <: AggregationMethod end
struct Last <: AggregationMethod end

# Extrapolation methods
struct FillConstant <: ExtrapolationMethod
    value::Float64
end
struct NearestExtrapolate <: ExtrapolationMethod end
struct MissingExtrapolate <: ExtrapolationMethod end
struct NaNExtrapolate <: ExtrapolationMethod end


function retime(ta, new_dt::Dates.Period; kwargs...)
    new_timestamps = timestamp(ta)[1]:new_dt:timestamp(ta)[end]
    return retime(ta, new_timestamps; kwargs...)
end

function retime(ta, period::Function; kwargs...)
    new_timestamps = map(i -> first(i), _split(timestamp(ta), period))
    return retime(ta, new_timestamps; kwargs...)
end

function retime(
    ta::TimeArray{T,N,D,A},
    new_timestamps::AbstractVector{DN};
    upsample::InterpolationMethod=Previous(),
    downsample::AggregationMethod=Mean(),
    extrapolate::ExtrapolationMethod=NearestExtrapolate(),
    # TODO: handle missing and NaN
) where {T,N,D,A,DN}
    new_values = __get_new_values(T, length(new_timestamps), size(values(ta), 2), extrapolate) # zeros(T, length(new_timestamps), size(values(ta), 2))
    old_timestamps = convert(Vector{DN}, timestamp(ta))
    old_values = values(ta)
    @views begin
        for col_i in 1:size(old_values, 2)
            _retime!(new_values[:, col_i], old_timestamps, old_values[:, col_i], new_timestamps, upsample, downsample, extrapolate)
        end
    end
    return TimeArray(new_timestamps, new_values, colnames(ta), meta(ta))
end

function _retime!(
    new_values::AbstractVector{AN},
    old_timestamps::AbstractVector{D},
    old_values::AbstractVector{A},
    new_timestamps::AbstractVector{D},
    upsample::InterpolationMethod,
    downsample::AggregationMethod,
    extrapolate::ExtrapolationMethod,
) where {D,AN,A}

    x = Dates.value.(old_timestamps)
    x_min, x_max = extrema(x)
    x_new = Dates.value.(new_timestamps)

    N = length(x_new)

    # check each interval between i and i+1 if there is no or one sample (upsample), more than one sample (downsample)
    for i in 1:N
        if x_new[i] < x_min || x_new[i] > x_max
            # Handle extrapolation
            new_values[i] = _extrapolate(extrapolate, x_new[i], x, old_values)
        else
            idx = if i < N
                _get_idx(x, x_new[i], x_new[i+1])
            else
                # assume that the last interval is the same length as the second to last one
                _get_idx(x, x_new[i], x_new[i] + (x_new[i] - x_new[i-1]))
            end

            if isempty(idx)
                # No original samples lie between x_new[i] and x_new[i+1] --> Upsampling
                new_values[i] = _upsample(upsample, x, old_values, x_new[i])
            elseif length(idx) == 1
                if x_new[i] == x[idx[1]] # directly hit the sample, do not try the upsample method
                    new_values[i] = old_values[idx[1]]
                else
                    # Only one sample found in the interval x_new[i] and x_new[i+1] --> use the upsample method
                    new_values[i] = _upsample(upsample, x, old_values, x_new[i])
                end
            else
                # Multiple samples were found in the interval [x_new[i], x_new[i+1]) --> use the downsample method to get the agglomeration
                new_values[i] = _downsample(downsample, old_values[idx])
            end
        end
    end
    return
end

function __get_new_values(T, N, n, extrapolate)
    zeros(T, N, n)
end
function __get_new_values(T, N, n, extrapolate::MissingExtrapolate)
    zeros(Union{Missing,T}, N, n)
end


function _get_idx(x::AbstractVector{<:Real}, x_left::Real, x_right::Real)
    idx_left = searchsortedfirst(x, x_left) # greater or equal to x_left
    idx_right = searchsortedlast(x, prevfloat(Float64(x_right))) # less to x_right
    return idx_left:idx_right
end

# Extrapolation dispatch
function _extrapolate(m::FillConstant, t_new, x, y)
    return m.value
end

function _extrapolate(::NearestExtrapolate, t_new, x, y)
    idx = argmin(abs.(x .- t_new))
    return y[idx]
end

function _extrapolate(::MissingExtrapolate, t_new, x, y)
    return missing
end

function _extrapolate(::NaNExtrapolate, t_new, x, y)
    return NaN
end

# Interpolation dispatch
function _upsample(::Linear, x_old, old_values, x)
    idx_next = searchsortedfirst(x_old, x) # greater or equal to x
    idx_prev = searchsortedlast(x_old, x) # less or equal to x
    y = if idx_prev == idx_next # avoid division by zero
        old_values[idx_prev]
    else
        old_values[idx_prev] + (x - x_old[idx_prev]) * (old_values[idx_next] - old_values[idx_prev]) / (x_old[idx_next] - x_old[idx_prev])
    end

end

function _upsample(::Previous, x_old, old_values, x)
    idx_prev = searchsortedlast(x_old, x) # less or equal to x
    return old_values[idx_prev]
end

function _upsample(::Next, x_old, old_values, x)
    idx_next = searchsortedfirst(x_old, x) # greater or equal to x
    return old_values[idx_next]
end

function _upsample(::Nearest, x_old, old_values, x)
    idx_next = searchsortedfirst(x_old, x) # greater or equal to x
    idx_prev = searchsortedlast(x_old, x)  # less or equal to x
    y = if idx_prev == idx_next # avoid division by zero
        old_values[idx_prev]
    else
        pos = (x - x_old[idx_prev]) / (x_old[idx_next] - x_old[idx_prev])
        if pos < 0.5
            old_values[idx_prev]
        else
            old_values[idx_next]
        end
    end
    return y
end

# Aggregation dispatch
function _downsample(::Mean, values_in_range)
    return mean(values_in_range)
end

function _downsample(::Min, values_in_range)
    return minimum(values_in_range)
end

function _downsample(::Max, values_in_range)
    return maximum(values_in_range)
end

function _downsample(::Count, values_in_range)
    return count(!ismissing, values_in_range)
end

function _downsample(::Sum, values_in_range)
    return sum(values_in_range)
end

function _downsample(::Median, values_in_range)
    return median(values_in_range)
end

function _downsample(::First, values_in_range)
    return first(values_in_range)
end

function _downsample(::Last, values_in_range)
    return last(values_in_range)
end