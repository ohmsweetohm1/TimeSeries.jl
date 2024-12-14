using Test
using MarketData
using TimeSeries
using Dates
using Statistics

@testset "retime" begin

    @testset "single column" begin
        new_timestamps = collect(Dates.Date(2000):Dates.Week(1):Dates.Date(2001))

        funcs = [mean, sum, minimum, maximum, last]
        downsamples = [TimeSeries.Mean(), TimeSeries.Sum(), TimeSeries.Min(), TimeSeries.Max(), TimeSeries.Last()]
        @testset for (func, downsample) in zip(funcs, downsamples)
            cl_new = retime(cl, new_timestamps; upsample=TimeSeries.Linear(), downsample)

            @test timestamp(cl_new) == new_timestamps

            # extrapolation
            @test values(cl_new[1, :Close]) == values(cl[1, :Close])

            # aggregation
            idx = new_timestamps[2] .<= timestamp(cl) .< new_timestamps[3]
            @test func(values(cl[:Close][idx])) == values(cl_new[:Close][2])[1]
        end
    end

    @testset "single column interpolation" begin
        new_timestamps = collect(Dates.DateTime(2000):Dates.Hour(1):Dates.DateTime(2001))

        upsamples = [TimeSeries.Linear(), TimeSeries.Previous(), TimeSeries.Next(), TimeSeries.Nearest()]
        @testset for upsample in upsamples
            cl_new = retime(cl, new_timestamps; upsample)

            @test timestamp(cl_new) == new_timestamps

            # TODO: real tests
        end
    end

    @testset "single column extrapolate" begin
        new_timestamps = collect(Dates.DateTime(2000):Dates.Hour(1):Dates.DateTime(2001))

        cl_new = retime(cl, new_timestamps; extrapolate=TimeSeries.FillConstant(0.0))
        @test timestamp(cl_new) == new_timestamps
        @test values(cl_new[:Close][1])[1] == 0.0

        cl_new = retime(cl, new_timestamps; extrapolate=TimeSeries.NearestExtrapolate())
        @test timestamp(cl_new) == new_timestamps
        @test values(cl_new[:Close][1])[1] == values(cl[:Close][1])[1]

        cl_new = retime(cl, new_timestamps; extrapolate=TimeSeries.MissingExtrapolate())
        @test timestamp(cl_new) == new_timestamps
        @test all(ismissing.(values(cl_new[:Close][1])))

        cl_new = retime(cl, new_timestamps; extrapolate=TimeSeries.NaNExtrapolate())
        @test timestamp(cl_new) == new_timestamps
        @test all(isnan.(values(cl_new[:Close][1])))
    end

    @testset "multi column" begin
        new_timestamps = collect(Dates.Date(2000):Dates.Week(1):Dates.Date(2001))

        funcs = [mean, sum, minimum, maximum, last]
        downsamples = [TimeSeries.Mean(), TimeSeries.Sum(), TimeSeries.Min(), TimeSeries.Max(), TimeSeries.Last()]
        @testset for (func, downsample) in zip(funcs, downsamples)

            ohlc_new = retime(ohlc, new_timestamps; upsample=TimeSeries.Linear(), downsample=TimeSeries.Mean())

            @test timestamp(ohlc_new) == new_timestamps

            # extrapolation
            @test values(ohlc_new[1]) == values(ohlc_new[1])

            idx = new_timestamps[2] .<= timestamp(ohlc) .< new_timestamps[3]
            @test mean(values(ohlc[idx]); dims=1) == values(ohlc_new[2])
        end
    end

    @testset "multi column interpolation" begin
        new_timestamps = collect(Dates.DateTime(2000):Dates.Hour(1):Dates.DateTime(2001))

        upsamples = [TimeSeries.Linear(), TimeSeries.Previous(), TimeSeries.Next(), TimeSeries.Nearest()]
        @testset for upsample in upsamples
            ohlc_new = retime(ohlc, new_timestamps; upsample)

            @test timestamp(ohlc_new) == new_timestamps

            # TODO: real tests
        end
    end

    @testset "multi column extrapolate" begin
        new_timestamps = collect(Dates.DateTime(2000):Dates.Hour(1):Dates.DateTime(2001))

        ohlc_new = retime(ohlc, new_timestamps; extrapolate=TimeSeries.FillConstant(0.0))
        @test timestamp(ohlc_new) == new_timestamps
        @test values(ohlc_new[1]) == zeros(1, 4)

        ohlc_new = retime(ohlc, new_timestamps; extrapolate=TimeSeries.NearestExtrapolate())
        @test timestamp(ohlc_new) == new_timestamps
        @test values(ohlc_new[1]) == values(ohlc[1])

        ohlc_new = retime(ohlc, new_timestamps; extrapolate=TimeSeries.MissingExtrapolate())
        @test timestamp(ohlc_new) == new_timestamps
        @test all(ismissing.(values(ohlc_new[1])))

        ohlc_new = retime(ohlc, new_timestamps; extrapolate=TimeSeries.NaNExtrapolate())
        @test timestamp(ohlc_new) == new_timestamps
        @test all(isnan.(values(ohlc_new[1])))
    end

    @testset "single column with missing" begin
        new_timestamps = collect(Dates.Date(2000):Dates.Week(1):Dates.Date(2001))
        # corrupt some values
        cl_missing = TimeArray(
            timestamp(cl),
            let vals = convert(Vector{Union{Float64,Missing}}, copy(values(cl)))
                vals[rand(1:length(vals), 100)] .= missing
                vals
            end,
            colnames(cl),
        )

        cl_new = retime(cl_missing, new_timestamps; upsample=TimeSeries.Linear(), downsample=TimeSeries.Mean())
    end

    @testset "single column with NaN" begin
        new_timestamps = collect(Dates.Date(2000):Dates.Week(1):Dates.Date(2001))
        # corrupt some values
        cl_missing = TimeArray(
            timestamp(cl),
            let vals = copy(values(cl))
                vals[rand(1:length(vals), 100)] .= NaN
                vals
            end,
            colnames(cl),
        )

        cl_new = retime(cl_missing, new_timestamps; upsample=TimeSeries.Linear(), downsample=TimeSeries.Mean())
    end

end