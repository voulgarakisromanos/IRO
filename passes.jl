using Dates
using SatelliteToolbox
using SatelliteToolbox: Propagators, rv_to_kepler, r_eci_to_ecef, GCRF, ITRF, julian2datetime, true_to_mean_anomaly, mean_to_true_anomaly, KeplerianElements
using StaticArrays
using Geodesy: LLA, ECEF, wgs84
using LinearAlgebra: norm
using Interpolations

include("poi.jl")
include("poi_visibility.jl")

const EopType{T<:Number} = EopIau1980{
    Interpolations.Extrapolation{
        T,
        1,
        Interpolations.GriddedInterpolation{
            T,
            1,
            Vector{T},
            Interpolations.Gridded{
                Interpolations.Linear{Interpolations.Throw{Interpolations.OnGrid}}
            },
            Tuple{Vector{T}},
        },
        Interpolations.Gridded{
            Interpolations.Linear{Interpolations.Throw{Interpolations.OnGrid}}
        },
        Interpolations.Flat{Nothing},
    },
}

const IERS_EOP = Base.OncePerProcess{EopType{Float64}}() do
    # We're not using IAU2000A because of sun_position_mod, etc.
    # TODO: consider https://github.com/JuliaSpace/SatelliteToolboxCelestialBodies.jl#rationale
    eop = fetch_iers_eop()
    return eop
end

"""
    get_passes(ke, jd, poi, duration, dt=60)

Returns `Vector{Vector{@NamedTuple{start_time::Float64, duration::Float64}}}`
- Outer vector: One element per point of interest
- Inner vector: All completed passes for that specific POI
- Each pass contains `start_time` (standard JD) and `duration` (seconds)
"""
function get_passes(ke, jd, poi, duration, dt=60)
    orbp = Propagators.init(Val(:J2), ke)
    passes_vec = [
        @NamedTuple{
            start_time::DateTime,
            duration::Float64,
            ke::KeplerianElements{Float64,Float64},
            ritrf::SVector{3,Float64},
        }[] for _ in poi
    ]
    prev_availability = [false for _ in poi]
    current_pass_start = [0.0 for _ in poi]
    current_pass_ke = [ke for _ in poi]
    current_pass_ritrf = [zero(SVector{3,Float64}) for _ in poi]

    step_jd = dt/86400
    final_jd = jd + (duration/86400)
    while jd < final_jd
        r, v = Propagators.step!(orbp, dt)
        R_gcrf2itrf = r_eci_to_ecef(GCRF(), ITRF(), jd, IERS_EOP())
        r_itrf = R_gcrf2itrf * r
        for i in eachindex(poi)
            current_availability = POI_is_available(poi[i], r_itrf)
            if current_availability && !prev_availability[i]
                current_pass_start[i] = jd
                current_pass_ke[i] = rv_to_kepler(r, v, jd)
                current_pass_ritrf[i] = r_itrf
            end
            if !current_availability && prev_availability[i]
                pass_duration = jd - current_pass_start[i]
                push!(
                    passes_vec[i],
                    (
                        start_time=julian2datetime(current_pass_start[i]),
                        duration=round(pass_duration * 86400),
                        ke=current_pass_ke[i],
                        ritrf=current_pass_ritrf[i],
                    ),
                )
            end
            prev_availability[i] = current_availability
        end
        jd += step_jd
    end
    return passes_vec
end
