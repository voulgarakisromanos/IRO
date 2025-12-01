abstract type AbstractPOI end

"""
    POI

Type alias for `Union` of all available point of interest composite types.
"""

JSON3.StructType(::Type{AbstractPOI}) = JSON3.AbstractType()
JSON3.subtypekey(::Type{AbstractPOI}) = :type

struct GroundStation{T<:Number} <: AbstractPOI
    min_elevation_angle::T
    latitude::T
    longitude::T
end
JSON3.StructType(::Type{GroundStation}) = JSON3.Struct()

struct Vertex{T<:Number}
    latitude::T
    longitude::T
end
JSON3.StructType(::Type{Vertex}) = JSON3.Struct()
function Base.getindex(v::Vertex, i::Int)
    if i == 1
        return v.latitude
    elseif i == 2
        return v.longitude
    else
        throw(BoundsError(v, i))
    end
end

struct Area{T<:Number} <: AbstractPOI
    vertices::Vector{Vertex{T}}
end
JSON3.StructType(::Type{Area}) = JSON3.Struct()

function validate_poly(p)
    if first(p) != last(p)
        error("Polygon should have first and last elements equal")
    end
end

const POI{T<:Number} = Union{GroundStation{T},Area{T}}

struct POIConfig{T<:Number}
    pointsofinterest::Vector{POI{T}}
end

function POI_is_available(gs::GroundStation{T}, env) where {T}
    gs_lla = LLA(gs.latitude, gs.longitude, T(0))
    gs_itrf = ECEF(gs_lla, wgs84)
    return is_ground_facility_visible(env.r_itrf, gs_itrf, deg2rad(gs.min_elevation_angle))
end

function POI_is_available(gs::GroundStation{T}, r_itrf::AbstractVector) where {T}
    gs_lla = LLA(gs.latitude, gs.longitude, T(0))
    gs_itrf = ECEF(gs_lla, wgs84)
    return is_ground_facility_visible(r_itrf, gs_itrf, deg2rad(gs.min_elevation_angle))
end

function POI_is_available(area::Area{T}, env) where {T} # https://github.com/JuliaGeometry/PolygonOps.jl/blob/master/src/inpolygon.jl
    polygon = area.vertices
    r_geodetic = env.r_geodetic # POIs are given in geodetic coordinates
    groundtrack = SVector{2,T}(r_geodetic.lat, r_geodetic.lon) # in rad
    return inpolygon(groundtrack, polygon) != T(0)
end

function POI_is_available(area::Area{T}, r_itrf::AbstractVector) where {T}
    r_geodetic = ecef_to_geodetic(r_itrf)
    polygon = area.vertices
    groundtrack = SVector{2,T}(r_geodetic[1], r_geodetic[2]) # lat, lon, in rad
    return inpolygon(groundtrack, polygon) != T(0)
end

"""
Takes an `ODEProblem` type.

Returns a `Vector{Bool}` indicating whether each point of interest is available.

Ground stations are available when θ > min elevation angle.
Areas are available when the satellite ground track is inside the polygon.
"""
function POI_availability(int)
    sctx = get_sctx(int)
    poi = sctx.pointsofinterest
    env = sctx.env
    return [POI_is_available(p, env) for p in poi]
end
