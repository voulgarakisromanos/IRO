"""
    is_ground_facility_visible(sat_r_e::AbstractVector, gf_lat::Number, gf_lon::Number, gf_h::Number, θ::Number) -> Bool
    is_ground_facility_visible(sat_r_e::AbstractVector, gf_r_e::AbstractVector, θ::Number) -> Bool
    is_ground_facility_visible(sat_r_ned::AbstractVector, θ::Number) -> Bool

Check if the satellite with position vector `sat_r_e` (ECEF) is inside the visibility circle
of a ground facility with latitude `gf_lat` [rad], longitude `gf_lon` [rad], altitude `gf_h`
(WGS-84) or ECEF position `gf_r_e` [m]. The algorithm considers that the ground station has
visibility to the satellite if its elevation angle is larger than `θ` [rad].

The user can also pass the satellite position represented in the NED (North-East-Down)
reference frame `sat_r_ned` [m] at the ground station location, which increases the
performance since the algorithm performs no reference frame conversion.

# Returns

- `Bool`: `true` if the satellite is inside the visibility circle, or `false` otherwise.
"""
function is_ground_facility_visible(
    sat_r_e::AbstractVector, gf_lat::Number, gf_lon::Number, gf_h::Number, θ::Number
)
    r_ned = ecef_to_ned(sat_r_e, gf_lat, gf_lon, gf_h; translate=true)
    return is_ground_facility_visible(r_ned, θ)
end

function is_ground_facility_visible(
    sat_r_e::AbstractVector, gf_r_e::AbstractVector, θ::Number
)
    gf_wgs84 = ecef_to_geodetic(gf_r_e)
    return is_ground_facility_visible(sat_r_e, gf_wgs84..., θ)
end

function is_ground_facility_visible(r_ned::AbstractVector, minimum_elevation::Number)
    # Check if the satellite is within the minimum elevation supported by the facility.
    # Using the NED vector of the satellite, w.r.t. to the current ground facility, it is
    # sufficient to check the angle θ between `r_ned` and the local vertical (-Z axis),
    # then: el = π / 2 - θ.
    z = r_ned[3]
    r = norm(r_ned)
    θ = acos(-z / r)
    return (π / 2 - θ) > minimum_elevation
end

"""
Check the membership of `p` in `poly` where `poly` is an `AbstractVector` of `AbstractVector`s.
`poly` should have the first and last elements equal.

Returns:
- in = 1
- on = -1
- out = 0

Algorithm by Hao and Sun (2018):
https://doi.org/10.3390/sym10100477
"""
function inpolygon(p, poly)
    k = 0
    xp = p[1]
    yp = p[2]
    PT = eltype(p)
    @inbounds for i in firstindex(poly):(lastindex(poly) - 1)
        v1 = poly[i][2] - yp
        v2 = poly[i + 1][2] - yp
        if v1 < zero(PT) && v2 < zero(PT) || v1 > zero(PT) && v2 > zero(PT)
            continue
        end
        u1 = poly[i][1] - xp
        u2 = poly[i + 1][1] - xp
        f = (u1 * v2) - (u2 * v1)
        if v2 > zero(PT) && v1 <= zero(PT)
            if f > zero(PT)
                k += 1
            elseif iszero(f)
                return -1
            end
        elseif v1 > zero(PT) && v2 <= zero(PT)
            if f < zero(PT)
                k += 1
            elseif iszero(f)
                return -1
            end
        elseif iszero(v2) && v1 < zero(PT)
            iszero(f) && return -1
        elseif iszero(v1) && v2 < zero(PT)
            iszero(f) && return -1
        elseif iszero(v1) && iszero(v2)
            if u2 <= zero(PT) && u1 >= zero(PT)
                return -1
            elseif u1 <= zero(PT) && u2 >= zero(PT)
                return -1
            end
        end
    end
    iszero(k % 2) && return 0
    return 1
end
