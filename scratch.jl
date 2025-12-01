using JSON3

include("passes.jl")

# Points of Interest
pointsofinterest = JSON3.read("pointsofinterest.json", POIConfig{Float64}).pointsofinterest

jd = 2458850.0
RAAN = 5.360315179750939 # Radians
argument_of_perigee = 5.4548663422376515 # Radians
eccentricity = 0.001220726589790339
inclination = 1.710422666954443 # Radians
mean_anomaly = -3.054957385063848 # Radians
semi_major_axis = 6863.698878774085 # km

ke = KeplerianElements(
        jd,
        semi_major_axis * 1e3,
        eccentricity,
        inclination,
        RAAN,
        argument_of_perigee,
        mean_to_true_anomaly(eccentricity, mean_anomaly)
    )

passes = get_passes(ke, jd, pointsofinterest, 10 * 86_400, 60)