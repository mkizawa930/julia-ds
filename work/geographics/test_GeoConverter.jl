
include("GeoConverter.jl")

using .GeoConverter

lat, lon = GeoConverter.xy2latlon(11573.375, 22694.980, 33.0, 131.0)

x, y = GeoConverter.latlon2xy(36.103774791666666, 140.08785504166664, 36., 139+50/60)
