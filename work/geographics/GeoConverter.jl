module GeoConverter

"""
参考記事
- https://qiita.com/sw1227/items/e7a590994ad7dcd0e8ab
- https://sw1227.hatenablog.com/entry/2018/11/30/200702#3-%E7%B5%8C%E7%B7%AF%E5%BA%A6%E3%81%8B%E3%82%89%E5%B9%B3%E9%9D%A2%E7%9B%B4%E8%A7%92%E5%BA%A7%E6%A8%99%E3%81%B8%E3%81%AE%E5%A4%89%E6%8F%9B%E5%BC%8F
"""

using OffsetArrays
using PyCall

# numpyで書き換える用
# np = pyimport("numpy")
# sin = np.sin
# cos = np.cos
# tan = np.tan
# sinh = np.sinh
# cosh = np.cosh
# atan = np.arctan
# atanh = np.arctanh

const m0 = 0.9999
const a = 6378137.0 
const F = 298.257222101

const n = 1.0 / (2F - 1)

const Α = OffsetArray([
    1 + n^2 / 4 + n^4 / 64,
    - (3/2) * (n - n^3 / 8 - n^5 / 64),
    (15/16) * (n^2 - n^4 / 4),
    -(35/48) * (n^3 - 5/16 * n^5),
    (315/512) * n^4,
    -(693/1280) * n^5,
], 0:5)

const α = [
    (1/2)*n - (2/3)*n^2 + (5/16)*n^3 + (41/180)*n^4 - (127/288)*n^5,
    (13/48)*n^2 - (3/5)*n^3 + (557/1440)*n^4 + (281/630)*n^5,
    (61/240)*n^3 - (103/140)*n^4 + (15061/26880)*n^5,
    (49561/161280)*n^4 - (179/168)*n^5,
    (34729/80640)*n^5,
]

const β = [
    (1/2)*n - (2/3)*n^2 + (37/96)*n^3 - (1/360)*n^4 - (81/512)*n^5,
    (1/48)*n^2 + (1/15)*n^3 - (437/1440)*n^4 + (46/105)*n^5,
    (17/480)*n^3 - (37/840)*n^4 - (209/4480)*n^5,
    (4397/161280)*n^4 - (11/504)*n^5,
    (4583/161280)*n^5,
]

const δ = [
    2*n - (2/3)*n^2 - 2*n^3 + (116/45)*n^4 + (26/45)*n^5 - (2854/675)*n^6,
    (7/3)*n^2 - (8/5)*n^3 - (227/45)*n^4 + (2704/315)*n^5 + (2323/945)*n^6,
    (56/15)*n^3 - (136/35)*n^4 - (1262/105)*n^5 + (73814/2835)*n^6,
    (4279/630)*n^4 - (332/35)*n^5 ^ (399572/14175)*n^6,
    (4174/315)*n^5 - (144838/6237)*n^6,
    (601676/22275)*n^6,
]


"""
XY平面座標系から地理座標系(緯度経度)に変換する
"""
function xy2latlon(x, y, ϕ0, λ0)
    global Α, β, δ

    ϕ0 = deg2rad(ϕ0)
    λ0 = deg2rad(λ0)

    Ᾱ = (m0 * a) / (1 + n) * Α[0]
    S̄ = ((m0 * a) / (1 + n)) * (Α[0] * ϕ0 + sum(Α[j] * sin(2j * ϕ0) for j in 1:5) )

    ξ = (x + S̄) / Ᾱ
    η = y / Ᾱ

    ξ′ = ξ - sum( β[j] * sin(2j*ξ) * cosh(2j*η) for j in 1:5 ) 
    η′ = η - sum( β[j] * cos(2j*ξ) * sinh(2j*η) for j in 1:5 )

    χ = asin( sin(ξ′) / cosh(η′) )

    lat = χ + sum( δ[j] * sin(2*χ*j) for j in 1:6 )
    lon = λ0 + atan( sinh(η′) / cos(ξ′) )

    rad2deg(lat), rad2deg(lon)
end



"""
地理座標系(緯度経度)からXY平面座標系に変換する

(ϕ0, λ0): 基準平面座標の中心を表す緯度と経度 (degree)
"""
function latlon2xy(ϕ, λ, ϕ0, λ0)
    global Α, α
    
    ϕ = deg2rad(ϕ)
    λ = deg2rad(λ)
    ϕ0 = deg2rad(ϕ0)
    λ0 = deg2rad(λ0)

    Ᾱ = (m0 * a) / (1 + n) * Α[0]
    S̄ = ((m0 * a) / (1 + n)) * (Α[0] * ϕ0 + sum(Α[j] * sin(2j * ϕ0) for j in 1:5) )

    λc = cos(λ - λ0)
    λs = sin(λ - λ0)

    t = sinh( atanh(sin(ϕ)) - ((2*sqrt(n)) / (1 + n)) * atanh((2*sqrt(n)/(1+n)) * sin(ϕ)) )
    t̄ = sqrt(1 + t^2)
    
    ξ′ = atan(t / λc)
    η′ = atanh(λs / t̄)

    x = Ᾱ * ( ξ′ + sum(α[j] * sin(2j*ξ′) * cosh(2j*η′) for j in 1:5) ) - S̄
    y = Ᾱ * ( η′ + sum(α[j] * cos(2j*ξ′) * sinh(2j*η′) for j in 1:5) )

    x, y
end


end # module