using DataFrames
using Gadfly
using Distributions

function g(E)
  return (1./(sqrt(2*pi)*1.6e-9^3))*exp(-(E^2)/2.)
end

function s(g)
  prefac = (1./(sqrt(2*pi)*1.6e-9^3))
  return sqrt(2.*log(prefac/g))
end

Erand = DataFrame()

Erand[:evals] = collect(-5:0.001:5)
Erand[:gvals] = Array(Float64,size(Erand[:evals],1))
Erand[:svals] = Array(Float64,size(Erand[:evals],1))
Erand[:rands] = rand(Normal(0,1), size(Erand[:evals],1))
Erand[:ehist] = Array(Float64,size(Erand[:evals],1))
Erand[:label] = "E_rand"

prefac = (1./(sqrt(2*pi)*1.6e-9^3))

for i = 1:size(Erand[:evals],1)
  Erand[:gvals][i] = g(Erand[:evals][i])
  Erand[:svals][i] = prefac*exp(-(Erand[:rands][i]^2)/2)
  Erand[:ehist][i] = Erand[:evals]
end

p1 = plot(Erand,x="evals",y="gvals",Geom.line,Coord.Cartesian(xmin=-5.,xmax=5.))
p2 = plot(Erand,x="svals",color="label",Geom.histogram)
p3 = plot(Erand,x="rands",y="svals",Geom.line)
push!(p1,layer(Erand,x="rands",y="svals",color="label",Geom.point))

println("Ending script")
