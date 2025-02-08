using SoilDifferentialEquations
using Ipaper, RTableTools
import Ipaper


function Ipaper.approx(x::AbstractVector, A::AbstractMatrix, xout::AbstractVector)
  ntime = size(A, 1)
  yout = zeros(ntime, length(xout))
  for i = 1:ntime
    y = @view A[i, :]
    yout[i, :] .= approx(x, y, xout)
  end
  yout
end


begin
  d = fread("./test/SM_China/SM_J1193.csv")

  dates = d[:, 1]
  _depths = [10, 20, 30, 40, 50, 60, 80, 100]
  _A = d[:, 2:end] ./ 100 |> Matrix |> drop_missing
  depths = 10:10:100.
  A = approx(_depths, _A, depths)
end
