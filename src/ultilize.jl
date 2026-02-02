using Ipaper: approx


function interp_data_depths(A::M, z::V, zout::V) where {
    T<:Real,M<:AbstractMatrix{T},V<:AbstractVector{T}}

    ntime = size(A, 1)
    yout = zeros(ntime, length(zout))
    @inbounds for i in 1:ntime
        yout[i, :] .= approx(z, view(A, i, :), zout)
    end
    yout
end

export interp_data_depths
