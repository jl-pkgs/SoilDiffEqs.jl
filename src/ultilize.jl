using Ipaper: approx


function find_ibeg(z_bound_top::FT, zs_center::V) where {FT<:Real,V<:AbstractVector{FT}}
    surface_layer_idx = findfirst(==(z_bound_top), zs_center)
    isnothing(surface_layer_idx) && error("z_bound_top=$z_bound_top not found in zs_center")
    surface_layer_idx + 1 # ibeg
end


function interp_data_depths(A::M, z::V, zout::V) where {
    T<:Real,M<:AbstractMatrix{T},V<:AbstractVector{T}}

    ntime = size(A, 1)
    yout = zeros(ntime, length(zout))
    @inbounds for i in 1:ntime
        yout[i, :] .= approx(z, view(A, i, :), zout)
    end
    yout
end

export find_ibeg, interp_data_depths
