module GlobalOptions
# 用 Dict 来存储参数
using Parameters
@with_kw mutable struct Options
  method_retention::String = "Campbell"
  method_solve::String = "Bonan"
  same_layer::Bool = true
  ibeg::Int = 2
  yobs::Union{Nothing, AbstractMatrix{Float64}} = nothing
  θ_surf::Union{Nothing, Vector{Float64}} = nothing
end

options = Options()

# 函数用于修改选项
function set_option!(; kw...)
  for (key, value) in kw
    setfield!(options, key, value)
    # setindex!(options, value, key)
    # options[key] = value
  end
end

# 函数用于获取选项
function get_option(key)
  return options[key]
end

export Options, set_option!, get_option
end
