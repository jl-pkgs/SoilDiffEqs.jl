module GlobalOptions
# 用 Dict 来存储参数
using Parameters
@with_kw mutable struct Options
  method_retention::String = "Campbell" # "van_Genuchten", "Campbell"
  method_solve::String = "Bonan"
  same_layer::Bool = true
  ibeg::Int = 2
  yobs::Union{Nothing, AbstractMatrix{Float64}} = nothing
  
  θ0::Union{Nothing, Vector{Float64}} = nothing     # 初始状态
  θ_surf::Union{Nothing, Vector{Float64}} = nothing # 边界层条件
end

options = Options()

# 函数用于修改选项
function set_option!(; kw...)
  for (key, value) in kw
    setfield!(options, key, value)
    # setindex!(options, value, key)
    # options[key] = value
  end
  options
end

# 函数用于获取选项
function get_option(key)
  return options[key]
end

export Options, options, set_option!, get_option

end
