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
function set_option!(opt::Options; kw...)
  for (key, value) in kw
    hasproperty(opt, key) || throw(ArgumentError("未知选项: $(key)"))
    setproperty!(opt, key, value)
  end
  opt
end

set_option!(; kw...) = set_option!(options; kw...)

# 函数用于获取选项
get_option(key::Symbol) = getproperty(options, key)
get_option(key::AbstractString) = getproperty(options, Symbol(key))

export Options, options, set_option!, get_option

end
