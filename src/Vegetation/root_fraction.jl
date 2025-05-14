## 设置根系分布情况
"""
- z: in cm, 向下为正
"""
root_fraction(z_cm::AbstractVector; β) = β .^ z_cm # 

root_fraction(soil::Soil; β=0.943) = root_fraction(-soil.z_cm[0:N]; β)

## 根据最大根系深度，推求beta
