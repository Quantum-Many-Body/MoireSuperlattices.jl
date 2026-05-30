# Moiré Honeycomb 统一框架设计

## 概述

为 MoireSuperlattices 添加对**有效蜂窝晶格**的支持，并将 Wannier 函数、hopping 系数、库伦相互作用统一为适用于三角晶格和蜂窝晶格的框架。

## 动机

- 当前 `MoireTriangular` 只支持单原子 per unitcell 的有效三角晶格（MM 堆叠位点）
- 转角 TMD 超晶格中 MX 和 XM 堆叠位点构成有效的蜂窝晶格，需要双原子 per unitcell 的支持
- `MoireTriangularWannier` 只处理三角晶格，需要扩展为统一的 `MoireWannier`
- `CoulombIntegral`、hopping 计算需要一起统一

## 设计

### 1. MoireHoneycomb

参照 `MoireTriangular{N, D}` 的设计，新增 `MoireHoneycomb{N, D}`：

```julia
struct MoireHoneycomb{N, D<:Number} <: MoireSuperlattice{D}
    name::Symbol
    coordinates::Matrix{D}        # 2×2, 列分别对应 MX 和 XM 位置
    vectors::SVector{2, SVector{2, D}}
    neighbors::NTuple{N, Vector{SVector{2, D}}}
end
```

**坐标**：MM 堆叠为原点，MX 在 `(v₁+v₂)/3`，XM 在 `(2v₁-v₂)/3`（以超晶格基矢为单位）。

**neighbors**：从 kind=1 到 kind=N，收集逻辑与 `MoireTriangular` 一致——按距离分壳层，不区分 intracell/intercell。双原子带来的所有 bond（MX→MX, MX→XM, XM→MX, XM→XM）由 `bonds(lattice, truncation)` 算法自动处理。

### 2. MoireWannier（统一）

`MoireTriangularWannier` 重命名为统一的 `MoireWannier`，泛型约束从 `L<:MoireTriangular` 扩大为 `L<:MoireSuperlattice`：

```julia
struct MoireWannier{L<:MoireSuperlattice, G<:MoireReciprocalLattice, B<:BrillouinZone}
    aₘ::Float64
    lattice::L
    reciprocallattice::G
    brillouinzone::B
    eigenvalues::Matrix{Float64}    # (nband, Nₖ) — 原始能带能量
    bloch::Array{ComplexF64, 3}    # (nG×nlayer, nband, Nₖ) — 原始 Bloch 态 a_{Gν}^{kl}
    U::Array{ComplexF64, 3}        # (nband, nband, Nₖ) — gauge 变换矩阵
end
```

**bloch 是 gauge 变换之前的原始本征态**，直接从 continuum model 对角化得到。`U` 独立存储变换矩阵，两者互补：

- `(w)(r, sublattice)`：内部对 `bloch * U` 做傅里叶变换
- `HoppingIntegral`：使用 `U ε U†` 计算 hopping

#### 构造流程

| 步骤 | Triangular (nband=1) | Honeycomb (nband=2) |
|------|---------------------|---------------------|
| 取本征态 | 单个 band 的原始本征矢 | 两个相邻 bands 的原始本征矢 |
| SU(2) rotation | 不需要 | 对角化 layer projection，最大化 ψ̃₁ 在 bottom layer、ψ̃₂ 在 top layer |
| U(1) gauge fix | ψ̃_k(r_MM) 实正 | ψ̃₁k(r_XM) 实正，ψ̃₂k(r_MX) 实正 |
| 组装 U(k) | `[e^{iφ_k}]` | `Ũ(k) × diag(e^{iφ₁k}, e^{iφ₂k})` |

#### 调用接口

```julia
(w::MoireWannier)(r, sublattice::Int) -> SVector{nlayer, ComplexF64}
```

- Triangular：`sublattice = 1`
- Honeycomb：`sublattice = 1` → W_XM，`sublattice = 2` → W_MX

#### 构造函数签名

```julia
# Triangular
MoireWannier(system::MoireSystem, lattice::MoireTriangular, bz::BrillouinZone; band::Int)

# Honeycomb
MoireWannier(system::MoireSystem, lattice::MoireHoneycomb, bz::BrillouinZone; bands::UnitRange{Int})
```

### 3. HoppingIntegral（新增）

与 `CoulombIntegral` 对称的设计模式，封装 hopping 振幅 \( t_{im,jn}(R) \) 的计算：

```julia
struct HoppingIntegral{W<:MoireWannier}
    wannier::W
end
```

**callable 接口**：

```julia
(h::HoppingIntegral)(R::AbstractVector) -> SMatrix{nband, nband, ComplexF64}
```

公式（来自 PRR 2, 033087 补充材料）：

\[
t(R) = \frac{1}{N} \sum_k e^{-ik\cdot R}\, U(k)\,\text{diag}(\varepsilon_{1k}, \varepsilon_{2k})\,U^\dagger(k)
\]

其中 `U(k)` 和 `ε(k)` 取自 `wannier.U` 和 `wannier.eigenvalues`。

- Triangular (nband=1)：返回 1×1 矩阵
- Honeycomb (nband=2)：返回 2×2 矩阵，sublattice 索引对应 MX/XM

### 4. CoulombIntegral（更新）

泛型约束从 `W<:MoireTriangularWannier` 更新为 `W<:MoireWannier`，内部维度处理适配 `nband`。

### 5. terms 函数

Hopping 和 Coulomb 拆分为两个函数：

```julia
terms(h::HoppingIntegral; order::Int=truncation(h.wannier.lattice)) -> NTuple{...}, Term}
terms(c::CoulombIntegral; order::Int=truncation(c.wannier.lattice)) -> NTuple{...}, Term}
```

`order` 控制截断壳层数。对称性等价逻辑（C₃ 旋转）：
- AA 和 BB 组各自闭合
- AB 与 BA 在 C₃ 下混合为一组

详细 Term 生成算法（含 amplitude 函数）在实现阶段进一步细化。

### 6. 系数函数 coefficients

**删除**。用户通过 `HoppingIntegral` callable 直接获取 hopping：

```julia
h = HoppingIntegral(wannier)
t_R = h(R)          # 任意 R
t_0 = h(zero(R))    # onsite
```

### 7. 文件/模块结构

从单一文件 `src/MoireSuperlattices.jl`（约 650 行）拆分为 3 个文件：

```
src/
├── MoireSuperlattices.jl    # 主模块，include 子文件
├── lattices.jl              # 晶格类型 (~230行)
│   CommensurateBilayerHoneycomb, MoireReciprocalLattice,
│   MoireSuperlattice, MoireTriangular, MoireHoneycomb, RealZone
├── systems.jl               # 内部自由度 + 连续模型 (~200行)
│   MoireSpinor, MoireSpace, MoireSystem, BLTMD, bltmd!, bltmdmap
└── analysis.jl              # Wannier + Integrals + terms (~220行)
    MoireWannier, HoppingIntegral, CoulombIntegral,
    BareCoulomb, ImageCoulomb, TanhCoulomb, terms
```

### 8. 向后兼容

- `MoireTriangularWannier` → `MoireWannier{MoireTriangular}`（通过 const 别名兼容或直接重命名）
- `terms(bltmd, lattice, bz)` → 替换为 `terms(HoppingIntegral(wannier))`
- `coefficients` → 删除，迁移到 `HoppingIntegral` callable

## 参考文献

- Pan, Wu, Das Sarma, *Band topology, Hubbard model, Heisenberg model, and Dzyaloshinskii-Moriya interaction in twisted bilayer WSe₂*, Phys. Rev. Research 2, 033087 (2020)
- Zhou, Dong, Gu, Li, *Itinerant topological magnons and spin excitons in twisted transition metal dichalcogenides: Mapping electron topology to spin counterpart*, Supplemental Material
