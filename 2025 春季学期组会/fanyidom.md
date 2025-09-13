<!-- Meanless: Chapter 16-->

# 离散坐标法 $\left( {{S}_{N}\text{-近似}}\right)$

### 16.1 引言

与球谐函数法类似，离散坐标法是一种将辐射传输方程（对于灰体介质，或在光谱基础上）转化为一组联立偏微分方程的工具。与 ${P}_{N}$ 法相似，离散坐标法或 ${S}_{N}$ 法可以达到任意阶数和精度，尽管高阶 ${S}_{N}$ 格式的数学公式要简单得多。${S}_{N}$ 法最初由Chandrasekhar [1] 在其关于恒星和大气辐射的工作中提出，但在传热学界起初并未受到太多关注。同样，与 ${P}_{N}$ 法一样，离散坐标法首先被系统地应用于中子输运理论问题，特别是Lee [2] 和Lathrop [3,4] 的工作。早期有一些未经优化的尝试，将该方法应用于一维平面热辐射问题（Love等人 [5,6]，Hottel等人 [7]，Roux和Smith [8,9]）。然而，直到过去40年，离散坐标法才被应用于并优化于一般的辐射传热问题，这主要归功于Fiveland [10-13] 和Truelove [14-16] 的开创性工作。

离散坐标法（DOM）基于辐射强度的方向变化的离散表示。通过求解预先选定的一组横跨 ${4\pi }$ 总立体角范围的离散方向上的辐射传输方程来找到输运问题的解。因此，离散坐标法仅仅是对辐射传输方程方向依赖性的有限差分。对立体角的积分通过数值求积来近似（例如，用于评估辐射源项、辐射热流等）。

如今，许多数值传热模型使用有限体积法而非有限差分法。类似地，人们也可以使用有限立体角进行方向（或角度）离散化。这种离散坐标法的变体在历史上被称为（辐射传输的）有限体积法，并且越来越受欢迎。在本书中，为了与前面提到的使用有限立体角的理念保持一致（并避免与空间离散的有限体积法混淆），它将被称作有限角法（FAM）。由于高阶实现方案的公式相对直接，DOM及其有限角法的近亲FAM受到了极大的关注，如今与 ${P}_{1}$ 近似一起，可能是最受欢迎的RTE求解器。它们的一些版本被集成在大多数商业计算流体动力学（CFD）代码中。Charest等人 [17] 和Coelho [18] 对DOM和FAM的能力和缺点进行了详细的综述。后者为通用几何形状提供了该方法最完整的描述，远远超出了我们在本书中能提供的细节。

在本章中，我们将首先为标准DOM建立偏微分方程组及其边界条件。随后的一节将描述如何将该方法应用于一维平面平行介质，另一节则处理球形和圆柱形几何。然后，将概述其在一般多维问题中的应用，并针对二维笛卡尔几何给出具体的简化和示例。接下来是FAM的推导和演示。最后，本章将以对其他相关方法的简要介绍结束。

### 16.2 一般关系

对于吸收、发射和各向异性散射介质的一般辐射传输方程，根据方程(9.21)，为

$$
\frac{dI}{ds} = \widehat{\mathbf{s}} \cdot  \nabla I\left( {\mathbf{r},\widehat{\mathbf{s}}}\right)  = \kappa \left( \mathbf{r}\right) {I}_{b}\left( \mathbf{r}\right)  - \beta \left( \mathbf{r}\right) I\left( {\mathbf{r},\widehat{\mathbf{s}}}\right)  + \frac{{\sigma }_{s}\left( \mathbf{r}\right) }{4\pi }{\int }_{4\pi }I\left( {\mathbf{r},{\widehat{\mathbf{s}}}^{\prime }}\right) \Phi \left( {\mathbf{r},{\widehat{\mathbf{s}}}^{\prime },\widehat{\mathbf{s}}}\right) d{\Omega }^{\prime }. \tag{16.1}
$$

<!-- Meanless: Radiative Heat Transfer. https://doi.org/10.1016/B978-0-12-818143-0.00024-9 563 Copyright © 2022 Elsevier Inc. All rights reserved.-->




<!-- Meanless: 564 Radiative Heat Transfer-->

方程(16.1)对灰体介质有效，或在光谱基础上对非灰体介质有效，并受边界条件

$$
I\left( {{\mathbf{r}}_{w},\widehat{\mathbf{s}}}\right)  = \epsilon \left( {\mathbf{r}}_{w}\right) {I}_{b}\left( {\mathbf{r}}_{w}\right)  + \frac{\rho \left( {\mathbf{r}}_{w}\right) }{\pi }{\int }_{\widehat{\mathbf{n}} \cdot  {\widehat{\mathbf{s}}}^{\prime } < 0}I\left( {{\mathbf{r}}_{w},{\widehat{\mathbf{s}}}^{\prime }}\right) \left| {\widehat{\mathbf{n}} \cdot  {\widehat{\mathbf{s}}}^{\prime }}\right| d{\Omega }^{\prime }, \tag{16.2}
$$

的约束，这里我们仅限于具有不透明、漫射发射和漫反射壁面的外壳。将方程(16.2)推广到更复杂的边界条件是直接的。

## 离散坐标方程

在离散坐标法中，方程(16.1)是为一组 $n$ 个不同的方向 ${\widehat{\mathbf{s}}}_{i},i = 1,2,\ldots ,n$ 求解的，并且方向上的积分由数值求积代替，即

$$
{\int }_{4\pi }f\left( \widehat{\mathbf{s}}\right) {d\Omega } \simeq  \mathop{\sum }\limits_{{i = 1}}^{n}{w}_{i}f\left( {\widehat{\mathbf{s}}}_{i}\right) , \tag{16.3}
$$

其中 ${w}_{i}$ 是与方向 ${\widehat{\mathbf{s}}}_{i}$ 相关联的求积权重。因此，方程(16.1)被一组 $n$ 个方程近似，

$$
{\widehat{\mathbf{s}}}_{i} \cdot  \nabla I\left( {\mathbf{r},{\widehat{\mathbf{s}}}_{i}}\right)  = \kappa \left( \mathbf{r}\right) {I}_{b}\left( \mathbf{r}\right)  - \beta \left( \mathbf{r}\right) I\left( {\mathbf{r},{\widehat{\mathbf{s}}}_{i}}\right)  + \frac{{\sigma }_{s}\left( \mathbf{r}\right) }{4\pi }\mathop{\sum }\limits_{{j = 1}}^{n}{w}_{j}I\left( {\mathbf{r},{\widehat{\mathbf{s}}}_{j}}\right) \Phi \left( {\mathbf{r},{\widehat{\mathbf{s}}}_{j},{\widehat{\mathbf{s}}}_{i}}\right) ,\;i = 1,2,\ldots ,n, \tag{16.4}
$$

受边界条件约束

$$
I\left( {{\mathbf{r}}_{w},{\widehat{\mathbf{s}}}_{i}}\right)  = \epsilon \left( {\mathbf{r}}_{w}\right) {I}_{b}\left( {\mathbf{r}}_{w}\right)  + \frac{\rho \left( {\mathbf{r}}_{w}\right) }{\pi }\mathop{\sum }\limits_{{\widehat{\mathbf{n}} \cdot  {\widehat{\mathbf{s}}}_{j} < 0}}{w}_{j}I\left( {{\mathbf{r}}_{w},{\widehat{\mathbf{s}}}_{j}}\right) \left| {\widehat{\mathbf{n}} \cdot  {\widehat{\mathbf{s}}}_{j}}\right| ,\;\widehat{\mathbf{n}} \cdot  {\widehat{\mathbf{s}}}_{i} > 0. \tag{16.5}
$$

每束沿 ${\widehat{\mathbf{s}}}_{i}$ 方向传播的光线与外壳表面相交两次：一次是光束从壁面射出时 $\left( {\widehat{\mathbf{n}} \cdot  {\widehat{\mathbf{s}}}_{i} > 0}\right)$ ，另一次是它撞击壁面被吸收或反射时 $\left( {\widehat{\mathbf{n}} \cdot  {\widehat{\mathbf{s}}}_{i} < 0}\right)$ 。控制方程是一阶的，只需要一个边界条件（对于射出强度，$\widehat{\mathbf{n}} \cdot  {\widehat{\mathbf{s}}}_{i} > 0$）。方程(16.4)及其边界条件(16.5)构成了一组关于未知量 ${I}_{i}\left( \mathbf{r}\right)  = I\left( {\mathbf{r},{\widehat{\mathbf{s}}}_{i}}\right)$ 的 $n$ 个联立一阶线性偏微分方程。${I}_{i}$ 的解可以使用任何标准技术（解析或数值）找到。如果存在散射 $\left( {{\sigma }_{s} \neq  0}\right)$ ，和/或边界壁是反射的，这些方程会以一种通常需要迭代程序的方式耦合。即使在没有散射和表面反射的情况下，如果存在辐射平衡，温度场也可能未知，必须从强度场计算得出，这同样需要迭代。只有在没有散射和壁面反射，并且温度场给定的情况下，强度 ${I}_{i}$ 的求解才是直接的（正如精确解一样）。

一旦确定了强度，就可以很容易地计算出所需的方向积分量。介质内部或表面的辐射热流可以从其定义，方程(9.55)中找到，

$$
\mathbf{q}\left( \mathbf{r}\right)  = {\int }_{4\pi }I\left( {\mathbf{r},\widehat{\mathbf{s}}}\right) \widehat{\mathbf{s}}{d\Omega } \simeq  \mathop{\sum }\limits_{{i = 1}}^{n}{w}_{i}{I}_{i}\left( \mathbf{r}\right) {\widehat{\mathbf{s}}}_{i}. \tag{16.6}
$$

入射辐射 $G$ [以及通过方程(9.62)得到的辐射热流散度]同样可以确定为

$$
G\left( \mathbf{r}\right)  = {\int }_{4\pi }I\left( {\mathbf{r},\widehat{\mathbf{s}}}\right) {d\Omega } \simeq  \mathop{\sum }\limits_{{i = 1}}^{n}{w}_{i}{I}_{i}\left( \mathbf{r}\right) . \tag{16.7}
$$

在表面上，热流也可以通过表面能量平衡[方程(4.1)和(3.16)]确定为

$$
\mathbf{q} \cdot  \widehat{\mathbf{n}}\left( {\mathbf{r}}_{w}\right)  = \epsilon \left( {\mathbf{r}}_{w}\right) \left\lbrack  {\pi {I}_{b}\left( {\mathbf{r}}_{w}\right)  - H\left( {\mathbf{r}}_{w}\right) }\right\rbrack   \simeq  \epsilon \left( {\mathbf{r}}_{w}\right) \left( {\pi {I}_{b}\left( {\mathbf{r}}_{w}\right)  - \mathop{\sum }\limits_{{\widehat{\mathbf{n}} \cdot  {\widehat{\mathbf{s}}}_{i} < 0}}{w}_{i}{I}_{i}\left( {\mathbf{r}}_{w}\right) \left| {\widehat{\mathbf{n}} \cdot  {\widehat{\mathbf{s}}}_{i}}\right| }\right) . \tag{16.8}
$$




<!-- Meanless: The Method of Discrete Ordinates ( ${S}_{N}$ -Approximation) Chapter-->

如果我们将分析限制在线性各向异性散射，即散射相函数为

$$
\Phi \left( {\mathbf{r},\widehat{\mathbf{s}},{\widehat{\mathbf{s}}}^{\prime }}\right)  = 1 + {A}_{1}\left( \mathbf{r}\right) {\widehat{\mathbf{s}}}^{\prime } \cdot  \widehat{\mathbf{s}}. \tag{16.9}
$$

的情况下，方程(16.4)和(16.5)可以写成更紧凑的形式。那么，使用方程(16.6)和(16.7)和/或方程(13.15)可得

$$
{\widehat{\mathbf{s}}}_{i} \cdot  \nabla {I}_{i} + \beta {I}_{i} = \kappa {I}_{b} + \frac{{\sigma }_{s}}{4\pi }\left( {G + {A}_{1}\mathbf{q} \cdot  {\widehat{\mathbf{s}}}_{i}}\right) ,\;i = 1,2,\ldots ,n, \tag{16.10}
$$

边界条件为

$$
{I}_{i} = \frac{{J}_{w}}{\pi } = {I}_{bw} - \frac{1 - \epsilon }{\epsilon \pi }\mathbf{q} \cdot  \widehat{\mathbf{n}},\;\widehat{\mathbf{n}} \cdot  {\widehat{\mathbf{s}}}_{i} > 0 \tag{16.11}
$$

在壳体表面。当然，辐射热流和入射辐射是待定未知数，需要从方程(16.6)和(16.7)的级数中根据方向强度确定。

## 离散坐标方向的选择

求积格式的选择是任意的，尽管为了保持对称性和满足某些条件，可能会对方向 ${\widehat{\mathbf{s}}}_{i}$ 和求积权重 ${w}_{i}$ 产生限制。通常选择完全对称（即，在任何 ${90}^{ \circ  }$ 旋转后不变的集合）并且满足零阶、一阶和二阶矩的方向和权重集合，即

$$
{\int }_{4\pi }{d\Omega } = {4\pi } = \mathop{\sum }\limits_{{i = 1}}^{n}{w}_{i} \tag{16.12a}
$$

$$
{\int }_{4\pi }\widehat{\mathbf{s}}{d\Omega } = \mathbf{0} = \mathop{\sum }\limits_{{i = 1}}^{n}{w}_{i}{\widehat{\mathbf{s}}}_{i} \tag{16.12b}
$$

$$
{\int }_{4\pi }\widehat{\mathbf{s}}\widehat{\mathbf{s}}{d\Omega } = \frac{4\pi }{3}\delta  = \mathop{\sum }\limits_{{i = 1}}^{n}{w}_{i}{\widehat{\mathbf{s}}}_{i}{\widehat{\mathbf{s}}}_{i}, \tag{16.12c}
$$

其中 $\delta$ 是单位张量[参见方程(15.30)]。满足所有这些标准的不同方向和权重集已被列表，例如，由Lee [2]和Lathrop和Carlson [19]提供。Fiveland [12]和Truelove [15]观察到，不同的坐标集可能导致精度有相当大的差异。他们注意到(i)强度在壁面处可能存在方向不连续性，以及(ii)壁面处重要的辐射热流是通过强度在 ${2\pi }$ 半空间上的一阶矩来评估的[方程(16.8)]。他们得出结论，坐标和权重集也应满足半空间上的一阶矩，即

$$
{\int }_{\widehat{\mathbf{n}} \cdot  \widehat{\mathbf{s}} < 0}\left| {\widehat{\mathbf{n}} \cdot  \widehat{\mathbf{s}}}\right| {d\Omega } = {\int }_{\widehat{\mathbf{n}} \cdot  \widehat{\mathbf{s}} > 0}\widehat{\mathbf{n}} \cdot  \widehat{\mathbf{s}}{d\Omega } = \pi  = \mathop{\sum }\limits_{{\widehat{\mathbf{n}} \cdot  {\widehat{\mathbf{s}}}_{i} > 0}}{w}_{i}\widehat{\mathbf{n}} \cdot  {\widehat{\mathbf{s}}}_{i}. \tag{16.13}
$$

虽然对于任意方向的表面法线不可能满足方程(16.13)，但如果 $\widehat{\mathbf{n}} = \widehat{\mathbf{i}},\widehat{\mathbf{j}}$ 或 $\widehat{\mathbf{k}}$ ，则可以为主方向满足它。Lathrop和Carlson [19]给出了满足(i)对称性要求，(ii)矩方程(16.12)，和(iii)半矩方程(16.13)（对于 $\widehat{\mathbf{n}}$ 的三个主方向）$^{1}$ 的坐标和权重集。前四组标记为 ${S}_{2} - ,{S}_{4} - ,{S}_{6} -$ 和 ${S}_{8}$ -近似的集合在表16.1中再现。在表中，${\xi }_{i},{\eta }_{i}$ 和 ${\mu }_{i}$ 是 ${\widehat{\mathbf{s}}}_{i}$ 的方向余弦，或

$$
{\widehat{\mathbf{s}}}_{i} = \left( {{\widehat{\mathbf{s}}}_{i} \cdot  \widehat{\mathbf{i}}}\right) \widehat{\mathbf{i}} + \left( {{\widehat{\mathbf{s}}}_{i} \cdot  \widehat{\mathbf{j}}}\right) \widehat{\mathbf{j}} + \left( {{\widehat{\mathbf{s}}}_{i} \cdot  \widehat{\mathbf{k}}}\right) \widehat{\mathbf{k}} = {\xi }_{i}\widehat{\mathbf{i}} + {\eta }_{i}\widehat{\mathbf{j}} + {\mu }_{i}\widehat{\mathbf{k}}. \tag{16.14}
$$

表16.1中只给出了正方向余弦，覆盖了总立体角范围4π的八分之一。为了覆盖整个 ${4\pi }$，${\xi }_{i},{\eta }_{i}$ 和 ${\mu }_{i}$ 的任何或所有值都可以是正的或负的。因此，每一行坐标包含八个不同的方向。例如，对于 ${S}_{2}$ -近似，不同的方向是 ${\widehat{\mathbf{s}}}_{1} = {0.577350}\left( {\widehat{\mathbf{i}} + \widehat{\mathbf{j}} + \widehat{\mathbf{k}}}\right) ,{\widehat{\mathbf{s}}}_{2} = {0.577350}\left( {\widehat{\mathbf{i}} + \widehat{\mathbf{j}} - \widehat{\mathbf{k}}}\right) ,\ldots ,{\widehat{\mathbf{s}}}_{8} =  - {0.577350}\left( {\widehat{\mathbf{i}} + \widehat{\mathbf{j}} + \widehat{\mathbf{k}}}\right)$。由于对称的 ${S}_{2}$ -近似不满足半矩条件，表16.1中还包括了一个非对称的 ${S}_{2}$ -近似，由Truelove [15]提出。这个近似满足两个主方向的方程(16.13)，应该应用于一维和二维问题，非对称项从中消失（如在下一节的例16.1中所示）。“${S}_{N}$ -近似”这个名称表示每个主方向使用 $N$ 个不同的方向余弦。例如，对于 ${S}_{4}$ -近似，${\xi }_{i} =  \pm  {0.295876}$ 和 $\pm  {0.908248}$（或 ${\eta }_{i}$ 或 ${\mu }_{i}$）。总共总是有 $n = N\left( {N + 2}\right)$ 个不同的方向需要考虑（由于对称性，对于一维和二维问题，其中许多可能是不必要的）。文献中可以找到其他几种求积方案。Carlson [20]提出了一组具有相等权重 ${w}_{i}$ 的方案（如表16.1中的 ${S}_{2}$ 和 ${S}_{4}$ 集）。Fiveland [21]给出了另外两个求积法和对所有离散坐标集适用性的良好回顾。其他记录求积集生成程序的出版物有Sánchez和Smith [22]以及El-Wakil和Sacadura [23]的。Thurgood及其同事[24]给出了一系列新的求积集，像 ${S}_{n}$ 集一样在 ${90}^{ \circ  }$ 旋转中对称，但方向排列不同，作者将其称为 ${T}_{n}$ 集。这些集总是产生正权重，并声称可以减少所谓的“射线效应”（稍后将在第588页讨论）。这些集已由$\mathrm{{Li}}$及其同事[25]进一步完善。Koch和Becker [26]对方向求积方案进行了全面的回顾，包括对其准确性的评估。以上坐标集都不能准确处理准直（即单向）辐照。为了解决这个问题，Li及其同事[27]开发了ISW方案，在常规求积集中增加了一个“无限小权重”的单坐标。

---

<!-- Footnote -->

1. 对称的 ${S}_{2}$ 近似除外。

<!-- Footnote -->

---




<!-- Meanless: 566 Radiative Heat Transfer-->

<!-- Media -->

表16.1 ${S}_{N}$ -近似的离散坐标 $\left( {N = 2,4,6,8}\right)$，引自[19]。

<table><tr><td rowspan="2">近似阶数</td><td colspan="3">坐标</td><td>权重</td></tr><tr><td>$\xi$</td><td>$\eta$</td><td>$\mu$</td><td>$w$</td></tr><tr><td>${S}_{2}$ (对称)</td><td>0.5773503</td><td>0.5773503</td><td>0.5773503</td><td>1.5707963</td></tr><tr><td>${S}_{2}$ (非对称)</td><td>0.5000000</td><td>0.7071068</td><td>0.5000000</td><td>1.5707963</td></tr><tr><td rowspan="3">${S}_{4}$</td><td>0.2958759</td><td>0.2958759</td><td>0.9082483</td><td>0.5235987</td></tr><tr><td>0.2958759</td><td>0.9082483</td><td>0.2958759</td><td>0.5235987</td></tr><tr><td>0.9082483</td><td>0.2958759</td><td>0.2958759</td><td>0.5235987</td></tr><tr><td rowspan="6">${S}_{6}$</td><td>0.1838670</td><td>0.1838670</td><td>0.9656013</td><td>0.1609517</td></tr><tr><td>0.1838670</td><td>0.6950514</td><td>0.6950514</td><td>0.3626469</td></tr><tr><td>0.1838670</td><td>0.9656013</td><td>0.1838670</td><td>0.1609517</td></tr><tr><td>0.6950514</td><td>0.1838670</td><td>0.6950514</td><td>0.3626469</td></tr><tr><td>0.6950514</td><td>0.6950514</td><td>0.1838670</td><td>0.3626469</td></tr><tr><td>0.9656013</td><td>0.1838670</td><td>0.1838670</td><td>0.1609517</td></tr><tr><td rowspan="10">${S}_{8}$</td><td>0.1422555</td><td>0.1422555</td><td>0.9795543</td><td>0.1712359</td></tr><tr><td>0.1422555</td><td>0.5773503</td><td>0.8040087</td><td>0.0992284</td></tr><tr><td>0.1422555</td><td>0.8040087</td><td>0.5773503</td><td>0.0992284</td></tr><tr><td>0.1422555</td><td>0.9795543</td><td>0.1422555</td><td>0.1712359</td></tr><tr><td>0.5773503</td><td>0.1422555</td><td>0.8040087</td><td>0.0992284</td></tr><tr><td>0.5773503</td><td>0.5773503</td><td>0.5773503</td><td>0.4617179</td></tr><tr><td>0.5773503</td><td>0.8040087</td><td>0.1422555</td><td>0.0992284</td></tr><tr><td>0.8040087</td><td>0.1422555</td><td>0.5773503</td><td>0.0992284</td></tr><tr><td>0.8040087</td><td>0.5773503</td><td>0.1422555</td><td>0.0992284</td></tr><tr><td>0.9795543</td><td>0.1422555</td><td>0.1422555</td><td>0.1712359</td></tr></table>

<!-- Media -->

### 16.3 一维平板

我们将首先演示如何将 ${S}_{N}$ 离散坐标法应用于一个简单的一维平面平行板，该板由两个漫射发射和反射的等温板界定。与前面的章节一样，我们将自己限制在线性各向异性散射，尽管推广到任意各向异性散射是直接的。我们在这里避免这样做，是为了使推导步骤更容易理解。如果我们选择 $z$ 作为两板之间的空间坐标 $\left( {0 \leq  z \leq  L}\right)$ ，并引入光学坐标 $\tau$ ，其中 ${d\tau } = {\beta dz}\left( {0 \leq  \tau  \leq  {\tau }_{L}}\right)$ ，则方程(16.4)变换为


<!-- Meanless: The Method of Discrete Ordinates (S ${}_{N}$ -Approximation) Chapter | :-->

$$
{\mu }_{i}\frac{d{I}_{i}}{d\tau } = \left( {1 - \omega }\right) {I}_{b} - {I}_{i} + \frac{\omega }{4\pi }\mathop{\sum }\limits_{{j = 1}}^{n}{w}_{j}{I}_{j}\left\lbrack  {1 + {A}_{1}\left( {{\mu }_{i}{\mu }_{j} + {\xi }_{i}{\xi }_{j} + {\eta }_{i}{\eta }_{j}}\right) }\right\rbrack  ,\;i = 1,2,\ldots ,n. \tag{16.15}
$$

<!-- Media -->

表16.2 一维 ${S}_{N}$ -近似的离散坐标 $\left( {N = 2,4,6,8}\right)$。

<table><tr><td rowspan="2">近似阶数</td><td>坐标</td><td>权重</td></tr><tr><td>$\mu$</td><td>${w}^{\prime }$</td></tr><tr><td>${S}_{2}$ (对称)</td><td>0.5773503</td><td>6.2831853</td></tr><tr><td>${S}_{2}$ (非对称)</td><td>0.5000000</td><td>6.2831853</td></tr><tr><td rowspan="2">${S}_{4}$</td><td>0.2958759</td><td>4.1887902</td></tr><tr><td>0.9082483</td><td>2.0943951</td></tr><tr><td rowspan="3">${S}_{6}$</td><td>0.1838670</td><td>2.7382012</td></tr><tr><td>0.6950514</td><td>2.9011752</td></tr><tr><td>0.9656013</td><td>0.6438068</td></tr><tr><td rowspan="4">${S}_{8}$</td><td>0.1422555</td><td>2.1637144</td></tr><tr><td>0.5773503</td><td>2.6406988</td></tr><tr><td>0.8040087</td><td>0.7938272</td></tr><tr><td>0.9795543</td><td>0.6849436</td></tr></table>

<!-- Media -->

对于一维平板，强度与方位角无关。由于对于每个坐标 $j$（具有给定的 ${\mu }_{j}$）且 ${\xi }_{j}$ 为正值，都有另一个具有相同但为负值的坐标，并且两种坐标的强度相同，因此方程(16.15)中涉及 ${\xi }_{j}$ 的项相加为零。涉及 ${\eta }_{j}$ 的项也是如此，但涉及 ${\mu }_{j}$ 的项则不然（因为强度确实依赖于极角 $\theta$，且 $\mu  = \cos \theta$）。然而，涉及 ${\mu }_{j}$ 的项被重复了多次：表16.1中一行显示的每个 $\mu$ 值（分开计算正负 $\mu$ 值）对应四个不同的坐标（$\xi$ 和 $\eta$ 的正负值组合）。此外，一个特定的 $\mu$ 值可能出现在表16.1的多行中。如果将对应于单个 $\mu$ 值的所有求积权重相加，方程(16.15)简化为

$$
{\mu }_{i}\frac{d{I}_{i}}{d\tau } = \left( {1 - \omega }\right) {I}_{b} - {I}_{i} + \frac{\omega }{4\pi }\mathop{\sum }\limits_{{j = 1}}^{N}{w}_{j}^{\prime }{I}_{j}\left( {1 + {A}_{1}{\mu }_{i}{\mu }_{j}}\right) ,\;i = 1,2,\ldots ,N, \tag{16.16}
$$

其中 ${w}_{i}^{\prime }$ 是求和后的求积权重，求和到 $N$，即方法的阶数或不同方向余弦的数量。例如，对于 ${S}_{4}$ -近似中的 $\mu  = {0.2958759}$，求和后的求积权重是 ${w}^{\prime } = 4 \times  \left( {{0.5235987} + {0.5235987}}\right)  = {4\pi }/3$，依此类推。一维平板的坐标和求积权重列于表16.2。方程(16.16)本可以通过使用方程(16.10)而非(16.4)更不费力地找到，直接导致

$$
{\mu }_{i}\frac{d{I}_{i}}{d\tau } + {I}_{i} = \left( {1 - \omega }\right) {I}_{b} + \frac{\omega }{4\pi }\left( {G + {A}_{1}q{\mu }_{i}}\right) ,\;i = 1,2,\ldots ,N. \tag{16.17}
$$

在继续讨论方程(16.17)的边界条件之前，我们应该认识到，$N$ 个不同强度中，一半从 $\tau  = 0$ 处的壁面射出（${\mu }_{i} > 0$），另一半从 $\tau  = {\tau }_{L}$ 处的壁面射出（${\mu }_{i} < 0$）。遵循第13章的符号，我们将 $N$ 个不同的 ${I}_{i}$ 替换为

$$
{I}_{1}^{ + },{I}_{2}^{ + },\ldots ,{I}_{N/2}^{ + }\text{和}{I}_{1}^{ - },{I}_{2}^{ - },\ldots ,{I}_{N/2}^{ - }\text{.}
$$




<!-- Meanless: 568 Radiative Heat Transfer-->

然后方程(16.17)可以重写为

$$
{\mu }_{i}\frac{d{I}_{i}^{ + }}{d\tau } + {I}_{i}^{ + } = \left( {1 - \omega }\right) {I}_{b} + \frac{\omega }{4\pi }\left( {G + {A}_{1}q{\mu }_{i}}\right) , \tag{16.18a}
$$

$$
 - {\mu }_{i}\frac{d{I}_{i}^{ - }}{d\tau } + {I}_{i}^{ - } = \left( {1 - \omega }\right) {I}_{b} + \frac{\omega }{4\pi }\left( {G - {A}_{1}q{\mu }_{i}}\right) , \tag{16.18b}
$$

$$
i = 1,2,\ldots ,N/2;\;{\mu }_{i} > 0.
$$

使用这个符号，方程(16.18)的边界条件由方程(16.5)或(16.11)得出

$$
\tau  = 0 : \;{I}_{i}^{ + } = {J}_{1}/\pi  = {I}_{b1} - \frac{1 - {\epsilon }_{1}}{{\epsilon }_{1}\pi }{q}_{1}, \tag{16.19a}
$$

$$
\tau  = {\tau }_{L} : \;{I}_{i}^{ - } = {J}_{2}/\pi  = {I}_{b2} + \frac{1 - {\epsilon }_{2}}{{\epsilon }_{2}\pi }{q}_{2}, \tag{16.19b}
$$

$$
i = 1,2,\ldots ,N/2,\;{\mu }_{i} > 0.
$$

（对于在 ${\tau }_{L}$ 处的边界条件，符号会改变，因为 $\widehat{n}$ 指向与 $z$ 相反的方向。）

辐射热流 $q$ 和入射辐射 $G$ 通过方程(16.6)和(16.7)与方向强度相关，或

$$
q = \mathop{\sum }\limits_{{i = 1}}^{{N/2}}{w}_{i}^{\prime }{\mu }_{i}\left( {{I}_{i}^{ + } - {I}_{i}^{ - }}\right) , \tag{16.20a}
$$

$$
G = \mathop{\sum }\limits_{{i = 1}}^{{N/2}}{w}_{i}^{\prime }\left( {{I}_{i}^{ + } + {I}_{i}^{ - }}\right) . \tag{16.20b}
$$

在两个表面上，辐射热流更方便地从方程(16.8)评估为

$$
\tau  = 0 : \;{q}_{1} = q\left( 0\right)  = {\epsilon }_{1}\left( {{E}_{b1} - \mathop{\sum }\limits_{{i = 1}}^{{N/2}}{w}_{i}^{\prime }{\mu }_{i}{I}_{i}^{ - }}\right) , \tag{16.21a}
$$

$$
\tau  = {\tau }_{L} : \;{q}_{2} =  - q\left( {\tau }_{L}\right)  =  - {\epsilon }_{2}\left( {{E}_{b2} - \mathop{\sum }\limits_{{i = 1}}^{{N/2}}{w}_{i}^{\prime }{\mu }_{i}{I}_{i}^{ + }}\right) . \tag{16.21b}
$$

例16.1. 考虑两块大的、平行的、灰体-漫射且等温的板，相距为 $L$。一块板的温度为 ${T}_{1}$，发射率为 ${\epsilon }_{1}$，另一块板的温度为 ${T}_{2}$，发射率为 ${\epsilon }_{2}$。两板之间的介质是灰色的、吸收/发射和线性各向异性散射的气体 $\left( {n = 1}\right)$，具有恒定的消光系数 $\beta$ 和单次散射反照率 $\omega$。假设存在辐射平衡，使用 ${S}_{2}$ -近似确定两板之间的辐射热流。

## 解

对于辐射平衡，我们有，从方程(9.62)得到 ${I}_{b} = G/{4\pi }$ 且 $q =$ 常数；方程(16.18)和(16.19)变为

$$
{\mu }_{1}\frac{d{I}_{1}^{ + }}{d\tau } + {I}_{1}^{ + } = \frac{1}{4\pi }\left( {G + {A}_{1}\omega {\mu }_{1}q}\right) ,
$$

$$
 - {\mu }_{1}\frac{d{I}_{1}^{ - }}{d\tau } + {I}_{1}^{ - } = \frac{1}{4\pi }\left( {G - {A}_{1}\omega {\mu }_{1}q}\right) ,
$$

$$
\tau  = 0 : \;{I}_{1}^{ + } = {J}_{1}/\pi ,\tau  = {\beta L} = {\tau }_{L} : \;{I}_{1}^{ - } = {J}_{2}/\pi .
$$

对于 ${S}_{2}$ -近似，我们只有一个坐标方向 ${\mu }_{1}$（对于 ${I}_{1}^{ + }$ 指向 ${\tau }_{L}$，对于 ${I}_{1}^{ - }$ 指向 0），其中对于对称 ${S}_{2}$ -近似，${\mu }_{1} = {0.57735}$，对于非对称 ${S}_{2}$ -近似[满足半空间矩，方程(16.13)]，${\mu }_{1} = {0.5}$。对于简单的 ${S}_{2}$ -近似，联立方程（在这种情况下只有两个）可以被分离。我们在这里通过消除 ${I}_{1}^{ + }$ 和 ${I}_{1}^{ - }$，转而使用 $G$ 和 $q$ 来实现。从方程(16.20)，其中 ${w}_{i}^{\prime } = {2\pi }$，


<!-- Meanless: The Method of Discrete Ordinates ( ${S}_{N}$ -Approximation) Chapter $\mid$-->

$$
G = {2\pi }\left( {{I}_{1}^{ + } + {I}_{1}^{ - }}\right) ,
$$

$$
q = {2\pi }{\mu }_{1}\left( {{I}_{1}^{ + } - {I}_{1}^{ - }}\right) .
$$

因此，将两个微分方程相加和相减，并乘以 ${2\pi }$，得到

$$
\frac{dq}{d\tau } + G = G,\;\text{ 或 }\;\frac{dq}{d\tau } = 0,
$$

$$
{\mu }_{1}\frac{dG}{d\tau } + \frac{1}{{\mu }_{1}}q = {A}_{1}\omega {\mu }_{1}q,\;\text{ 或 }\;\frac{dG}{d\tau } =  - \left( {\frac{1}{{\mu }_{1}^{2}} - {A}_{1}\omega }\right) q.
$$

第一个方程仅仅是辐射平衡的重述，而第二个方程可以积分（因为 $q =$ 常数），或

$$
G = C - \left( {\frac{1}{{\mu }_{1}^{2}} - {A}_{1}\omega }\right) {q\tau }.
$$

这个关系包含两个未知常数（C和q），必须从边界条件确定，即

$$
\tau  = 0 : \;{I}_{1}^{ + } = \frac{1}{4\pi }\left( {G + \frac{q}{{\mu }_{1}}}\right)  = {J}_{1}/\pi ,
$$

$$
\tau  = {\tau }_{L} : \;{I}_{1}^{ - } = \frac{1}{4\pi }\left( {G - \frac{q}{{\mu }_{1}}}\right)  = {J}_{2}/\pi ,
$$

或

$$
\tau  = 0 : \;4{J}_{1} = G + \frac{q}{{\mu }_{1}} = C + \frac{q}{{\mu }_{1}},
$$

$$
\tau  = {\tau }_{L} : \;4{J}_{2} = G - \frac{q}{{\mu }_{1}} = C - \left( {\frac{1}{{\mu }_{1}^{2}} - {A}_{1}\omega }\right) q{\tau }_{L} - \frac{q}{{\mu }_{1}}.
$$

相减，我们得到，

$$
\Psi  = \frac{q}{{J}_{1} - {J}_{2}} = \frac{2{\mu }_{1}}{1 + \left( {1/{\mu }_{1}^{2} - {A}_{1}\omega }\right) {\mu }_{1}{\tau }_{L}/2},
$$

由此可以通过方程(13.48)消去辐射出射度。对于对称 ${S}_{2}$ -近似，${\mu }_{1} =$ ${0.57735} = 1/\sqrt{3}$，且对于各向同性散射，${A}_{1} = 0$，此表达式变为

$$
{\Psi }_{\text{对称 }} = \frac{1}{\sqrt{3}/2 + 3{\tau }_{L}/4}.
$$

另一方面，对于非对称 ${S}_{2}$ -近似 $\left( {{\mu }_{1} = {0.5}}\right)$，同样对于各向同性散射，

$$
{\Psi }_{\text{非对称 }} = \frac{1}{1 + {\tau }_{L}}.
$$

两种 ${S}_{2}$ -近似的结果在表16.3中与 ${P}_{1}$ -近似和精确解进行了比较。可以看出，${S}_{2}$ -方法的精度大约与 ${P}_{1}$ -近似相当。非对称 ${S}_{2}$ -近似优于对称近似，因为对称 ${S}_{2}$ 不满足半矩条件，方程(16.13)，并在光学薄极限下导致显著误差。

${S}_{2}$ -近似与第14.3节中讨论的双通量法相同，而非对称 ${S}_{2}$ -法不过是Schuster-Schwarzschild近似。

作为一维离散坐标法的第二个例子，我们将重复例15.4，该例最初旨在演示 ${P}_{3}$ -近似的使用。

例16.2. 考虑一个温度为 $T$ 的等温介质，被限制在两块大的、平行的、等温（温度相同为 ${T}_{w}$）的黑板之间。该介质是灰色的，吸收和发射，但不散射。使用 ${S}_{2}$ 和 ${S}_{4}$ 离散坐标近似确定介质内的传热率表达式。解


<!-- Meanless: 70 Radiative Heat Transfer-->

<!-- Media -->

表16.3 通过处于辐射平衡的一维平面平行介质的辐射热流；${S}_{2}$- 和 ${P}_{1}$-近似的比较。

<table><tr><td rowspan="2">${\tau }_{L}$</td><td colspan="4">$\Psi  = q/\left( {{J}_{1} - {J}_{2}}\right)$</td></tr><tr><td>精确解</td><td>${S}_{2}\left( {\mathbf{对称}}\right)$</td><td>${S}_{2}$ (非对称)</td><td>${P}_{1}$</td></tr><tr><td>0.0</td><td>1.0000</td><td>1.1547</td><td>1.0000</td><td>1.0000</td></tr><tr><td>0.1</td><td>0.9157</td><td>1.0627</td><td>0.9091</td><td>0.9302</td></tr><tr><td>0.5</td><td>0.7040</td><td>0.8058</td><td>0.6667</td><td>0.7273</td></tr><tr><td>1.0</td><td>0.5532</td><td>0.6188</td><td>0.5000</td><td>0.5714</td></tr><tr><td>5.0</td><td>0.2077</td><td>0.2166</td><td>0.1667</td><td>0.2105</td></tr></table>

<!-- Media -->

对于这个特别简单的情况，方程(16.18)简化为

$$
{\mu }_{i}\frac{d{I}_{i}^{ + }}{d\tau } + {I}_{i} = {I}_{b}
$$

$$
 - {\mu }_{i}\frac{d{I}_{i}^{ - }}{d\tau } + {I}_{i} = {I}_{b}
$$

由于 ${I}_{b} =$ 常数，这些方程可以立即积分，得到

$$
{I}_{i}^{ + } = {I}_{b} + {C}^{ + }{e}^{-\tau /{\mu }_{i}},
$$

$$
{I}_{i}^{ - } = {I}_{b} + {C}^{ - }{e}^{\tau /{\mu }_{i}}.
$$

积分常数 ${C}^{ + }$ 和 ${C}^{ - }$ 可以从边界条件(16.19)中找到

$$
\tau  = 0 : \;{I}_{i}^{ + } = {I}_{bw} = {I}_{b} + {C}^{ + },\;\text{ 或 }\;{C}^{ + } = {I}_{bw} - {I}_{b};
$$

$$
\tau  = {\tau }_{L} : \;{I}_{i}^{ - } = {I}_{bw} = {I}_{b} + {C}^{ - }{e}^{{\tau }_{L}/{\mu }_{i}},\;\text{ 或 }\;{C}^{ - } = \left( {{I}_{bw} - {I}_{b}}\right) {e}^{-{\tau }_{L}/{\mu }_{i}}.
$$

因此，

$$
{I}_{i}^{ + } = {I}_{b} + \left( {{I}_{bw} - {I}_{b}}\right) {e}^{-\tau /{\mu }_{i}},
$$

$$
{I}_{i}^{ - } = {I}_{b} + \left( {{I}_{bw} - {I}_{b}}\right) {e}^{-\left( {{\tau }_{L} - \tau }\right) /{\mu }_{i}}.
$$

辐射热流然后由方程(16.20)得到

$$
q = \mathop{\sum }\limits_{{i = 1}}^{{N/2}}{w}_{i}^{\prime }{\mu }_{i}\left( {{I}_{bw} - {I}_{b}}\right) \left( {{e}^{-\tau /{\mu }_{i}} - {e}^{-\left( {{\tau }_{L} - \tau }\right) /{\mu }_{i}}}\right) ,
$$

或者，以无量纲形式，

$$
\Psi  = \frac{q}{{n}^{2}\sigma \left( {{T}_{w}^{4} - {T}^{4}}\right) } = \frac{1}{\pi }\mathop{\sum }\limits_{{i = 1}}^{{N/2}}{w}_{i}^{\prime }{\mu }_{i}\left( {{e}^{-\tau /{\mu }_{i}} - {e}^{-\left( {{\tau }_{L} - \tau }\right) /{\mu }_{i}}}\right) .
$$

对于非对称 ${S}_{2}$ -近似，我们有 ${w}_{1}^{\prime } = {2\pi }$ 和 ${\mu }_{1} = {0.5}$，或

$$
{\Psi }_{{S}_{2}} = {e}^{-{2\tau }} - {e}^{-2\left( {{\tau }_{L} - \tau }\right) }.
$$

对于 ${S}_{4}$ -近似，${w}_{1}^{\prime } = {4\pi }/3,{w}_{2}^{\prime } = {2\pi }/3,{\mu }_{1} = {0.2958759},{\mu }_{2} = {0.9082483}$，且 $\sum {w}_{i}^{\prime }{\mu }_{i} = \pi$，所以

$$
{\Psi }_{{S}_{4}} = {0.3945012}\left( {{e}^{-\tau /{0.2958759}} - {e}^{-\left( {{\tau }_{L} - \tau }\right) /{0.2958759}}}\right)  + {0.6054088}\left( {{e}^{-\tau /{0.9082483}} - {e}^{-\left( {{\tau }_{L} - \tau }\right) /{0.9082483}}}\right) .
$$

结果应与例15.2和15.4中 ${P}_{1}$ -和 ${P}_{3}$ -近似的结果进行比较。注意，${S}_{N}$ -方法在壁面处趋于正确的光学厚极限 $\left( {{\tau }_{L} \rightarrow  \infty }\right)$，即 $\Psi  \rightarrow  1$[如果满足方程(16.13)的半矩]。另一方面，${P}_{N}$ -近似对于这个特定的例子，会高估光学厚极限。


<!-- Meanless: The Method of Discrete Ordinates ( ${S}_{N}$ -Approximation) Chapter-->

应该强调的是，最后一个例子——处理一个非散射的、等温的介质——特别适合于离散坐标法。不应期望对于一个普遍问题，${S}_{4}$ -法比 ${P}_{3}$ -近似更准确。

许多研究人员已经通过离散坐标法解决了更复杂的一维问题。Fiveland [12] 考虑了与本节所介绍的相同情况，但允许任意各向异性散射。他通过有限差分法求解方程组，注意到更高阶的 ${S}_{N}$ -法需要更小的数值步长 ${\Delta \tau }$，以便获得稳定的解。Kumar及其同事 [28] 不仅允许任意各向异性散射，还考虑了具有镜面反射的边界以及具有准直辐照的边界（如第18章所讨论）。Stamnes及其同事 [29,30] 研究了与Kumar及其同事相同的问题，但还允许变化的辐射性质和表面上的通用双向反射函数。他们使用线性代数的方法解耦了联立方程组，并找到了以特征值和特征向量表示的精确解析解。使用一维离散坐标模型作为解决更复杂问题的工具的其他例子可以在[31-40]中找到。

### 16.4 一维同心球体和圆柱体

对于同心球体和圆柱体，应用离散坐标法并利用一维问题中的对称性，比平面平行板要困难得多。原因是在穿过这类外壳的直视线上传播时，局部方向余弦会发生变化。

## 同心球体

考虑两个半径分别为 ${R}_{1}$ 和 ${R}_{2}$ 的同心球体。内球表面的发射率为 ${\epsilon }_{1}$，并保持在等温温度 ${T}_{1}$，而外球的温度为 ${T}_{2}$，发射率为 ${\epsilon }_{2}$。如果介质内的温度仅是半径的函数，则辐射传输方程由方程(13.85)给出，

$$
\mu \frac{\partial I}{\partial r} + \frac{1 - {\mu }^{2}}{r}\frac{\partial I}{\partial \mu } + {\beta I} = {\beta S}, \tag{16.22a}
$$

或者，

$$
\frac{\mu }{{r}^{2}}\frac{\partial }{\partial r}\left( {{r}^{2}I}\right)  + \frac{1}{r}\frac{\partial }{\partial \mu }\left\lbrack  {\left( {1 - {\mu }^{2}}\right) I}\right\rbrack   + {\beta I} = {\beta S}, \tag{16.22b}
$$

其中 $\mu  = \cos \theta$ 是极角的余弦，从径向方向测量（见图13.8）。$S$ 是辐射源函数，

$$
S\left( {r,\mu }\right)  = \left( {1 - \omega }\right) {I}_{b} + \frac{\omega }{2}{\int }_{-1}^{1}I\left( {r,{\mu }^{\prime }}\right) \Phi \left( {\mu ,{\mu }^{\prime }}\right) d{\mu }^{\prime }. \tag{16.23}
$$

额外的困难在于，方程(16.22)包含一个关于方向余弦 $\mu$ 的导数，这个导数在离散坐标法中需要被离散化。将 ${S}_{N}$ -法应用于方程(16.22)，我们得到

$$
\frac{{\mu }_{i}}{{r}^{2}}\frac{d}{dr}\left( {{r}^{2}{I}_{i}}\right)  + \frac{1}{r}{\left\{  \frac{\partial }{\partial \mu }\left\lbrack  \left( 1 - {\mu }^{2}\right) I\right\rbrack  \right\}  }_{\mu  = {\mu }_{i}} + \beta {I}_{i} = \beta {S}_{i},\;i = 1,2,\ldots ,N, \tag{16.24}
$$

其中 ${S}_{i}$ 很容易从方程(16.23)确定（除非介质各向异性散射，否则与坐标方向无关）。方程(16.24)仅应用于 $N$ 个主坐标，因为与平板类似，没有方位角依赖性。由于方向向量 $\mu$ 是离散的，它的导数必须用有限差分来近似。我们可以写成

$$
{\left\{  \frac{\partial }{\partial \mu }\left\lbrack  \left( 1 - {\mu }^{2}\right) I\right\rbrack  \right\}  }_{\mu  = {\mu }_{i}} \simeq  \frac{{\alpha }_{i + 1/2}{I}_{i + 1/2} - {\alpha }_{i - 1/2}{I}_{i - 1/2}}{{w}_{i}^{\prime }}, \tag{16.25}
$$




<!-- Meanless: 72 Radiative Heat Transfer-->

这是一个中心差分，其中 ${I}_{i \pm  1/2}$ 在两个坐标之间的边界处评估，如图16.1所示。由于任何两个连续的 ${\mu }_{i}$ 之间的差是非均匀的，几何系数 $\alpha$ 是非恒定的，需要被确定。$\alpha$ 的值仅取决于差分格式，因此与强度无关，可以通过检查一个特别简单的强度场来确定。例如，如果两个球体都处于相同的温度，那么 ${I}_{b1} = {I}_{b2} = {I}_{b} =$ 常数，并且 $I = {I}_{b} =$ 常数。这导致

$$
{\alpha }_{i + 1/2} - {\alpha }_{i - 1/2} = {w}_{i}^{\prime }{\left\lbrack  \frac{\partial }{\partial \mu }\left( 1 - {\mu }^{2}\right) \right\rbrack  }_{\mu  = {\mu }_{i}} =  - 2{w}_{i}^{\prime }{\mu }_{i},\;i = 1,2,\ldots ,N. \tag{16.26}
$$

<!-- Media -->

<!-- figureText: $N + 1/2$ $N - 1/2$ ${\mu }_{N - 1}$ $N - 3/2$ 5/2 ${\mu }_{2}$ ${\mu }_{N}$ 。 ${\theta }_{1}$ 3/2 ${\mu }_{1}$ 1/2 -->

<img src="https://cdn.noedgeai.com/bo_d141kgv7aajc73bva480_9.jpg?x=699&y=222&w=444&h=872&r=0"/>

图 16.1 一维问题的方向离散化和离散坐标值。

<!-- Media -->

如果可以确定 ${\alpha }_{1/p}$ 的值，这个表达式可以作为 ${\alpha }_{i + 1/p}$ 的递推公式。这个值是通过注意到 ${I}_{1{/}_{2}}$ 是在 $\mu  =  - 1$ 处评估的（图16.1），这里 $\left( {1 - {\mu }^{2}}\right) I = 0$，因此 ${\alpha }_{1{/}_{2}} = 0$。类似地，${I}_{N + 1/2}$ 是在 $\mu  =  + 1$ 处评估的，并且 ${\alpha }_{N + 1/2} = 0$。方程(16.25)和(16.26)的有限差分格式满足关系[4]

$$
{\int }_{-1}^{+1}\frac{\partial }{\partial \mu }\left\lbrack  {\left( {1 - {\mu }^{2}}\right) I}\right\rbrack  {d\mu } = {\left. \left( 1 - {\mu }^{2}\right) I\right| }_{-1}^{+1} = 0
$$

$$
 = \mathop{\sum }\limits_{{i = 1}}^{N}{w}_{i}^{\prime }{\left\{  \frac{\partial }{\partial \mu }\left\lbrack  \left( 1 - {\mu }^{2}\right) I\right\rbrack  \right\}  }_{\mu  = {\mu }_{i}} = \mathop{\sum }\limits_{{i = 1}}^{N}\left( {{\alpha }_{i + 1/2}{I}_{i + 1/2} - {\alpha }_{i - 1/2}{I}_{i - 1/2}}\right) 
$$

$$
 = {\alpha }_{3/2}{I}_{3/2} - {\alpha }_{1/2}{I}_{1/2} + {\alpha }_{5/2}{I}_{5/2} - {\alpha }_{3/2}{I}_{3/2} + \cdots  \cdot  {\alpha }_{N + 1/2}{I}_{N + 1/2} - {\alpha }_{N - 1/2}{I}_{N - 1/2}
$$

最后，节点边界处的强度 ${I}_{i \pm  {V}_{i}}$ 需要用节点中心值 ${I}_{i}$ 来表示。我们这里将使用简单的线性平均，即 ${I}_{i + 1/2} \simeq  \frac{1}{2}\left( {{I}_{i} + {I}_{i + 1}}\right)$。现在方程(16.24)可以重写为

$$
\frac{{\mu }_{i}}{{r}^{2}}\frac{d}{dr}\left( {{r}^{2}{I}_{i}}\right)  + \frac{{\alpha }_{i + 1/2}{I}_{i + 1} + \left( {{\alpha }_{i + 1/2} - {\alpha }_{i - 1/2}}\right) {I}_{i} - {\alpha }_{i - 1/2}{I}_{i - 1}}{{2r}{w}_{i}^{\prime }} + \beta {I}_{i} = \beta {S}_{i},
$$




<!-- Meanless: The Method of Discrete Ordinates (S ${}_{N}$ -Approximation) Chapter $\mid$-->

或者，进行微分并使用方程(16.26)，

$$
{\mu }_{i}\frac{d{I}_{i}}{dr} + \frac{{\mu }_{i}}{r}{I}_{i} + \frac{{\alpha }_{i + 1/2}{I}_{i + 1} - {\alpha }_{i - 1/2}{I}_{i - 1}}{{2r}{w}_{i}^{\prime }} + \beta {I}_{i} = \beta {S}_{i}, \tag{16.27a}
$$

$$
{\alpha }_{i + 1/2} = {\alpha }_{i - 1/2} - 2{w}_{i}^{\prime }{\mu }_{i},\;{\alpha }_{1/2} = {\alpha }_{N + 1/2} = 0,\;i = 1,2,\ldots ,N. \tag{16.27b}
$$

方程(16.27)构成了关于 $N$ 个未知强度 ${I}_{i}$ 的一组 $N$ 个联立微分方程，受边界条件[参见方程(16.19)]的约束

$$
r = {R}_{1} : {I}_{i} = {J}_{1}/\pi  = {I}_{b1} - \frac{1 - {\epsilon }_{1}}{{\epsilon }_{1}\pi }{q}_{1},\;i = \frac{N}{2} + 1,\frac{N}{2} + 2,\ldots ,N\left( {{\mu }_{i} > 0}\right) , \tag{16.28a}
$$

$$
r = {R}_{2} : {I}_{i} = {J}_{2}/\pi  = {I}_{b2} + \frac{1 - {\epsilon }_{2}}{{\epsilon }_{2}\pi }{q}_{2},\;i = 1,2,\ldots ,\frac{N}{2}\left( {{\mu }_{i} < 0}\right) . \tag{16.28b}
$$

对于一维平板，辐射热流和入射辐射的评估[参见方程(16.20)和(16.21)]由下式给出

$$
G\left( r\right)  = \mathop{\sum }\limits_{{i = 1}}^{N}{w}_{i}^{\prime }{I}_{i}\left( r\right)  \tag{16.29a}
$$

$$
q\left( r\right)  = \mathop{\sum }\limits_{{i = 1}}^{N}{w}_{i}^{\prime }{\mu }_{i}{I}_{i}\left( r\right)  \tag{16.29b}
$$

和

$$
q\left( {R}_{1}\right)  = {q}_{1} = {\epsilon }_{1}\left( {{E}_{b1} + \mathop{\sum }\limits_{\substack{{i = 1} \\  \left( {{\mu }_{i} < 0}\right)  }}^{{N/2}}{w}_{i}^{\prime }{\mu }_{i}{I}_{i}}\right) , \tag{16.29c}
$$

$$
 - q\left( {R}_{2}\right)  = {q}_{2} = {\epsilon }_{2}\left( {{E}_{b2} - \mathop{\sum }\limits_{\substack{{i = N/2 + 1} \\  \left( {{\mu }_{i} > 0}\right)  }}^{N}{w}_{i}^{\prime }{\mu }_{i}{I}_{i}}\right) . \tag{16.29d}
$$

例16.3. 考虑一个处于辐射平衡状态的非散射介质，该介质包含在两个等温的灰色球体之间。介质的吸收系数可以假定为灰色且恒定。使用 ${S}_{2}$ -近似确定两个同心球体之间的辐射热流。

## 解

从方程(16.27)我们发现，当 $N = 2$ 时，${\alpha }_{{1}_{2}} = {\alpha }_{{5}_{2}} = 0,{\alpha }_{{3}_{2}} =  - 2{w}_{1}^{\prime }{\mu }_{1} = 2{w}_{2}^{\prime }{\mu }_{2} = {4\pi \mu }$（因为 ${\mu }_{2} =  - {\mu }_{1} > 0$；我们保留 $\mu  = {\mu }_{2}$ 作为一个非数值，以便在对称和非对称 ${S}_{2}$ -近似之间进行比较）。对于处于辐射平衡的灰色非散射介质，我们有 $\beta  = \kappa$ 和 $\nabla  \cdot  \mathbf{q} = 0$，源函数由方程(9.64)和(16.39)给出，$S = {I}_{b} = G/{4\pi }$。

$$
i = 1 : \; - \mu \frac{d{I}_{1}}{d\tau } - \frac{\mu }{\tau }{I}_{1} + \frac{\mu }{\tau }{I}_{2} + {I}_{1} = \frac{G}{4\pi } = \frac{1}{2}\left( {{I}_{1} + {I}_{2}}\right) ,
$$

$$
 - \mu \frac{d{I}_{1}}{d\tau } - \left( {\frac{\mu }{\tau } - \frac{1}{2}}\right) \left( {{I}_{1} - {I}_{2}}\right)  = 0,
$$

$$
i = 2 : \;\mu \frac{d{I}_{2}}{d\tau } + \frac{\mu }{\tau }{I}_{2} - \frac{\mu }{\tau }{I}_{1} + {I}_{2} = \frac{1}{2}\left( {{I}_{1} + {I}_{2}}\right) ,
$$

$$
\mu \frac{d{I}_{2}}{d\tau } - \left( {\frac{\mu }{\tau } + \frac{1}{2}}\right) \left( {{I}_{1} - {I}_{2}}\right)  = 0.
$$

两个方程相加仅得到辐射平衡的重述（如例16.1），但将它们相减（并乘以 ${w}_{i}^{\prime } = {2\pi }$）得到

$$
\mu \frac{d}{d\tau }\left\lbrack  {{2\pi }\left( {{I}_{1} + {I}_{2}}\right) }\right\rbrack   + {2\pi }\left( {{I}_{1} - {I}_{2}}\right)  = 0,
$$

或


<!-- Meanless: 574 Radiative Heat Transfer-->

$$
\frac{dG}{d\tau } =  - \frac{q}{{\mu }^{2}} =  - \frac{{\tau }^{2}q}{{\mu }^{2}}\frac{1}{{\tau }^{2}}.
$$

由于对于同心球体之间处于辐射平衡的介质，$Q = {4\pi }{r}^{2}q =$ 常数，因此 ${\tau }^{2}q =$ 常数，入射辐射可以通过积分找到，

$$
G\left( \tau \right)  = \frac{{\tau }^{2}{q1}}{{\mu }^{2}}\frac{1}{\tau } + C,
$$

其中两个常数 $\left( {{\tau }^{2}q}\right)$ 和 $C$ 仍是未知的，必须由边界条件，方程(16.28)来确定：

$$
{I}_{2}\left( {\tau }_{1}\right)  = {J}_{1}/\pi ,\;{I}_{1}\left( {\tau }_{2}\right)  = {J}_{2}/\pi .
$$

使用 $q$ 和 $G$ 的定义，方程(16.29)，

$$
q = {2\pi \mu }\left( {{I}_{2} - {I}_{1}}\right) \;\text{ 和 }\;G = {2\pi }\left( {{I}_{2} + {I}_{1}}\right) ,
$$

或

$$
{I}_{1} = \frac{1}{4\pi }\left( {G - \frac{q}{\mu }}\right) ,{I}_{2} = \frac{1}{4\pi }\left( {G + \frac{q}{\mu }}\right) ,
$$

边界条件可以重写为关于 $q$ 和 $G$ 的形式

$$
\tau  = {\tau }_{1} : \;4{J}_{1} = G + \frac{{q}_{1}}{\mu } = \frac{{\tau }_{1}{q}_{1}}{{\mu }^{2}} + C + \frac{{q}_{1}}{\mu } = \frac{{\tau }^{2}q}{{\mu }^{2}}\left( {\frac{1}{{\tau }_{1}} + \frac{\mu }{{\tau }_{1}^{2}}}\right)  + C,
$$

$$
\tau  = {\tau }_{2} : \;4{J}_{2} = G - \frac{{q}_{2}}{\mu } = \frac{{\tau }_{2}{q}_{2}}{{\mu }^{2}} + C - \frac{{q}_{2}}{\mu } = \frac{{\tau }^{2}q}{{\mu }^{2}}\left( {\frac{1}{{\tau }_{2}} - \frac{\mu }{{\tau }_{2}^{2}}}\right)  + C.
$$

从第一个边界条件中减去第二个边界条件，我们得到

$$
\Psi  = \frac{{\tau }^{2}}{{\tau }_{1}^{2}}\frac{q}{{J}_{1} - {J}_{2}} = \frac{1}{\frac{1}{4\mu }\left( {1 + \frac{{\tau }_{1}^{2}}{{\tau }_{2}^{2}}}\right)  + \frac{{\tau }_{1}}{4{\mu }^{2}}\left( {1 - \frac{{\tau }_{1}}{{\tau }_{2}}}\right) }.
$$

对于对称 ${S}_{2}$ -近似，$\mu  = 1/\sqrt{3}$，此方程变为

$$
{\Psi }_{\text{对称 }} = \frac{1}{\frac{\sqrt{3}}{4}\left( {1 + \frac{{\tau }_{1}^{2}}{{\tau }_{2}^{2}}}\right)  + \frac{3{\tau }_{1}}{4}\left( {1 - \frac{{\tau }_{1}}{{\tau }_{2}}}\right) },
$$

对于非对称近似，$\mu  = {0.5}$，

$$
{\Psi }_{\text{非对称 }} = \frac{1}{\frac{1}{2}\left( {1 + \frac{{\tau }_{1}^{2}}{{\tau }_{2}^{2}}}\right)  + {\tau }_{1}\left( {1 - \frac{{\tau }_{1}}{{\tau }_{2}}}\right) }.
$$

${S}_{2}$ -近似的精度与 ${P}_{1}$ -近似非常相似，后者为

$$
{\Psi }_{{P}_{1}} = \frac{1}{\frac{1}{2}\left( {1 + \frac{{\tau }_{1}^{2}}{{\tau }_{2}^{2}}}\right)  + \frac{3{\tau }_{1}}{4}\left( {1 - \frac{{\tau }_{1}}{{\tau }_{2}}}\right) }.
$$

注意，该方法对于大的 ${\tau }_{1}$（大光学厚度）非常精确，但在光学薄条件下 $\left( {\kappa  \rightarrow  0}\right)$ 会失效，特别是对于小的半径比 ${R}_{1}/{R}_{2}$。在极限 $\left( {\kappa  \rightarrow  0,{R}_{1}/{R}_{2} \rightarrow  0}\right)$ 下，我们发现 ${\Psi }_{{P}_{1}} = {\Psi }_{{S}_{2}\text{,非对称 }} \rightarrow  2$，而正确的极限应趋于 ${\Psi }_{\text{精确 }} \rightarrow  1$。

Tsai及其同事[41]使用 ${S}_{8}$ 离散坐标法和Fiveland [12]的等权重坐标，报告了方程(16.27)的数值解，允许各向异性散射、可变属性和外部辐照。Jones和Bayazitoğlu [42,43]使用相同的方法来确定通过球壳的传导和辐射的综合效应。


<!-- Meanless: The Method of Discrete Ordinates (S ${}_{N}$ -Approximation) Chapter-->

## 同心圆柱体

两个同心圆柱体的分析遵循类似的思路。我们再次考虑一个被两个等温圆柱体所包含的吸收、发射和散射介质，其半径分别为 ${R}_{1}$（温度 ${T}_{1}$，漫射发射率 ${\epsilon }_{1}$）和 ${R}_{2}$（温度 ${T}_{2}$，发射率 ${\epsilon }_{2}$）。对于这种情况，辐射传输方程由方程(13.104)给出，

$$
\sin \theta \cos \psi \frac{\partial I}{\partial r} - \frac{\sin \theta \sin \psi }{r}\frac{\partial I}{\partial \psi } + {\beta I} = {\beta S}, \tag{16.30}
$$

其中极角 $\theta$ 是从 $z$ 轴测量的，方位角 $\psi$ 是从局部径向方向测量的（参见图13.9）。$S$ 是辐射源函数，已由方程(16.23)给出。引入方向余弦 $\xi  = \widehat{\mathbf{s}} \cdot  {\widehat{\mathbf{e}}}_{z} = \cos \theta ,\mu  = \widehat{\mathbf{s}} \cdot  {\widehat{\mathbf{e}}}_{r} = \sin \theta \cos \psi$ 和 $\eta  = \widehat{\mathbf{s}} \cdot  {\widehat{\mathbf{e}}}_{{\psi }_{c}} = \sin \theta \sin \psi$，我们可以将方程(16.30)重写为

$$
\frac{\mu }{r}\frac{\partial }{\partial r}\left( {rI}\right)  - \frac{1}{r}\frac{\partial }{\partial \psi }\left( {\eta I}\right)  + {\beta I} = {\beta S}. \tag{16.31}
$$

对于一维圆柱形介质，对称条件不像平板和球体那样直接。这里我们有

$$
I\left( {r,\theta ,\psi }\right)  = I\left( {r,\pi  - \theta ,\psi }\right)  = I\left( {r,\theta , - \psi }\right) . \tag{16.32}
$$

因此，对于正负值的 $\xi$ 以及正负值的 $\eta$，强度是相同的。因此，我们只需要考虑表16.1中 ${\xi }_{i}$ 和 ${\eta }_{i}$ 的正值，这导致 ${S}_{N}$ -近似有 ${N}_{c} = N\left( {N + 2}\right) /4$ 个不同的坐标，其求积权重为 ${w}_{i}^{\prime \prime } = 4{w}_{i}$。方程(16.31)可以写成离散坐标形式

$$
\frac{{\mu }_{i}}{r}\frac{d}{dr}\left( {r{I}_{i}}\right)  - \frac{1}{r}{\left\{  \frac{\partial }{\partial \psi }\left( \eta I\right) \right\}  }_{\psi  = {\psi }_{i}} + \beta {I}_{i} = \beta {S}_{i},\;i = 1,2,\ldots ,{N}_{c}. \tag{16.33}
$$

与同心球体情况一样，大括号中的项近似为

$$
{\left\{  \frac{\partial }{\partial \psi }\left( \eta I\right) \right\}  }_{\psi  = {\psi }_{i}} \simeq  \frac{{\alpha }_{i + 1/2}{I}_{i + 1/2} - {\alpha }_{i - 1/2}{I}_{i - 1/2}}{{w}_{i}^{\prime \prime }},\;i = 1,2,\ldots ,{N}_{i},{\xi }_{i}\text{ 固定. } \tag{16.34}
$$

在这个关系中，下标 $i + 1{/}_{2}$ 意味着“朝向下一个更高的 ${\psi }_{i}$ 值，保持 ${\xi }_{i}$ 不变”。${N}_{i}$ 的值取决于 ${\xi }_{i}$ 的值。例如，对于 ${S}_{4}$ -近似，我们从表16.1中得到，当 ${\xi }_{i} = {0.2958759}$ 时 ${N}_{i} = 4$（四个不同的 ${\mu }_{i}$ 值，两个正两个负），当 ${\xi }_{i} = {0.9082483}$ 时 ${N}_{i} = 2$。在同心圆柱体的情况下，$\alpha$ 的递推公式是通过在方程(16.31)中令 $I = S =$ 常数得到的

$$
{\alpha }_{i + 1/2} - {\alpha }_{i - 1/2} = {\left. {w}_{i}^{\prime \prime }\frac{\partial \eta }{\partial \psi }\right| }_{\psi  = {\psi }_{i}} = {w}_{i}^{\prime \prime }{\mu }_{i},\;i = 1,2,\ldots ,{N}_{i},{\xi }_{i}\text{ 固定. } \tag{16.35}
$$

再次，${\alpha }_{{1}_{2}} = 0$ 因为在该位置 ${\psi }_{{1}_{2}} = 0$，因此 $\eta  = 0$。类似地，${\alpha }_{{N}_{i} + {1}_{2}} = 0$ 因为 ${\psi }_{{N}_{i} + {1}_{2}} = \pi$ 且 $\eta  = 0$。最后，对半节点强度使用线性平均得到

$$
{\mu }_{i}\frac{d{I}_{i}}{dr} + \frac{{\mu }_{i}}{2r}{I}_{i} - \frac{{\alpha }_{i + 1/2}{I}_{i + 1} - {\alpha }_{i - 1/2}{I}_{i - 1}}{{2r}{w}_{i}^{\prime \prime }} + \beta {I}_{i} = \beta {S}_{i},\;i = 1,2,\ldots ,{N}_{c}, \tag{16.36a}
$$

$$
{\alpha }_{i + 1/2} = {\alpha }_{i - 1/2} + {w}_{i}^{\prime \prime }{\mu }_{i},\;{\alpha }_{1/2} = {\alpha }_{N + 1/2} = 0,\;i = 1,2,\ldots {N}_{i},{\xi }_{i}\text{ 固定. } \tag{16.36b}
$$

方程(16.36)是关于 ${N}_{c} = N\left( {N + 2}\right) /4$ 个未知方向强度 ${I}_{i}$ 的同心圆柱体方程组，并且等价于同心球体的方程组，方程(16.27)。圆柱体和球体的边界条件基本相同[方程(16.28)]，除了一些重新编号，入射强度和辐射热流的表达式也相同[方程(16.29)]，即


<!-- Meanless: 576 Radiative Heat Transfer-->

$$
r = {R}_{1} : \;{I}_{i} = \frac{{J}_{1}}{\pi } = {I}_{b1} - \frac{1 - {\epsilon }_{1}}{{\epsilon }_{1}\pi }{q}_{1},\;i = \frac{{N}_{c}}{2} + 1,\frac{{N}_{c}}{2} + 2,\ldots ,{N}_{c}\left( {{\mu }_{i} > 0}\right) , \tag{16.37a}
$$

$$
r = {R}_{2} : \;{I}_{i} = \frac{{J}_{2}}{\pi } = {I}_{b2} + \frac{1 - {\epsilon }_{2}}{{\epsilon }_{2}\pi }{q}_{2},\;i = 1,2,\ldots ,\frac{{N}_{c}}{2}\left( {{\mu }_{i} < 0}\right) , \tag{16.37b}
$$

$$
G\left( r\right)  = \mathop{\sum }\limits_{{i = 1}}^{{N}_{c}}{w}_{i}^{\prime \prime }{I}_{i}\left( r\right)  \tag{16.37c}
$$

和

$$
q\left( r\right)  = \mathop{\sum }\limits_{{i = 1}}^{{N}_{c}}{w}_{i}^{\prime \prime }{\mu }_{i}{I}_{i}\left( r\right)  \tag{16.37d}
$$

$$
q\left( {R}_{1}\right)  = {q}_{1} = {\epsilon }_{1}\left( {{E}_{b1} + \mathop{\sum }\limits_{\substack{{i = 1} \\  \left( {{\mu }_{i} < 0}\right)  }}^{{{N}_{c}/2}}{w}_{i}^{\prime \prime }{\mu }_{i}{I}_{i}}\right) , \tag{16.37e}
$$

$$
q\left( {R}_{2}\right)  = {q}_{2} = {\epsilon }_{2}\left( {{E}_{b2} - \mathop{\sum }\limits_{\substack{{i = {N}_{c}/2 + 1} \\  \left( {{\mu }_{i} > 0}\right)  }}^{{N}_{c}}{w}_{i}^{\prime \prime }{\mu }_{i}{I}_{i}}\right) . \tag{16.37f}
$$

Krishnaprakas [44] 的工作是一维介质中离散坐标法使用的一个例子，他考虑了在灰色、恒定属性介质中具有各种散射行为的传导和辐射的组合。

### 16.5 多维问题

虽然离散坐标法很容易推广到多维构型，但该方法会产生一组联立的一阶偏微分方程，通常必须用数值方法求解。与一维几何体一样，无论采用笛卡尔坐标系、圆柱坐标系还是球坐标系，辐射传输方程都略有不同。我们将首先描述笛卡尔坐标系下的方法，然后简要概述该方法在多维笛卡尔、圆柱和球几何体中的应用。

## 笛卡尔坐标系下的一般公式

对于笛卡尔坐标，方程(16.4)利用方程(16.14)变为

$$
{\xi }_{i}\frac{\partial {I}_{i}}{\partial x} + {\eta }_{i}\frac{\partial {I}_{i}}{\partial y} + {\mu }_{i}\frac{\partial {I}_{i}}{\partial z} + \beta {I}_{i} = \beta {S}_{i},\;i = 1,2,\ldots ,n, \tag{16.38}
$$

其中 ${S}_{i}$ 再次是辐射源函数的简写

$$
{S}_{i} = \left( {1 - \omega }\right) {I}_{b} + \frac{\omega }{4\pi }\mathop{\sum }\limits_{{j = 1}}^{n}{w}_{j}{\Phi }_{ij}{I}_{j},\;i = 1,2,\ldots ,n. \tag{16.39}
$$

方程(16.38)受每个表面上的边界条件方程(16.5)的约束。例如，对于一个平行于 $y - z$ 平面的表面，其法向量为 $\widehat{\mathbf{n}} = \widehat{\mathbf{t}}$ 且 $\widehat{\mathbf{n}} \cdot  {\widehat{\mathbf{s}}}_{j} = {\widehat{\mathbf{s}}}_{j} \cdot  \widehat{\mathbf{t}} = {\xi }_{j}$，我们对于所有 ${\xi }_{i} > 0$ 的 $i$（$n/2$ 个边界条件）有

$$
{I}_{i} = {J}_{w}/\pi  = {\epsilon }_{w}{I}_{bw} + \frac{1 - {\epsilon }_{w}}{\pi }\mathop{\sum }\limits_{\substack{{\xi }_{j} < 0 \\  \left( {n/2\text{ 个值 }}\right)  }}{w}_{j}{I}_{j}\left| {\xi }_{j}\right| . \tag{16.40}
$$




<!-- Meanless: The Method of Discrete Ordinates (S ${}_{N}$ -Approximation) Chapter $\mid$-->

虽然方程(16.38)的数值解可以通过标准的有限差分法找到，但对于一般的非正交（曲线）网格，使用后向差分来近似强度的导数是繁琐的，因为弯曲的网格线不一定遵循射线（强度）的方向。因此，更常见的是采用Carlson和Lathrop [4]的有限体积法，该方法被用于下面描述的一般公式。此外，大多数现代CFD代码都是基于有限体积的，这使得基于空间有限体积离散化的离散坐标公式易于与此类代码集成，用于通用传热计算。

<!-- Media -->

<!-- figureText: ${N}_{1} \bullet$ ${N}_{4}$ ${\widehat{\mathbf{n}}}_{3}$ ${}^{ \bullet  }{N}_{3}$ -->

<img src="https://cdn.noedgeai.com/bo_d141kgv7aajc73bva480_14.jpg?x=602&y=226&w=630&h=519&r=0"/>

图 16.2 一个通用的控制体积（或单元格），$P$，显示了编号为1到4的单元格面，编号为${N}_{1}$到${N}_{4}$的相邻单元格，以及表面法线。单元格中心用圆圈表示，而面中心用方块表示。还显示了一个参考的坐标方向。

<!-- Media -->

为了推导有限体积方程，我们首先将方程(16.4)重写如下：

$$
{\widehat{\mathbf{s}}}_{i} \cdot  \nabla {I}_{i} = \nabla  \cdot  \left( {{I}_{i}{\widehat{\mathbf{s}}}_{i}}\right)  =  - \beta {I}_{i} + \beta {S}_{i},\;i = 1,2,\ldots ,n, \tag{16.41}
$$

其中 ${I}_{i} = I\left( {\mathbf{r},{\widehat{\mathbf{s}}}_{i}}\right)$，而 ${S}_{i}$ 由方程(16.39)给出。接下来，我们将方程(16.41)在一个任意的控制体积（或单元格）上积分，其体积为 ${V}_{P}$，如图16.2所示。因此，

$$
{\int }_{{V}_{P}}\nabla  \cdot  \left( {{I}_{i}{\widehat{\mathbf{s}}}_{i}}\right) {dV} =  - {\int }_{{V}_{P}}\beta {I}_{i}{dV} + {\int }_{{V}_{P}}\beta {S}_{i}{dV},\;i = 1,2,\ldots ,n. \tag{16.42}
$$

任意量 $\phi$ 在体积 ${V}_{P}$ 上的体积平均定义为

$$
\bar{\phi } = \frac{1}{{V}_{P}}{\int }_{{V}_{P}}{\phi dV}
$$

使用这个定义，方程(16.42)可以写成

$$
{\int }_{{V}_{P}}\nabla  \cdot  \left( {{I}_{i}{\widehat{\mathbf{s}}}_{i}}\right) {dV} =  - \overline{\beta {I}_{i}}{V}_{P} + \overline{\beta {S}_{i}}{V}_{P},\;i = 1,2,\ldots ,n. \tag{16.43}
$$

在传统的有限体积公式中，假设任何量的体积平均值是该量在体积几何中心处的值。根据这个假设，$\overline{\beta {I}_{i}} \approx  {\left( \beta {I}_{i}\right) }_{P} \approx  {\beta }_{P}{I}_{i,P}$ 和 $\overline{\beta {S}_{i}} \approx  {\left( \beta {S}_{i}\right) }_{P} \approx  {\beta }_{P}{S}_{i,P}$。这些近似的准确性取决于被近似的函数（或量）的非线性和进行平均的体积大小。对于线性变化的函数，这个近似实际上是一个恒等式。如果体积很小，任何非线性函数都可以被认为是局部线性的，因此近似所产生的误差很小。换句话说，当网格细化时，这个近似所产生的误差会变得非常小。使用这些近似，方程(16.43)可以简化为


<!-- Meanless: 578 Radiative Heat Transfer-->

$$
{\int }_{{V}_{P}}\nabla  \cdot  \left( {{I}_{i}{\widehat{\mathbf{s}}}_{i}}\right) {dV} =  - {\beta }_{P}{I}_{i,P}{V}_{P} + {\beta }_{P}{S}_{i,P}{V}_{P},\;i = 1,2,\ldots ,n. \tag{16.44}
$$

接下来，我们将高斯散度定理应用于方程(16.44)的左侧，得到

$$
{\int }_{{A}_{P}}{I}_{i}\left( {{\widehat{\mathbf{s}}}_{i} \cdot  \widehat{\mathbf{n}}}\right) {dA} =  - {\beta }_{P}{I}_{i,P}{V}_{P} + {\beta }_{P}{S}_{i,P}{V}_{P},\;i = 1,2,\ldots ,n, \tag{16.45}
$$

其中 ${A}_{P}$ 是其质心由 $P$ 表示的单元格的边界面积。为了进一步简化方程(16.45)，我们假设所考虑的体积是一个任意的凸多面体（在二维中为凸多边形），它由平面（在二维中为线段）界定，每个平面都有唯一的表面法线。由此可得

$$
\mathop{\sum }\limits_{{f = 1}}^{{N}_{f,P}}{\int }_{{A}_{f}}{I}_{i}\left( {{\widehat{\mathbf{s}}}_{i} \cdot  {\widehat{\mathbf{n}}}_{f}}\right) {dA} =  - {\beta }_{P}{I}_{i,P}{V}_{P} + {\beta }_{P}{S}_{i,P}{V}_{P},\;i = 1,2,\ldots ,n, \tag{16.46}
$$

其中 ${N}_{f,P}$ 是所讨论单元格（单元格 $P$）的面数，求和是对单元格的所有面进行的。${A}_{f}$ 是面 $f$ 的面积，而 ${\widehat{\mathbf{n}}}_{f}$ 是其向外的表面法线。由于每个面的表面法线是唯一的，并且在面上的位置不发生变化，所以量 $\left( {{\widehat{\mathbf{s}}}_{i} \cdot  {\widehat{\mathbf{n}}}_{f}}\right)$ 对每个面都是一个常数。因此，方程(16.46)可以重新整理为

$$
\mathop{\sum }\limits_{{f = 1}}^{{N}_{f,P}}\left( {{\widehat{\mathbf{s}}}_{i} \cdot  {\widehat{\mathbf{n}}}_{f}}\right) {\int }_{{A}_{f}}{I}_{i}{dA} =  - {\beta }_{P}{I}_{i,P}{V}_{P} + {\beta }_{P}{S}_{i,P}{V}_{P},\;i = 1,2,\ldots ,n. \tag{16.47}
$$

将体积平均的定义扩展到面积平均，可以写成 ${\bar{I}}_{i}{A}_{f} = {\int }_{{A}_{f}}{I}_{i}{dA}$，其中 $\overline{{I}_{i}}$ 是面 $f$ 上的面积平均强度。再次，根据前面的论证，它可以近似为在面质心（或面中心；见图16.2）处评估的强度，即 $\overline{{I}_{i}} \approx  {I}_{i,f}$。将这个近似代入方程(16.47)得到

$$
\mathop{\sum }\limits_{{f = 1}}^{{N}_{f,P}}\left( {{\widehat{\mathbf{s}}}_{i} \cdot  {\widehat{\mathbf{n}}}_{f}}\right) {I}_{i,f}{A}_{f} =  - {\beta }_{P}{I}_{i,P}{V}_{P} + {\beta }_{P}{S}_{i,P}{V}_{P},\;i = 1,2,\ldots ,n. \tag{16.48}
$$

接下来的任务是用单元格中心强度来表示面强度 ${I}_{i,f}$。最常用的方法是使用迎风格式 [45]。由于强度是向前传播的，用射线来源的单元格的强度值来近似面强度是合乎逻辑的。在数学上，这可以写成

$$
{I}_{i,f} \approx  \left\{  {\begin{array}{ll} {I}_{i,P} & \text{ 如果 }\left( {{\widehat{\mathbf{s}}}_{i} \cdot  {\widehat{\mathbf{n}}}_{f}}\right)  \geq  0 \\  {I}_{i,{N}_{f}} & \text{ 如果 }\left( {{\widehat{\mathbf{s}}}_{i} \cdot  {\widehat{\mathbf{n}}}_{f}}\right)  < 0 \end{array},}\right.  \tag{16.49}
$$

其中下标 ${N}_{f}$ 表示面 $f$ 的相邻单元格的单元格中心强度（或相邻单元格的 ${I}_{i,P}$），如图16.2所示。与其使用方程(16.49)中所示的逻辑“if”检查来确定面强度，不如使用所谓的切换函数 [45] 来用单元格中心强度表示面强度，这样更方便：

$$
\left( {{\widehat{\mathbf{s}}}_{i} \cdot  {\widehat{\mathbf{n}}}_{f}}\right) {I}_{i,f} \approx  \frac{\left| \left( {{\widehat{\mathbf{s}}}_{i} \cdot  {\widehat{\mathbf{n}}}_{f}}\right) \right|  + \left( {{\widehat{\mathbf{s}}}_{i} \cdot  {\widehat{\mathbf{n}}}_{f}}\right) }{2}{I}_{i,P} - \frac{\left| \left( {{\widehat{\mathbf{s}}}_{i} \cdot  {\widehat{\mathbf{n}}}_{f}}\right) \right|  - \left( {{\widehat{\mathbf{s}}}_{i} \cdot  {\widehat{\mathbf{n}}}_{f}}\right) }{2}{I}_{i,{N}_{f}} \tag{16.50}
$$

将方程(16.50)代入方程(16.48)，并重新整理所得方程，我们得到

$$
\left\lbrack  {\mathop{\sum }\limits_{{f = 1}}^{{N}_{f,P}}\left( {\frac{\left| \left( {{\widehat{\mathbf{s}}}_{i} \cdot  {\widehat{\mathbf{n}}}_{f}}\right) \right|  + \left( {{\widehat{\mathbf{s}}}_{i} \cdot  {\widehat{\mathbf{n}}}_{f}}\right) }{2}{A}_{f}}\right)  + {\beta }_{P}{V}_{P}}\right\rbrack  {I}_{i,P} - \left\lbrack  {\mathop{\sum }\limits_{{f = 1}}^{{N}_{f,P}}\left( {\frac{\left| \left( {{\widehat{\mathbf{s}}}_{i} \cdot  {\widehat{\mathbf{n}}}_{f}}\right) \right|  - \left( {{\widehat{\mathbf{s}}}_{i} \cdot  {\widehat{\mathbf{n}}}_{f}}\right) }{2}{A}_{f}}\right) }\right\rbrack  {I}_{i,{N}_{f}} = {\beta }_{P}{S}_{i,P}{V}_{P},i = 1,2,...,n \tag{16.51}
$$




<!-- Meanless: The Method of Discrete Ordinates (S ${}_{N}$ -Approximation) Chapter $\mid$ 579-->

方程(16.51)是RTE的一般有限体积形式（用于空间离散化），其中使用有限差分（通过DOM）来离散化角度。它适用于任何形状或方向的控制体积（或单元格），只要该单元格是凸的并且由一组平面界定。它既适用于结构化网格也适用于非结构化网格，因为它不依赖于(i,j,k)的排序。在稍后的一个小节中，我们将为二维结构化正交网格（具有(i,j)排序）进一步简化它，以突出某些要点。

## 空间离散格式

方程(16.49)中使用的近似，在辐射文献中通常被称为阶梯格式，其精度仅为一阶。这基本上意味着，由于此近似而产生的误差将随着网格的细化呈线性比例变化（如果网格尺寸减半，误差将大约减少一半）。这可以通过泰勒级数展开轻松证明，如标准数值方法教科书[45]所示。为了提高精度，通常需要更高阶的近似。离散坐标法的早期工作研究了几种不同的插值（从单元格到面）格式以实现这一目标。最早尝试的高阶格式是所谓的菱形格式。菱形格式在均匀网格上是二阶精确的格式，基本上利用了中心差分格式的思想[45]，但用于外插而不是内插。根据均匀网格上的中心差分格式（见下图16.3），例如

$$
{I}_{P} \approx  \frac{{I}_{e} + {I}_{w}}{2} \tag{16.52}
$$

与其直接使用方程(16.52)，不如将其重新排列，用单元中心强度和西面强度来表示东面强度，得到

$$
{I}_{e} \approx  2{I}_{P} - {I}_{w}. \tag{16.53}
$$

换句话说，不仅仅是将单元中心强度外插到下游面（如阶梯格式中那样），而是进行一种更高阶的外插，其中使用一个单元中心值和一个上游面值进行外插。从方程(16.52)到方程(16.53)的重排是基于所讨论的坐标方向具有正方向余弦的假设。如果方向余弦 ${\xi }_{i}$ 为负，则方程(16.53)中的 ${I}_{e}$ 和 ${I}_{w}$ 需要交换。同样，可以在北面和南面之间进行相同的外插，得到

$$
{I}_{n} \approx  2{I}_{P} - {I}_{s} \tag{16.54}
$$

对于非均匀网格，方程(16.53)和(16.54)可以推广为

$$
{I}_{e} \approx  \frac{{I}_{P} - \left( {1 - {\gamma }_{x}}\right) {I}_{w}}{{\gamma }_{x}};\;{I}_{n} \approx  \frac{{I}_{P} - \left( {1 - {\gamma }_{y}}\right) {I}_{s}}{{\gamma }_{y}}, \tag{16.55}
$$

其中 $\frac{1}{2} \leq  {\gamma }_{x},{\gamma }_{y} \leq  1$。Carlson和Lathrop [4]将此格式称为加权菱形格式。对于均匀网格（菱形格式），${\gamma }_{x} = {\gamma }_{y} = \frac{1}{2}$。通常，${\gamma }_{x}$ 和 ${\gamma }_{y}$ 可以使用适当的距离加权插值在局部确定。为了节省计算时间，即使对于非均匀网格，也通常使用 ${\gamma }_{x} = {\gamma }_{y} = \frac{1}{2}$。虽然菱形格式有望提高计算的准确性，但研究人员发现，使用此格式有时会导致物理上不切实际的负强度[4]。Fiveland [13] 提出了通过在某些约束条件下使用单元尺寸来最小化（但不能完全避免）这种负强度的方法。然而，Chai及其同事 [46] 证明了细网格并不能保证正强度，实际上可能导致负强度。他们进一步注意到，菱形格式可能导致“过冲”，即预测出物理上不合理的高强度（离开控制体积的强度大于进入的强度加上内部发射）。最后一点，菱形格式需要存在“对面”，这一要求只能由结构化网格满足。即使对于这样的网格，如果所讨论的单元格高度倾斜，“对面”也不一定是“上游”和“下游”面。尽管如此，菱形格式在使用离散坐标法进行辐射输运的早期计算中得到了广泛应用。除了阶梯格式和菱形格式外，还使用过其他相对简单的格式，如正值格式（Lathrop [47]）、指数格式（Carlson和Lathrop [4]）、可变权重格式（Jamaluddin和Smith [48]）、上游追踪格式（Chai及其同事 [46]）和混合格式（Kim和Kim [49]）。接下来，将通过一个例子展示如何对一个简单的一维问题使用阶梯格式和菱形格式。


<!-- Meanless: 80 Radiative Heat Transfer-->

<!-- Media -->

<!-- figureText: $k + 1$ $k = K$ -->

<img src="https://cdn.noedgeai.com/bo_d141kgv7aajc73bva480_17.jpg?x=585&y=230&w=664&h=177&r=0"/>

图 16.3 例16.4中考虑的计算域，显示了单元格、面和所使用的索引方案。

<!-- Media -->

例16.4. 考虑一个厚度为 $L = 1\mathrm{\;m}$ 的平面平行介质（限制在 $x = 0$ 和 $x = L$ 之间），辐射沿着与正 $x$ 轴对齐的坐标方向传播。介质的吸收系数为 $\kappa  = 1{\mathrm{\;m}}^{-1}$，使得光学厚度 ${\kappa L} = 1$，并且介质不散射。介质的发射功率（温度）是固定的，由 ${E}_{b} = \pi \sin \left( {3x}\right)$ 给出。使用阶梯格式和菱形格式，计算从 $x = 0$ 到 $x = L$ 的视线上的强度，已知进入强度（在 $x = 0$ 处）等于 $1\mathrm{\;W}/{\mathrm{m}}^{2}\mathrm{{sr}}$。使用5个单元格和10个单元格跨越计算域进行两次单独的计算。将您的结果与解析解进行比较，以评估两种格式的准确性。评论您的结果。

## 解

对于一个吸收-发射介质，沿着指向正 $x$ 方向的视线，RTE可以写成

$$
\frac{dI}{dx} + {\kappa I} = \kappa {I}_{b}
$$

代入 $\kappa  = 1$ 和 ${I}_{b} = \sin \left( {3x}\right)$，我们得到

$$
\frac{dI}{dx} + I = \sin \left( {3x}\right) 
$$

对上述一阶线性常微分方程进行解析求解，我们得到

$$
I\left( x\right)  = C\exp \left( {-x}\right)  + \frac{\sin \left( {3x}\right)  - 3\cos \left( {3x}\right) }{10},
$$

其中 $C$ 是积分常数。使用边界条件 $I\left( 0\right)  = 1$，我们得到 $C = {13}/{10}$。

接下来我们考虑使用有限体积法的数值解。首先，我们将计算域离散为 $K$ 个等尺寸的单元格，使得单元格宽度（或网格尺寸）由 ${\Delta x} = L/K$ 给出。图16.3显示了离散化域的示意图。注意，单元格的东面和西面分别用 $e$ 和 $w$ 表示。为了推导必要的有限体积方程，我们将上述控制微分方程在一个一维单元格上积分，该单元格的中心用 $P$ 表示（索引为 $k$）。单元格 $P$ 的东面和西面相邻单元格分别用 $E$ 和 $W$ 表示。这得到

$$
{\int }_{w}^{e}\frac{dI}{dx}{dx} + {\int }_{w}^{e}{Idx} = {\int }_{w}^{e}\sin \left( {3x}\right) {dx}.
$$

对第一项进行精确积分，并对其他两项调用前面讨论过的（与有限体积法相关的）近似，我们得到

$$
{I}_{e} - {I}_{w} + {I}_{P}{\Delta x} = \sin \left( {3{x}_{P}}\right) {\Delta x},
$$

其中 ${x}_{P}$ 是单元中心 $P$ 的 $x$ 坐标。上述方程是所有单元格的基本有限体积方程，并且无论使用何种离散化格式，都将从此以后使用。推导相同方程的另一种方法是从方程(16.48)开始，并将其简化为具有面 $e$ 和 $w$ 的一维单元格，使得 ${N}_{f,P} = 2,{\widehat{\mathbf{n}}}_{e} = \widehat{\mathbf{r}}$ 和 ${\widehat{\mathbf{n}}}_{w} =  - \widehat{\mathbf{r}}$。此外，使用 ${\widehat{\mathbf{s}}}_{i} = \widehat{\mathbf{r}},{\beta }_{P} = 1,{S}_{i,P} = \sin \left( {3{x}_{P}}\right)$ 和 ${V}_{P}/{A}_{w} = {V}_{P}/{A}_{e} = {\Delta x}$，可以很容易地得到上述方程。

在阶梯格式中，面的强度值被上游单元格的值所替代，如方程(16.49)所示。因此，${I}_{e} \approx  {I}_{P}$ 和 ${I}_{w} \approx  {I}_{W}$。代入上述方程，然后重新整理，我们得到

$$
{I}_{P} = {I}_{k} = \frac{{I}_{W} + \sin \left( {3{x}_{P}}\right) {\Delta x}}{1 + {\Delta x}} = \frac{{I}_{k - 1} + \sin \left( {3{x}_{k}}\right) {\Delta x}}{1 + {\Delta x}},\;k = 2,3,\ldots ,K.
$$




<!-- Meanless: The Method of Discrete Ordinates (S ${}_{N}$ -Approximation) Chapter $\mid$-->

<!-- Media -->

<!-- figureText: 1 阶梯格式, $K = 5$ 阶梯格式, $K = {10}$ 菱形格式, $K = 5$ 菱形格式, $K = {10}$ 0.6 0.8 $x,\mathrm{\;m}$ 0.95 强度, $\mathrm{W}/{\mathrm{m}}^{2}$ sr 0.9 0.85 精确解 0.8 0.75 0.2 0.4 -->

<img src="https://cdn.noedgeai.com/bo_d141kgv7aajc73bva480_18.jpg?x=590&y=228&w=662&h=651&r=0"/>

图 16.4 使用阶梯格式和菱形格式对例16.4中考虑的问题在两种不同网格尺寸下得到的结果。

<!-- Media -->

对于第一个单元格 $\left( {k = 1}\right)$，必须使用左端的边界条件，即 ${I}_{w} = I\left( 0\right)$。这导致

$$
{I}_{P} = {I}_{1} = \frac{I\left( 0\right)  + \sin \left( {3{x}_{P}}\right) {\Delta x}}{1 + {\Delta x}} = \frac{I\left( 0\right)  + \sin \left( {3{x}_{1}}\right) {\Delta x}}{1 + {\Delta x}}.
$$

上述两个方程代表一组解耦的（或显式的）方程，可以逐一求解。首先求解第二个方程以确定 ${I}_{1}$，然后将其代入另一个方程以确定 ${I}_{2},{I}_{3}$ 等。

在菱形格式中，下游面的强度值被上游单元中心值和上游面值的组合所替代，如方程(16.53)所示。因此，${I}_{e} \approx  2{I}_{P} - {I}_{w}$。代入控制有限体积方程，并切换到带编号的索引，我们得到

$$
{I}_{k} = \frac{2{I}_{w,k} + \sin \left( {3{x}_{k}}\right) {\Delta x}}{2 + {\Delta x}}\;k = 2,3,\ldots ,K.
$$

其中 ${I}_{w,k}$ 表示单元格 $k$ 的西面强度。对于第一个单元格 $\left( {k = 1}\right)$，使用左边界条件，我们得到

$$
{I}_{1} = \frac{{2I}\left( 0\right)  + \sin \left( {3{x}_{1}}\right) {\Delta x}}{2 + {\Delta x}}
$$

这里，上述两个方程也代表一组解耦的（或显式的）方程，可以逐一求解。首先求解第二个方程以确定 ${I}_{1}$。接下来，使用方程 ${I}_{e,1} \approx  2{I}_{1} - I\left( 0\right)$ 来确定第一个单元格东面的强度。注意到 ${I}_{w,k + 1} = {I}_{e,k}$，解可以向前推进 $k$，使用另一个方程来确定 ${I}_{2},{I}_{3}$ 等。

图16.4显示了使用两种不同格式和两种不同网格尺寸（5个单元格 vs. 10个单元格）得到的结果。仅用5个单元格，阶梯格式产生的结果不是很准确。然而，即使在如此粗糙的网格下，强度的定性变化也被遵循：下降-上升-下降。对于相同的网格，菱形格式产生的结果甚至更准确。当网格细化时（10个单元格），阶梯格式的结果紧密地跟随精确解。菱形格式的结果从粗网格到细网格的准确性有显著提高。菱形格式是二阶精确的，这意味着当网格间距减半时，误差大约减少四倍，这解释了当网格细化时准确性的急剧提高。

$\nabla  \cdot  \left( {I\widehat{\mathbf{s}}}\right)$ 项与流体流动中的对流通量项具有很强的相似性。从CFD文献中众所周知，为减少离散化（或截断）误差而寻求的高阶对流格式，通常会因缺乏有界性而受到影响[50,51]。物理上，一个面上的强度必须受限于（在没有任何其他辐射源的情况下）跨越该面的两个单元中心处的强度。不幸的是，当使用高阶插值时，可能会出现局部最大值或最小值，导致违反此约束。强度可能不会从一个单元到另一个单元局部单调。除了阶梯格式（一阶迎风）外，没有其他格式能保证有界性。为了寻求既能保持有界性又能实现高阶精度的格式，辐射领域的早期研究者转向了CFD文献，并探索了几种所谓的有界性保持格式。有界性保持格式基于Leonard [51]首次提出的归一化变量图（NVD）的思想。有界性保持格式的核心思想是将高阶格式与一阶迎风（或阶梯）格式混合，使得 ${I}_{f} = {I}_{P} + {\beta }_{b}\left( {{\left\lbrack  {I}_{f}\right\rbrack  }_{\text{高阶 }} - {I}_{P}}\right)$，其中 ${\beta }_{b}$ 是一个混合因子。如果混合因子为1，则采用纯高阶格式（最高精度），而如果混合因子为0，则格式默认为阶梯格式（最低精度）。介于两者之间，精度是中等的。混合因子根据局部解的单调性程度在局部计算，并使用NVD构建一组约束。不同有界性保持格式的区别在于确定 ${\beta }_{b}$ 的方式（精确的数学关系）。


<!-- Meanless: 582 Radiative Heat Transfer-->

Jessee和Fiveland [52]早期为辐射计算探索的一种有界性保持格式是CLAM（曲线对流法）格式[50]。在此格式中，

$$
{I}_{f} = \left\{  \begin{array}{ll} {I}_{P} + \widetilde{I}\left( {{I}_{d} - {I}_{P}}\right) , & 0 \leq  \widetilde{I} \leq  1, \\  {I}_{P}, & \text{ 否则,} \end{array}\right.  \tag{16.56}
$$

其中 $\widetilde{I}$ 由下式确定

$$
\widetilde{I} = \frac{{I}_{P} - {I}_{u}}{{I}_{d} - {I}_{u}} \tag{16.57}
$$

其中 ${I}_{u}$ 是上游单元格的强度，${I}_{d}$ 是下游单元格的强度。在此格式中，$\widetilde{I}$ 是前述的归一化变量，表示强度在局部的行为，即单调性程度。当 $\widetilde{I}$ 超出0和1的界限（非单调）时，格式默认为纯阶梯格式，正如我们所讨论的，它总是有界的。另一方面，如果 $0 \leq  \widetilde{I} \leq  1$（单调），格式将阶梯格式（迎风差分）与菱形格式（中心差分）混合。如果 $\widetilde{I} = 1/2$，则恢复二阶菱形格式。通过以这种方式将一阶阶梯格式与二阶菱形格式混合，可以潜在地提高精度。例如，在没有强发射源和散射的介质中，强度很可能从一个边界到另一个边界单调变化，格式将在大多数单元格中默认为二阶，从而在不违反有界性约束的情况下提高精度。尽管CLAM格式在精度和有界性之间提供了一个很好的折衷，但它有一个严重的缺点。方程(16.56)显示CLAM格式是非线性的，因为 $\widetilde{I}$ 是强度本身的函数。因此，即使在没有散射和/或壁面反射的情况下，沿单一方向传播的强度也不能再通过一次扫描来评估。方程必须线性化，并通过迭代找到解。对于第 $n$ 次迭代，方程(16.56)修改为

$$
{I}_{f}^{n} = {I}_{P}^{n} + {\widetilde{I}}^{n - 1}\left( {{I}_{d}^{n} - {I}_{P}^{n}}\right) , \tag{16.58}
$$

其中上标表示迭代指数。$\widetilde{I}$ 可能从一次迭代到下一次迭代发生变化，这一事实使得这种格式本质上有些不稳定。

许多其他有界性保持的高阶格式已经从CFD文献中被采纳并用于使用DOM求解RTE。Liu等人[53]应用了Gaskell和Lau [54]的SMART格式——一种三阶精确的格式。他们发现SMART格式产生了准确、无振荡和正的辐射强度分布。与阶梯、菱形和其他格式的比较表明，虽然辐射强度对空间差分格式非常敏感，但积分量，如入射辐射，则不那么敏感。Jessee和Fiveland [52]应用了MINMOD [55]、MUSCL [56]、CLAM [50]和SMART [54]格式，他们的发现与Liu等人[53]的相似。在一项广泛的比较研究中，Coelho [57]除了已经提到的格式外，还测试了大量其他格式，包括几种所谓的总变差递减（TVD）格式和基本无振荡（ENO）格式。考虑了几种有散射和无散射的多维问题。虽然没有一种格式表现出明显的优越性，但得出的结论是，CLAM格式不如大多数高阶格式准确，但它相对稳定且计算成本较低。


<!-- Meanless: The Method of Discrete Ordinates (S ${}_{N}$ -Approximation) Chapter $\mid$-->

<!-- Media -->

<!-- figureText: _____ $\left( {j,k + 1}\right)$ $P$ $E$ (j,k) $\left( {j + 1,k}\right)$ $S$ (j,k - 1) ${\Delta y}$ ${\Delta x}$ $x,j$ $W$ (j - 1,k) (1,2) (1,1) (2,1) -->

<img src="https://cdn.noedgeai.com/bo_d141kgv7aajc73bva480_20.jpg?x=582&y=226&w=667&h=601&r=0"/>

图 16.5 二维正交结构化网格的单元格、面和索引模式的排列。

<!-- Media -->

虽然高阶有界性保持格式无疑比其低阶对应格式产生更好的结果，但当坐标方向与网格线不对齐时，其精度往往会受到影响。为了解决这个问题，Coelho [58] 提出了所谓的斜高阶分辨率（SHR）格式。虽然这种格式在多维问题上的精度提高得到了确凿的证明，但计算仅在简单几何形状（正方形和立方体）上并使用正交结构化网格进行。至今，将高阶格式用于非结构化网格或结构化曲线网格的DOM计算仍然具有挑战性。最近，Coelho [59] 在二维非结构化（三角形）网格中实现了高阶TVD格式用于DOM计算。这项研究发现，由于在非结构化网格中需要做出额外的近似来定位“上游”单元格，结果不如同样格式在笛卡尔网格上的实现准确。此外，发现收敛性比笛卡尔网格差。更多关于空间差分格式的细节可以在[18,52,57,59]中找到。

## 二维正交网格

为清晰起见，我们将再次为二维几何（即 $\partial I/\partial z \equiv  0$）发展该方法，该几何使用如图16.5所示的结构化二维均匀正交网格进行离散化。每个控制体积（或单元格）现在有四个面，我们将它们表示为 $e,w,n$ 和 $s$，而相应的相邻单元格分别表示为 $E,W,N$ 和 $S$。为简单起见，首先考虑阶梯格式，因为已经证明阶梯格式最容易实现，并且总是产生有界的结果，尽管不是很准确。然后将发展使用菱形格式的方程。对于阶梯格式，我们发展的起点是由方程(16.51)给出的一般有限体积方程。我们注意到 ${\widehat{\mathbf{n}}}_{e} = \widehat{\mathbf{1}},{\widehat{\mathbf{n}}}_{w} =  - \widehat{\mathbf{1}},{\widehat{\mathbf{n}}}_{n} = \widehat{\mathbf{J}}$ 和 ${\widehat{\mathbf{n}}}_{s} =  - \widehat{\mathbf{J}}$。因此，${\widehat{\mathbf{s}}}_{i} \cdot  {\widehat{\mathbf{n}}}_{e} = {\xi }_{i},{\widehat{\mathbf{s}}}_{i} \cdot  {\widehat{\mathbf{n}}}_{w} =  - {\xi }_{i},{\widehat{\mathbf{s}}}_{i} \cdot  {\widehat{\mathbf{n}}}_{n} = {\eta }_{i}$ 和 ${\widehat{\mathbf{s}}}_{i} \cdot  {\widehat{\mathbf{n}}}_{s} =  - {\eta }_{i}$。此外，${A}_{e} = {A}_{w} = {\Delta y},{A}_{n} = {A}_{s} = {\Delta x}$ 和 ${V}_{P} = {\Delta x\Delta y}$。将所有这些量代入方程(16.51)并在将求和展开为四项后，得到

$$
\left\lbrack  {\frac{\left| {\xi }_{i}\right|  + {\xi }_{i}}{2}{\Delta y} + \frac{\left| {\xi }_{i}\right|  - {\xi }_{i}}{2}{\Delta y} + \frac{\left| {\eta }_{i}\right|  + {\eta }_{i}}{2}{\Delta x} + \frac{\left| {\eta }_{i}\right|  - {\eta }_{i}}{2}{\Delta x} + {\beta }_{P}{\Delta x\Delta y}}\right\rbrack  {I}_{i,P} - \left\lbrack  {\frac{\left| {\xi }_{i}\right|  - {\xi }_{i}}{2}{\Delta y}}\right\rbrack  {I}_{i,E} - \left\lbrack  {\frac{\left| {\xi }_{i}\right|  + {\xi }_{i}}{2}{\Delta y}}\right\rbrack  {I}_{i,W}
$$

$$
 - \left\lbrack  {\frac{\left| {\eta }_{i}\right|  - {\eta }_{i}}{2}{\Delta x}}\right\rbrack  {I}_{i,N} - \left\lbrack  {\frac{\left| {\eta }_{i}\right|  + {\eta }_{i}}{2}{\Delta x}}\right\rbrack  {I}_{i,S} = {\beta }_{P}{S}_{i,P}{\Delta x\Delta y},\;i = 1,2,\ldots ,n. \tag{16.59}
$$




<!-- Meanless: 584 Radiative Heat Transfer-->

将强度、消光系数和源项从方向性(P,E,W等)索引转换为编号(j,k)索引，方程(16.59)可以写成以下通用形式：

$$
{C}_{P}{I}_{i,j,k} + {C}_{E}{I}_{i,j + 1,k} + {C}_{W}{I}_{i,j - 1,k} + {C}_{N}{I}_{i,j,k + 1} + {C}_{S}{I}_{i,j,k - 1} = {R}_{i,j,k},i = 1,2,\ldots ,n;j = 2,3,\ldots ,J - 1;k = 2,3,\ldots ,K - 1\text{,}
$$

(16.60)

其中各个系数和右侧可以写为

$$
{C}_{P} = \frac{\left| {\xi }_{i}\right|  + {\xi }_{i}}{2}{\Delta y} + \frac{\left| {\xi }_{i}\right|  - {\xi }_{i}}{2}{\Delta y} + \frac{\left| {\eta }_{i}\right|  + {\eta }_{i}}{2}{\Delta x} + \frac{\left| {\eta }_{i}\right|  - {\eta }_{i}}{2}{\Delta x} + {\beta }_{j,k}{\Delta x\Delta y};
$$

$$
{C}_{E} =  - \frac{\left| {\xi }_{i}\right|  - {\xi }_{i}}{2}{\Delta y};{C}_{W} =  - \frac{\left| {\xi }_{i}\right|  + {\xi }_{i}}{2}{\Delta y};{C}_{N} =  - \frac{\left| {\eta }_{i}\right|  - {\eta }_{i}}{2}{\Delta x};{C}_{S} =  - \frac{\left| {\eta }_{i}\right|  + {\eta }_{i}}{2}{\Delta x};{R}_{i,j,k} = {\beta }_{j,k}{S}_{i,j,k}{\Delta x\Delta y}. \tag{16.61}
$$

初看方程(16.60)，似乎所有四个相邻单元格都对单元格 $P\left\lbrack  {\text{或}\left( {j,k}\right) }\right\rbrack$ 的强度有贡献。进一步检查发现，只有四个邻居中的两个会影响 $P$ 处的解。例如，如果一个坐标方向使得射线在第一象限向前传播，那么 ${\xi }_{i} > 0$ 且 ${\eta }_{i} > 0$。因此，方程(16.61)简化为

$$
{C}_{P} = {\xi }_{i}{\Delta y} + {\eta }_{i}{\Delta x} + {\beta }_{j,k}{\Delta x\Delta y}
$$

$$
{C}_{E} = 0;\;{C}_{W} =  - {\xi }_{i}{\Delta y};\;{C}_{N} = 0;\;{C}_{S} =  - {\eta }_{i}{\Delta x};\;{R}_{i,j,k} = {\beta }_{j,k}{S}_{i,j,k}{\Delta x\Delta y}. \tag{16.62}
$$

因此，方程(16.60)简化为

$$
{C}_{P}{I}_{i,j,k} + {C}_{W}{I}_{i,j - 1,k} + {C}_{S}{I}_{i,j,k - 1} = {R}_{i,j,k},\;i \in  \left( {{\xi }_{i} > 0,{\eta }_{i} > 0}\right) ;j = 2,3,\ldots ,J - 1;k = 2,3,\ldots ,K - 1\text{,} \tag{16.63}
$$

换句话说，方程(16.60)是一个通用方程，它自动考虑了坐标方向及其相对于坐标系的方向，并且其写法使得四个系数中的两个总是会消失。

如上所示，方程(16.60)仅对内部单元格有效。对于与边界相邻的单元格，方程(16.48)中对面的求和必须分成两部分：一部分是对内部面的求和，另一部分是对边界面。为了说明边界条件的应用，我们现在考虑一个角单元格，即单元格(1,1)。它的南面和西面位于边界上。因此，对于一个坐标方向为 ${\xi }_{i} > 0$ 和 ${\eta }_{i} > 0$ 的情况，这个单元格的有限体积方程，在简化方程(16.48)后，变为

$$
\left( {{\xi }_{i}{\Delta y} + {\eta }_{i}{\Delta x} + {\beta }_{1,1}{\Delta x\Delta y}}\right) {I}_{i,1,1} - {\xi }_{i}\frac{{J}_{w}}{\pi }{\Delta y} - {\eta }_{i}\frac{{J}_{s}}{\pi }{\Delta x} = {\beta }_{1,1}{S}_{i,1,1}{\Delta x\Delta y}, \tag{16.64}
$$

其中，西面和南面的边界条件，即 ${I}_{i,w} = {J}_{w}/\pi$ 和 ${I}_{i,s} = {J}_{s}/\pi$，已经与内部面(n和e)的阶梯格式一起使用。如果墙壁是黑体，那么 ${I}_{i,w} = {I}_{bw}$ 和 ${I}_{i,s} = {I}_{bs}$，方程(16.64)可以进一步简化为

$$
\left( {{\xi }_{i}{\Delta y} + {\eta }_{i}{\Delta x} + {\beta }_{1,1}{\Delta x\Delta y}}\right) {I}_{i,1,1} = {\beta }_{1,1}{S}_{i,1,1}{\Delta x\Delta y} + {\xi }_{i}{I}_{bw}{\Delta y} + {\eta }_{i}{I}_{bs}{\Delta x}. \tag{16.65}
$$

此外，如果介质不散射，源项 ${S}_{i,1,1} = {I}_{b,1,1}$ 是一个已知量，如果介质的温度是给定的。因此，方程(16.65)可以重新排列写成一个显式方程如下：

$$
{I}_{i,1,1} = \frac{{\beta }_{1,1}{I}_{b,1,1}{\Delta x\Delta y} + {\xi }_{i}{I}_{bw}{\Delta y} + {\eta }_{i}{I}_{bs}{\Delta x}}{\left( {\xi }_{i}\Delta y + {\eta }_{i}\Delta x + {\beta }_{1,1}\Delta x\Delta y\right) }, \tag{16.66}
$$

其中右侧的所有量都是已知的。一旦(1,1)和所有其他西南边界单元格中的强度被确定，其他下游内部单元格的强度可以很容易地使用方程(16.63)确定，方法也是将其写成显式形式，其中南部和西部单元格中心强度已经确定。

总之，如果所有墙壁都是黑体且没有散射，所有的单元中心强度都可以通过一次遍历计算出来，因为所有的墙壁辐射出射度和所有内部源项都是先验已知的（如果温度场是给定或假设的）。如果墙壁是反射的并且/或者介质是散射的，则需要迭代。在对所有方向和所有单元格完成一次遍历后，墙壁辐射出射度和辐射源项的值将被更新，然后重复此过程直到满足收敛准则。最后，入射辐射和辐射热流的内部值由方程(16.6)和(16.7)确定，而墙壁上的热流可以由方程(16.8)计算。对于高反射性墙壁 $\left( {{\epsilon }_{w} \ll  1}\right)$ 和强散射介质 $\left( {1 - \omega  \ll  1}\right)$，离散坐标法将变得极其低效。正如Chai及其同事[60]指出的，由散射引起的迭代次数可以通过从相函数中移除前向散射并将其视为透射来减少。这可以在方程(16.38)和(16.39)中通过定义一个修正的消光系数和一个修正的源项来实现


<!-- Meanless: The Method of Discrete Ordinates (S ${}_{N}$ -Approximation) Chapter $\mid$-->

$$
{\beta }_{mi} = \beta  - \frac{{\sigma }_{s}}{4\pi }{w}_{i}{\Phi }_{ii} \tag{16.67}
$$

$$
{S}_{mi} = \left( {1 - \omega }\right) {I}_{b} + \frac{\omega }{4\pi }\mathop{\sum }\limits_{\substack{{j = 1} \\  {j \neq  i} }}^{n}{w}_{j}{\Phi }_{ij}{I}_{j},\;i = 1,2,\ldots ,n. \tag{16.68}
$$

这导致更快的收敛，特别是如果相函数具有很强的前向峰值（这对于大颗粒通常是这种情况；另见第11.9节的讨论）。最后需要注意的一点是，对于本讨论中考虑的二维问题，正负 ${\mu }_{i}$ 值的强度是相同的。因此，我们只需要考虑 ${\mu }_{i}$ 的正值（并将求积权重 ${w}_{i}$ 加倍）。

各向异性散射的处理在DOM计算中受到了特别的关注。许多研究[61-64]表明，如果散射相函数具有很强的前向峰值，DOM中使用的常规方向求积格式不能保守散射能量。一个已经提出的纠正措施是归一化散射相函数[61]。虽然这种方法有效地保守了散射能量，但Boulet等人[63]表明，它改变了相函数的整体形状和不对称因子，与蒙特卡罗方法生成的结果相比导致了很大的差异。Hunter和Guo [65-68]提出并演示了一种新的散射相函数归一化技术，它不仅保守散射能量，而且还保留了相函数的形状。

在存在壁面反射和/或散射的情况下，离散坐标方程（方向RTEs）变得强耦合。这不仅对收敛不利，而且对求解器的并行化也存在问题（例如，如果对基于方向的分布式内存并行化感兴趣）。为了解决方向RTEs之间的强耦合问题，Murthy和Mathur [69,70]提出了所谓的耦合坐标法（COMET）。在这种方法中，由方向RTEs离散化产生的线性代数方程是同时为所有方向求解的，而不是一次一个。虽然这种方法的收敛性显著提高，但这种耦合求解器需要更多的内存，并且在计算机程序中实现起来很困难。最近，Wang等人[71]提出了所谓的散射降阶（ROS）方法。在这种方法中，强度被写成各种阶数强度的叠加——每个阶数的方程代表与散射相关的某个物理事件，并且大小逐渐减小。例如，零阶是最大的分量，它没有入散射，因此可以无迭代耦合地求解。研究表明，2-5阶足以恢复精确的强度。该方法提供的计算增益尚不清楚。

接下来，我们推导使用菱形格式更新强度所需的方程。我们从方程(16.48)开始推导。在与上述正交网格相同的简化条件下，我们得到

$$
{\xi }_{i}\left( {{I}_{i,e} - {I}_{i,w}}\right) {\Delta y} + {\eta }_{i}\left( {{I}_{i,n} - {I}_{i,s}}\right) {\Delta x} =  - {\beta }_{P}{I}_{i,P}{\Delta x\Delta y} + {\beta }_{P}{S}_{i,P}{\Delta x\Delta y}. \tag{16.69}
$$

再次，如同为阶梯格式推导方程(16.59)时所做的那样，使用切换函数，我们可以写出一个适用于任何任意方向的方程：

$$
\left\lbrack  {\left( {\left| {\xi }_{i}\right|  + {\xi }_{i}}\right) {\Delta y} + \left( {\left| {\xi }_{i}\right|  - {\xi }_{i}}\right) {\Delta y} + \left( {\left| {\eta }_{i}\right|  + {\eta }_{i}}\right) {\Delta x} + \left( {\left| {\eta }_{i}\right|  - {\eta }_{i}}\right) {\Delta x} + {\beta }_{P}{\Delta x\Delta y}}\right\rbrack  {I}_{i,P}
$$

$$
 - \left\lbrack  {\left( {\left| {\xi }_{i}\right|  - {\xi }_{i}}\right) {\Delta y}}\right\rbrack  {I}_{i,e} - \left\lbrack  {\left( {\left| {\xi }_{i}\right|  + {\xi }_{i}}\right) {\Delta y}}\right\rbrack  {I}_{i,w}
$$

$$
 - \left\lbrack  {\left( {\left| {\eta }_{i}\right|  - {\eta }_{i}}\right) {\Delta x}}\right\rbrack  {I}_{i,n} - \left\lbrack  {\left( {\left| {\eta }_{i}\right|  + {\eta }_{i}}\right) {\Delta x}}\right\rbrack  {I}_{i,s} = {\beta }_{P}{S}_{i,P}{\Delta x\Delta y},\;i = 1,2,\ldots ,n. \tag{16.70}
$$




<!-- Meanless: 586 Radiative Heat Transfer-->

其中，使用了由方程(16.53)和(16.54)给出的菱形格式。方程(16.70)可以重新排列，写出对任何 $\left( {{\xi }_{i},{\eta }_{i}}\right)$ 组合都成立的方程：

$$
{I}_{i,P} = \frac{{\beta }_{P}{S}_{i,P}{\Delta x\Delta y} + \left( {\left| {\xi }_{i}\right|  - {\xi }_{i}}\right) {\Delta y}{I}_{i,e} + \left( {\left| {\xi }_{i}\right|  + {\xi }_{i}}\right) {\Delta y}{I}_{i,w} + \left( {\left| {\eta }_{i}\right|  - {\eta }_{i}}\right) {\Delta x}{I}_{i,n} + \left( {\left| {\eta }_{i}\right|  + {\eta }_{i}}\right) {\Delta x}{I}_{i,s}}{{\beta }_{P}{\Delta x\Delta y} + \left( {\left| {\xi }_{i}\right|  + {\xi }_{i}}\right) {\Delta y} + \left( {\left| {\xi }_{i}\right|  - {\xi }_{i}}\right) {\Delta y} + \left( {\left| {\eta }_{i}\right|  + {\eta }_{i}}\right) {\Delta x} + \left( {\left| {\eta }_{i}\right|  - {\eta }_{i}}\right) {\Delta x}},\;i = 1,2,\ldots ,n.
$$

(16.71)

如前所述，对于阶梯格式，无论选择何种坐标方向，四个相邻贡献中的两个总是会消失。例如，对于在第一象限向前传播的射线 $\left( {{\xi }_{i} > 0\text{且}{\eta }_{i} > 0}\right)$，方程(16.71)简化为

$$
{I}_{i,P} = {I}_{i,j,k} = \frac{{\beta }_{j,k}{S}_{i,j,k}{\Delta x\Delta y} + 2{\xi }_{i}{I}_{i,w}{\Delta y} + 2{\eta }_{i}{I}_{i,s}{\Delta x}}{{\beta }_{j,k}{\Delta x\Delta y} + 2{\xi }_{i}{\Delta y} + 2{\eta }_{i}{\Delta x}}. \tag{16.72}
$$

方程(16.72)中的西面和南面值是上游单元格的东面和北面值。与阶梯格式一样，对于在第一象限向前传播的射线，计算从西南角开始。对于单元格(1,1)，西面和南面的强度值被边界条件替换，如上所述。这允许计算 ${I}_{i,1,1}$。一旦计算出 ${I}_{i,1,1}$，它和(1,1)的西面值可以用来使用方程(16.53)计算(1,1)的东面值。(2,1)的西面值等于(1,1)的东面值。因此，对于(2,1)，西面和南面（边界条件）都变得可用，(2,1)处的强度可以使用方程(16.72)计算。对所有其他下游单元格重复此过程，直到每个单元格的强度都被计算出来。这个过程将在接下来的例子中进一步演示。

对于像这里考虑的规则几何体，很明显，对于具有正方向余弦的坐标，强度计算必须从西南角开始。计算过程中有一个固有的假设，即南面和西面先被更新。在刚刚讨论的过程中，西南角对应于单元格(1,1)。通常，编号可能不从西南角开始。让我们考虑一个场景，其中单元格(1,1)是西北角单元格，而其南侧的单元格编号为(1,2)。在这种特殊情况下，前面的计算过程将失败，因为(1,2)中的强度将在单元格(1,1)之后更新，而计算过程要求它在单元格(1,1)之前更新。换句话说，单元格的编号模式和必须用来更新单元格的顺序可能不总是一致的。对于由多个块描述的复杂几何体，编号模式可能完全被打乱，并且可能与所需的更新顺序有非常复杂的关系。正是由于这个原因，大多数为DOM进行辐射计算而设计的现代多维代码不使用前述的显式逐点更新过程。而是使用迭代线性代数方程求解器来同时求解所有单元格的方程(16.60)并获得强度。这种方法也有利于在没有排序的非结构化网格上获得解。

关于方向耦合，最常见的方法是按顺序求解方向RTEs，然后使用外迭代循环将它们耦合起来。这种方法的优点是，为一个方向分配的用于存储系数矩阵的内存可以为其他方向重用。这种方法也适合并行处理，其中不同的方向RTEs可以分配给不同的处理器。如前所述，这种方向解耦（或分离）求解的缺点是，如果介质是强散射的和/或壁面是强反射的，收敛性可能会很差。

## 虚假散射

离散坐标法的一个严重缺点是虚假散射。当坐标方向与网格线不对齐时，会发生虚假散射。例如，如果通过离散坐标法追踪一束单一的准直光束穿过一个外壳，随着光束离其原点越来越远，它会逐渐变宽。这种即使在没有真实散射的情况下，辐射强度的非物理性涂抹现象被称为虚假散射，可以通过使用更细的网格或更高阶的方法来减少。这在接下来的例子中有所说明。

例16.5. 考虑一个由 $3 \times  3$ 网格离散化的正方形计算域。一束强度等于1的准直光束以 ${45}^{ \circ  }$ 的角度从正方形的底边进入，如图 ${16.6a}$ 所示。介质是完全


<!-- Meanless: The Method of Discrete Ordinates ( ${S}_{N}$ -Approximation) Chapter $\mid$-->

<!-- Media -->

<!-- figureText: $y$ $y$ 0.3125 0.5 0.5 (2,3) (3,3) (1,3) (2,3) (3,3) 0.6875 0.5 (2,2) (3,2) (1,2) (2,2) (3,2) 0.75 0.875 0.5 (2,1) (3,1) (1,1) (2,1) (3,1) (b) (c) $y$ 0.125 $I = 0$ $I = 0$ (1,3) 0.25 (1,2) 0.5 (1,1) 平行光束, $I = 1$ (a) -->

<img src="https://cdn.noedgeai.com/bo_d141kgv7aajc73bva480_24.jpg?x=196&y=222&w=1444&h=481&r=0"/>

图 16.6 例16.5中考虑的问题：(a) 几何形状和边界条件，以及精确解，(b) 使用阶梯格式的解，和 (c) 使用菱形格式的解。

<!-- Media -->

透明的。使用阶梯格式和菱形格式，计算计算域中9个单元格中心的强度。

## 解

由于光束以 ${45}^{ \circ  }$ 进入，${\xi }_{i} = {\eta }_{i}$。又因为介质完全透明，$\beta  = 0$。使用 ${\Delta x} = {\Delta y}$，方程(16.61)（阶梯格式）简化为

$$
{I}_{j,k} = \frac{{I}_{w,j,k} + {I}_{s,j,k}}{2} \approx  \frac{{I}_{j - 1,k} + {I}_{j,k - 1}}{2}.
$$

从西南角开始顺序应用上述公式，我们得到 ${I}_{1,1} = \left( {0 + 1}\right) /2 = {0.5}$，${I}_{2,1} = \left( {{I}_{1,1} + 1}\right) /2 = \left( {{0.5} + 1}\right) /2 = {0.75},{I}_{3,1} = \left( {{I}_{2,1} + 1}\right) /2 = \left( {{0.75} + 1}\right) /2 = {0.875},{I}_{1,2} = \left( {0 + {0.5}}\right) /2 = {0.25}$，等等。完整的单元中心强度值在图16.6b中以图形方式描绘。沿着对角线，我们有一个不连续性，任何连续的强度表示都无法解决。因此，0.5的值是不可避免的。然而，在其他单元格中，强度值偏离了1或0。在底部三角形中，强度低于1，而在顶部三角形中，它们大于0——几乎就像准直光束因散射而重新分配了其能量。因此，术语为虚假散射。如果使用 $6 \times  6$ 的网格，右下角单元格，即单元格(6,1)中的强度将为0.984375，而不是0.875。这表明虚假散射可以通过使用更细的网格来减轻，但不能完全消除。

为了推导菱形格式的相应更新方程，我们从方程(16.48)开始，并将其简化为具有 ${\xi }_{i} = {\eta }_{i},{\Delta x} = {\Delta y}$ 和 $\beta  = 0$ 的二维正交网格。这得到

$$
\left( {{I}_{e,P} - {I}_{w,P}}\right)  + \left( {{I}_{n,P} - {I}_{s,P}}\right)  = 0.
$$

将菱形格式，即 ${I}_{e,P} = 2{I}_{P} - {I}_{w,P}$ 和 ${I}_{n,P} = 2{I}_{P} - {I}_{s,P}$ 代入上述方程，然后重新排列，得到

$$
{I}_{P} = {I}_{j,k} = \frac{{I}_{w,j,k} + {I}_{s,j,k}}{2}.
$$

再次，从西南角开始顺序应用上述公式，我们得到 ${I}_{1,1} = \left( {0 + 1}\right) /2 =$ ${0.5},{I}_{e,1,1} = 2{I}_{1,1} - {I}_{w,1,1} = 2 \times  {0.5} - 0 = 1,{I}_{w,2,1} = {I}_{e,1,1} = 1,{I}_{2,1} = \left( {{I}_{w,2,1} + 1}\right) /2 = \left( {1 + 1}\right) /2 = 1,{I}_{e,2,1} = 2{I}_{2,1} - {I}_{w,2,1} = 2 \times  1 - 1 = 1,$ ${I}_{w,3,1} = {I}_{e,2,1} = 1,{I}_{3,1} = \left( {{I}_{w,3,1} + 1}\right) /2 = \left( {1 + 1}\right) /2 = 1$，等等。完整的单元格和面强度值在图16.6c中以图形方式描绘。在这种情况下，除了沿对角线的不连续性外，可以看出完全没有虚假散射的迹象。底部三角形中的所有面和单元中心强度都等于1，而顶部三角形中的所有强度都等于零。虽然菱形格式似乎完全消除了虚假散射，但这一发现仅对这种介质不参与的特殊情况成立。对于参与介质，Chai等人[72]表明菱形格式减轻了虚假散射，但并未完全消除它。最后一点，如果光束沿 $y$ 方向传播，例如，无论使用何种格式，各处的强度值都将等于1。换句话说，如果坐标方向与网格对齐，则消除了虚假散射。然而，在现实中，一些坐标方向总是会与网格错位，因此，虚假散射将是不可避免的。


<!-- Meanless: 588 Radiative Heat Transfer-->

<!-- Media -->

<!-- figureText: w,3 3 ${I}_{b} = 0$ e,3 $1\mathrm{w},4$ e,4 ${I}_{b} = 0$ 5,4 n,2 $e,1/w,2$ $e,2$ s,2 5,3 $n,1$ ${I}_{b} = 0$ w, 1 $s,1$ ${I}_{b} = {I}_{bw}$ -->

<img src="https://cdn.noedgeai.com/bo_d141kgv7aajc73bva480_25.jpg?x=554&y=224&w=723&h=561&r=0"/>

图 16.7 例16.6的方形外壳。

<!-- Media -->

总之，虚假散射是使用低阶格式（如阶梯格式）和强度沿与网格线成角度的方向传播的综合效应的表现。这种效应可以通过细化网格或使用更高阶的格式来减轻。

## 射线效应

离散坐标法的另一个严重缺点是所谓的“射线效应”，这是角度离散化的结果。考虑一个具有非常小区域（体积或表面积）且发射率非常高的外壳。来自该区域的强度将沿着离散坐标的方向传播出去。远离发射区，这些射线可能会变得如此分散，以至于一些控制体积和/或表面区域可能接收不到来自这个高发射区的任何能量，从而导致不物理的结果。显然，射线效应可以通过增加控制体积（更粗的网格）和表面区域的尺寸来减少。当网格粗糙时，过多的虚假散射实际上通过涂抹射线来帮助减轻射线效应。因此，当使用更精细的空间网格来减少虚假散射时，这应伴随着方法阶数的增加（即，更精细的角度求积）。射线效应和虚假散射已被许多研究人员调查。关于射线效应以及如何减轻它们的讨论可以在[72-78]中找到。在20世纪90年代初，提出了一种离散坐标法的新变体——所谓的辐射有限体积法，在本文中称为有限角法（FAM）。该方法在第16.6节中讨论。由于该方法与射线效应有很强的关系，对射线效应的进一步讨论将推迟到该方法被介绍之后。

例16.6. 一个灰色的、吸收/发射（但不散射）的介质包含在一个边长为 $L$ 的方形外壳内。该介质处于辐射平衡状态，并具有一个恒定的吸收系数，使得 ${\kappa L} = 1$。顶部和两侧的墙壁温度为零，而底部的墙壁是等温的，温度为 ${T}_{w}$（具有恒定的黑体强度 ${I}_{bw}$）；所有四个表面都是黑体。使用离散坐标法计算底部表面的局部热损失。

## 解

为了本例的说明目的，我们将自己限制在简单的非对称 ${S}_{2}$ -近似，并使用图16.7中所示的粗糙节点系统。对于非对称 ${S}_{2}$ -近似（在 $z$ 方向上没有依赖性），我们必须考虑四个离散坐标，其方向向量（投影到 $x$ -y平面上）为 ${\widehat{\mathbf{s}}}_{i} = {\xi }_{i}\widehat{\mathbf{i}} + {\eta }_{i}\widehat{\mathbf{j}} =  \pm  {0.5}\left( {\widehat{\mathbf{i}} \pm  \widehat{\mathbf{j}}}\right)$，如表16.1所示。每个方向的求积权重，由于是二维的，加倍后为 ${w}_{i} = \pi$。对于灰色、非散射介质中的辐射平衡，$\nabla  \cdot  \mathbf{q} = 0$，源函数由方程(9.64)和(16.39)给出，$S = {I}_{b} = G/{4\pi }$，它不是方向的函数。

我们首先用流行的菱形格式来解决这个问题，即 ${\gamma }_{x} = {\gamma }_{y} = \frac{1}{2}$。由于所有单元格面面积都是 $A = L/2$，所有 $\left| {\xi }_{i}\right|  = \left| {\eta }_{i}\right|  = {0.5}$，且 ${\beta V} = \kappa {\left( L/2\right) }^{2} = {0.25\kappa }{\mathrm{L}}^{2} = {0.25}\mathrm{\;L}$，方程(16.72)变为

$$
{I}_{i,P} = \frac{\frac{1}{4}{S}_{P} + \frac{1}{2}{I}_{i,w} + \frac{1}{2}{I}_{i,s}}{\frac{1}{4} + \frac{1}{2} + \frac{1}{2}} = \frac{1}{5}\left( {{S}_{P} + 2{I}_{i,w} + 2{I}_{i,s}}\right) .
$$




<!-- Meanless: The Method of Discrete Ordinates (S ${}_{N}$ -Approximation) Chapter $\mid$-->

我们从左下角开始，对所有方向为 ${\xi }_{i} > 0$ 和 ${\eta }_{i} > 0$（即，对于 ${S}_{2}$ -近似，只有一个方向）。在接下来的方程中，单元中心强度使用两个下标。例如，${I}_{2,4}$ 表示方向 ${\widehat{\mathbf{s}}}_{2}$ 中单元格4的强度。对于面上的强度，使用三个下标，其中面索引在单元格索引之前。例如，${I}_{1,w,2}$ 表示指向方向 ${\widehat{\mathbf{s}}}_{1}$ 的单元格2西面(w)的强度。$i = 1\left\lbrack  {{\widehat{\mathbf{s}}}_{1} = {0.5}\left( {\widehat{\mathbf{i}} + \widehat{\mathbf{j}}}\right) }\right\rbrack$：对于所有单元格

$$
{I}_{1,P} = \frac{1}{5}\left( {{S}_{P} + 2{I}_{1,w,P} + 2{I}_{1,s,P}}\right) ,
$$

$$
{I}_{1,e,P} = 2{I}_{1,P} - {I}_{1,w,P},
$$

$$
{I}_{1,n,P} = 2{I}_{1,P} - {I}_{1,s,P},\;P = 1,2,3,4.
$$

从单元格1开始，我们有 ${I}_{1,w,1} = 0,{I}_{1,s,1} = {I}_{bw}$，并且

$$
{I}_{1,1} = \frac{1}{5}\left( {{S}_{1} + 2{I}_{bw}}\right) ,
$$

$$
{I}_{1,e,1} = 2{I}_{1,1} = {I}_{1,w,2},{I}_{1,n,1} = 2{I}_{1,1} - {I}_{bw} = {I}_{1,s,3};
$$

$$
{I}_{1,2} = \frac{1}{5}\left( {{S}_{2} + 2{I}_{1,w,2} + 2{I}_{1,s,2}}\right)  = \frac{1}{5}\left( {{S}_{2} + 4{I}_{1,1} + 2{I}_{bw}}\right) ,
$$

$$
{I}_{1,n,2} = 2{I}_{1,2} - {I}_{bw} = {I}_{1,s,4}
$$

$$
{I}_{1,3} = \frac{1}{5}\left( {{S}_{3} + 2{I}_{1,s,3}}\right)  = \frac{1}{5}\left( {{S}_{3} + 4{I}_{1,1} - 2{I}_{bw}}\right) ,
$$

$$
{I}_{1,e,3} = 2{I}_{1,3} = {I}_{1,w,4};
$$

$$
{I}_{1,4} = \frac{1}{5}\left( {{S}_{4} + 2{I}_{1,w,4} + 2{I}_{1,s,4}}\right) 
$$

$$
 = \frac{1}{5}\left( {{S}_{4} + 4{I}_{1,3} + 4{I}_{1,2} - 2{I}_{bw}}\right) .
$$

$i = 2\left\lbrack  {{\widehat{\mathbf{s}}}_{2} = {0.5}\left( {-\widehat{\mathbf{t}} + \widehat{\mathbf{j}}}\right) }\right\rbrack$：在一个没有对称性的问题中，我们会从右下角开始，再次扫描所有单元格。然而，在这个问题中，我们可以通过对称性立即确定强度，如下所示

$$
{I}_{2,1} = {I}_{1,2},{I}_{2,2} = {I}_{1,1},{I}_{2,3} = {I}_{1,4},{I}_{2,4} = {I}_{1,3}\text{.}
$$

$i = 3\left\lbrack  {{\widehat{\mathbf{s}}}_{3} =  - {0.5}\left( {\widehat{\mathbf{i}} + \widehat{\mathbf{j}}}\right) }\right\rbrack$：从右上角开始，对于所有单元格，我们有

$$
{I}_{3,P} = \frac{1}{5}\left( {{S}_{P} + 2{I}_{3,e,P} + 2{I}_{3,n,P}}\right) ,
$$

$$
{I}_{3,w,P} = 2{I}_{3,P} - {I}_{3,e,P},
$$

$$
{I}_{3,s,P} = 2{I}_{3,P} - {I}_{3,n,P}.
$$

从单元格4开始，${I}_{3,e,4} = {I}_{3,n,4} = 0$，我们发现

$$
{I}_{3,4} = \frac{1}{5}{S}_{4}
$$

$$
{I}_{3,s,4} = 2{I}_{3,4} = {I}_{3,n,2},{I}_{3,w,4} = 2{I}_{3,4} = {I}_{3,e,3};
$$

$$
{I}_{3,3} = \frac{1}{5}\left( {{S}_{3} + 2{I}_{3,e,3}}\right)  = \frac{1}{5}\left( {{S}_{3} + 4{I}_{3,4}}\right) ,
$$

$$
{I}_{3,s,3} = 2{I}_{3,3} = {I}_{3,n,1};
$$

$$
{I}_{3,2} = \frac{1}{5}\left( {{S}_{2} + 2{I}_{3,n,2}}\right)  = \frac{1}{5}\left( {{S}_{2} + 4{I}_{3,4}}\right) ,
$$

$$
{I}_{3,w,2} = 2{I}_{3,2} = {I}_{3,e,1};
$$

$$
{I}_{3,1} = \frac{1}{5}\left( {{S}_{1} + 2{I}_{3,e,1} + 2{I}_{3,n,1}}\right)  = \frac{1}{5}\left( {{S}_{1} + 4{I}_{3,2} + 4{I}_{3,3}}\right) .
$$

同样

$$
{I}_{3,s,1} = 2{I}_{3,1} - {I}_{3,n,1} = 2\left( {{I}_{3,1} - {I}_{3,3}}\right) ,
$$

$$
{I}_{3,s,2} = 2{I}_{3,2} - {I}_{3,n,2} = 2\left( {{I}_{3,2} - {I}_{3,4}}\right) ,
$$

这将在稍后用于从方程(16.8)计算壁面热流。

$i = 4\left\lbrack  {{\widehat{\mathbf{s}}}_{4} = {0.5}\left( {\widehat{\mathbf{i}} - \widehat{\mathbf{j}}}\right) }\right\rbrack$：再次，通过对称性立即得出

$$
{I}_{4,1} = {I}_{3,2},{I}_{4,2} = {I}_{3,1},{I}_{4,3} = {I}_{3,4},{I}_{4,4} = {I}_{3,3}\text{,}
$$

并且

$$
{I}_{4,s,1} = {I}_{3,s,2},{I}_{4,s,2} = {I}_{3,s,1}\text{.}
$$




<!-- Meanless: 590 Radiative Heat Transfer-->

总结一下，我们有

$$
{I}_{1,1} = {I}_{2,2} = \frac{1}{5}\left( {{S}_{1} + 2{I}_{bw}}\right) ,
$$

$$
{I}_{1,2} = {I}_{2,1} = \frac{1}{5}\left( {{S}_{2} + 4{I}_{1,1} + 2{I}_{bw}}\right) ,
$$

$$
{I}_{1,3} = {I}_{2,4} = \frac{1}{5}\left( {{S}_{3} + 4{I}_{1,1} - 2{I}_{bw}}\right) ,
$$

$$
{I}_{1,4} = {I}_{2,3} = \frac{1}{5}\left( {{S}_{4} + 4{I}_{1,3} + 4{I}_{1,2} - 2{I}_{bw}}\right) ,
$$

$$
{I}_{3,1} = {I}_{4,2} = \frac{1}{5}\left( {{S}_{1} + 4{I}_{3,2} + 4{I}_{3,3}}\right) ,
$$

$$
{I}_{3,2} = {I}_{4,1} = \frac{1}{5}\left( {{S}_{2} + 4{I}_{3,4}}\right) ,
$$

$$
{I}_{3,3} = {I}_{4,4} = \frac{1}{5}\left( {{S}_{3} + 4{I}_{3,4}}\right) ,
$$

$$
{I}_{3,4} = {I}_{4,3} = \frac{1}{5}{S}_{4}
$$

$$
{I}_{3,s,1} = {I}_{4,s,2} = 2\left( {{I}_{3,1} - {I}_{3,3}}\right) ,
$$

$$
{I}_{3,s,2} = {I}_{4,s,1} = 2\left( {{I}_{3,2} - {I}_{3,4}}\right) .
$$

<!-- Media -->

表16.4 例16.6的节点强度作为迭代次数的函数，由 ${I}_{bw}$ 归一化。

<table><tr><td>迭代次数</td><td>${I}_{1,1}$</td><td>${I}_{1,2}$</td><td>${I}_{1,3}$</td><td>${I}_{1,4}$</td><td>${I}_{3,1}$</td><td>${I}_{3,2}$</td><td>${I}_{3,3}$</td><td>${I}_{3,4}$</td><td>${S}_{1}$</td><td>${S}_{3}$</td></tr><tr><td colspan="11">菱形格式</td></tr><tr><td>1</td><td>0.4000</td><td>0.7200</td><td>${0.0000}^{a}$</td><td>0.1760</td><td>0.0000</td><td>0.0000</td><td>0.0000</td><td>0.0000</td><td>0.2800</td><td>0.0440</td></tr><tr><td>2</td><td>0.4560</td><td>0.8208</td><td>0.0000</td><td>0.2654</td><td>0.1191</td><td>0.0630</td><td>0.0158</td><td>0.0088</td><td>0.3647</td><td>0.0725</td></tr><tr><td>3</td><td>0.4729</td><td>0.8513</td><td>0.0000</td><td>0.2955</td><td>0.1615</td><td>0.0846</td><td>0.0261</td><td>0.0145</td><td>0.3926</td><td>0.0840</td></tr><tr><td>≥9</td><td>0.4815</td><td>0.8667</td><td>0.0037</td><td>0.3148</td><td>0.1852</td><td>0.0963</td><td>0.0333</td><td>0.0185</td><td>0.4074</td><td>0.0926</td></tr><tr><td colspan="11">阶梯格式</td></tr><tr><td>1</td><td>0.3333</td><td>0.4444</td><td>0.1111</td><td>0.1852</td><td>0.0000</td><td>0.0000</td><td>0.0000</td><td>0.0000</td><td>0.1944</td><td>0.0741</td></tr><tr><td>2</td><td>0.3981</td><td>0.5309</td><td>0.1574</td><td>0.2541</td><td>0.1001</td><td>0.0730</td><td>0.0329</td><td>0.0247</td><td>0.2755</td><td>0.1173</td></tr><tr><td>3</td><td>0.4252</td><td>0.5669</td><td>0.1808</td><td>0.2883</td><td>0.1442</td><td>0.1049</td><td>0.0521</td><td>0.0391</td><td>0.3103</td><td>0.1401</td></tr><tr><td>≥10</td><td>0.4459</td><td>0.5946</td><td>0.2027</td><td>0.3198</td><td>0.1802</td><td>0.1306</td><td>0.0721</td><td>0.0541</td><td>0.3378</td><td>0.1622</td></tr></table>

* 负值设为零。

<!-- Media -->

源函数可以很容易地从方程(16.7)和对称性中评估出来

$$
{S}_{1} = {S}_{2} = \frac{1}{4}\left( {{I}_{1,1} + {I}_{2,1} + {I}_{3,1} + {I}_{4,1}}\right) ,
$$

$$
{S}_{3} = {S}_{4} = \frac{1}{4}\left( {{I}_{1,3} + {I}_{2,3} + {I}_{3,3} + {I}_{4,3}}\right) .
$$

由于方程是线性的，可以将 ${S}_{P}$ 的关系代入上述方程，并通过求解得到的线性方程组来确定未知的 ${I}_{i,P}$。通常，如果方程数量很大（细网格和/或多方向），应使用迭代线性代数求解器。这里，我们采用一个简单的逐点迭代求解器。我们从设置所有 ${S}_{P} = 0$ 开始，找出 ${I}_{i,P}$ 的值，更新 ${S}_{P}$，重新评估 ${I}_{i,P}$，依此类推，直到达到收敛。强度的变化值（用 ${I}_{bw}$ 归一化）作为迭代的函数在表16.4中给出。三次迭代后达到约 $5\%$ 的精度，九次迭代后获得完全收敛的值（到四位有效数字）。收敛的强度用于确定在 $x = L/4$ 和 $x = {3L}/4$ 处从底壁的净辐射热流。从方程(16.8)我们有

$$
q\left( {x = {0.25L}}\right)  = q\left( {x = {0.75L}}\right)  = \pi {I}_{bw} - \mathop{\sum }\limits_{{i = 3}}^{4}{w}_{i}{I}_{i,s,1}\left| {\eta }_{i}\right|  = \pi {I}_{bw} - \frac{\pi }{2}\left( {{I}_{3,s,1} + {I}_{4,s,1}}\right) ,
$$

或

$$
\Psi  = \frac{{q}_{0.25L}}{{E}_{bw}} = \frac{{q}_{0.75L}}{{E}_{bw}} = 1 - \frac{{I}_{3,1} - {I}_{3,3} + {I}_{3,2} - {I}_{3,4}}{{I}_{bw}} = {0.7704}.
$$

为了比较，我们也将使用更简单但更稳定的阶梯差分格式来解决这个例子。那么我们从方程(16.60)和(16.61)得到

$$
{I}_{i,P} = \frac{1}{3}\left( {{S}_{P} + {I}_{i,{Ux}} + {I}_{i,{Uy}}}\right) ,
$$




<!-- Meanless: The Method of Discrete Ordinates (SN-Approximation) Chapter | 10-->

<!-- Media -->

<!-- figureText: 1.0 ${S}_{2}^{\prime }$ ${S}_{4}^{\prime }$ 口 ${S}_{2}$ 0.5 1 无量纲距离, $x/L$ 无量纲辐射热流 $\Psi  = {q}_{\mathrm{v}}/\sigma {T}_{w}^{4}$ 0.9 0.8 0.7 例 16.6 阶梯格式 0.6 菱形格式 - 精确解 [79] 0.5 ${S}_{2},{S}_{4}$ Truelove [15] ${\mathrm{S}}_{2}^{\prime },{S}_{4}^{\prime }$ Fiveland [11] 0.4 -->

<img src="https://cdn.noedgeai.com/bo_d141kgv7aajc73bva480_28.jpg?x=606&y=226&w=622&h=607&r=0"/>

图 16.8 例16.6中方形外壳底壁的无量纲热流。

<!-- Media -->

其中下标 ${Ux}$ 和 ${Uy}$ 分别表示 $x$ 和 $y$ 方向上的上游单元格。然后，遵循相同的程序，我们得到（稍微容易一些）

$$
{I}_{1,1} = \frac{1}{3}\left( {{S}_{1} + 0 + {I}_{bw}}\right)  = {I}_{2,2}
$$

$$
{I}_{1,2} = \frac{1}{3}\left( {{S}_{2} + {I}_{1,1} + {I}_{bw}}\right)  = {I}_{2,1}
$$

$$
{I}_{1,3} = \frac{1}{3}\left( {{S}_{3} + 0 + {I}_{1,1}}\right)  = {I}_{2,4}
$$

$$
{I}_{1,4} = \frac{1}{3}\left( {{S}_{4} + {I}_{1,3} + {I}_{1,2}}\right)  = {I}_{2,3}
$$

$$
{I}_{3,4} = \frac{1}{3}\left( {{S}_{4} + 0 + 0}\right)  = {I}_{4,3}
$$

$$
{I}_{3,3} = \frac{1}{3}\left( {{S}_{3} + {I}_{3,4} + 0}\right)  = {I}_{4,4}
$$

$$
{I}_{3,2} = \frac{1}{3}\left( {{S}_{2} + 0 + {I}_{3,4}}\right)  = {I}_{4,1}
$$

$$
{I}_{3,1} = \frac{1}{3}\left( {{S}_{1} + {I}_{3,2} + {I}_{3,3}}\right)  = {I}_{4,2}
$$

源函数和热流的公式保持不变，但无量纲通量在用上游单元格强度替换面强度后变为，

$$
\Psi  = 1 - \frac{{I}_{3,1} + {I}_{4,1}}{2{I}_{bw}} = {0.844}
$$

阶梯格式的迭代结果也包含在表16.4中。两种格式的无量纲热流结果显示在图16.8中，同时还有Razzaque及其同事[79]报告的精确结果，以及Truelove [15]的 ${S}_{2}$ - 和 ${S}_{4}$ -计算结果（针对更精细的网格）。Truelove的结果表明了良好坐标集的重要性，至少对于低阶近似是如此：${S}_{2}^{\prime }$ 和 ${S}_{4}^{\prime }$ 的结果是使用不遵守方程(16.13)的半矩条件的集合获得的（如Fiveland [11]在对矩形外壳的首次研究中使用的），而 ${S}_{2}$ 和 ${S}_{4}$ 的结果是使用表16.1中给出的集合获得的（对于 ${S}_{2}$ 使用非对称坐标）。毫不奇怪，作为二阶精确格式的菱形格式比仅为一阶精确的阶梯格式更准确。阶梯格式显示出 ${I}_{P}$ 和 ${S}_{P}$ 更平滑的分布，并且总是稳定的。另一方面，菱形格式在前几次迭代中为 ${I}_{1,3}$ 产生了非物理的负强度，这些负强度被设为零。

射线效应虽然存在，但在本例中并不明显，因为为了方便手动计算而使用了大单元格（过度的虚假散射在这种特殊情况下有所帮助）。当用细网格重复此例时，它们变得非常明显，这将在下一节关于有限角法的主题中进行（见例16.9）。

许多早期工作都使用二维正交网格和简单几何形状进行DOM计算，主要目的是评估各种角度求积格式的效果。在他早期的计算中，Fiveland [11]将 ${S}_{2}{}^{ - },{S}_{4}{}^{ - }$ 和 ${S}_{6}$ -近似应用于纯散射矩形介质 $\left( {\omega  = 1}\right)$，以及由冷黑壁界定的等温、非散射介质。Truelove [15]重复了其中一些结果，以证明良好坐标集的重要性，并为方形外壳中的辐射平衡提供了一些新结果。Jamaluddin和Smith [48]将 ${S}_{4}$ -近似应用于具有已知温度分布的矩形、非散射外壳。Kim和Lee研究了强各向异性散射的影响，使用了高阶近似（高达 ${S}_{16}$）[61]，以及准直辐照的影响[80]。Baek和Kim [81]研究了线性各向异性散射矩形外壳中的传导和辐射组合。他们还使用相同的方法（灰色恒定属性，此处无散射）研究了辐射在越过后台阶的可压缩湍流中的影响[82]。最后，Lu及其同事[34]研究了二维填充床中的辐射，以及传导和对流。虽然他们也假设了灰色属性，但他们允许它们在局部变化；对于散射，他们使用了大漫射球相函数，方程(11.85)。离散坐标法二维笛卡尔形式的其他应用可以在[83-86]中找到，都涉及组合模式传热。这里特别值得注意的是Selçuk和Kayakol [87]的研究，他们比较了 ${S}_{4}$ 方法与相关的离散传递法[88]（见第604页）的性能，发现这些方法的精度相当，而 ${S}_{4}$ 解所需的计算机时间少了三个数量级。


<!-- Meanless: 592 Radiative Heat Transfer-->

## 三维正交网格

刚刚推导和演示的方程，可以通过给控制体增加前后面 ${A}_{f}$ 和 ${A}_{b}$，并重写方程(16.71)，立即扩展到由均匀正交网格离散化的三维几何体：

$$
{I}_{i,P} = \frac{{R}_{i,P}}{{D}_{i,P}},\;i = 1,2,\ldots ,n;
$$

$$
{R}_{i,P} = {\beta }_{P}{V}_{P}{S}_{i,P}
$$

$$
 + \left( {\left| {\xi }_{i}\right|  - {\xi }_{i}}\right) {A}_{e}{I}_{i,e} + \left( {\left| {\xi }_{i}\right|  + {\xi }_{i}}\right) {A}_{w}{I}_{i,w} + \left( {\left| {\eta }_{i}\right|  - {\eta }_{i}}\right) {A}_{n}{I}_{i,n} + \left( {\left| {\eta }_{i}\right|  + {\eta }_{i}}\right) {A}_{s}{I}_{i,s} + \left( {\left| {\mu }_{i}\right|  - {\mu }_{i}}\right) {A}_{f}{I}_{i,f} + \left( {\left| {\mu }_{i}\right|  + {\mu }_{i}}\right) {A}_{b};
$$

$$
{D}_{i,P} = {\beta }_{P}{V}_{P}
$$

$$
 + \left( {\left| {\xi }_{i}\right|  + {\xi }_{i}}\right) {A}_{e} + \left( {\left| {\xi }_{i}\right|  - {\xi }_{i}}\right) {A}_{w} + \left( {\left| {\eta }_{i}\right|  + {\eta }_{i}}\right) {A}_{n} + \left( {\left| {\eta }_{i}\right|  - {\eta }_{i}}\right) {A}_{s} + \left( {\left| {\mu }_{i}\right|  + {\mu }_{i}}\right) {A}_{f} + \left( {\left| {\mu }_{i}\right|  - {\mu }_{i}}\right) {A}_{b}. \tag{16.73}
$$

一个三维笛卡尔外壳有八个角，从每个角必须追踪 $\frac{1}{8}N\left( {N + 2}\right)$ 个方向（覆盖一个八分圆的方向），总共有 $N\left( {N + 2}\right)$ 个坐标。方程(16.73)仅对在第一八分圆中向前传播的射线有效。Jamaluddin和Smith [89]（具有规定温度的非散射介质）、Fiveland [13]和Truelove [16]（两者都研究了Mengüç和Viskanta [90]的理想化炉膛，考虑了在辐射平衡下具有内部热生成的线性各向异性散射介质）、Park和Yoon [91]（结合传导和辐射，使用逆分析来确定给定温度分布下的恒定、灰色值$\kappa$和${\sigma }_{s}$）以及Lacroix及其同事（激光焊接过程中形成的等离子体中的辐射）[92]等人已经进行了一些这样的计算。此外，Gonçalves和Coelho [93]展示了如何在并行计算机上实现离散坐标法。Fiveland和Jessee [94]讨论了几种用于光学厚几何体的加速方案，已知离散坐标法在这些几何体中收敛非常慢（或根本不收敛）。Balsara [95]从计算机科学的角度对离散坐标法进行了截至2000年的广泛综述，强调了收敛率以及多重网格和并行实现。

## 多维非笛卡尔几何

一些研究已将离散坐标法应用于二维和三维圆柱形外壳，最近，该方法也已应用于不规则几何形状。Fiveland [10] 首次考虑了二维轴对称外壳，他为已知温度分布的圆柱形炉计算了辐射热流率。Jamaluddin和Smith [89] 处理了一个非常类似的问题，他们稍后还解决了三维圆柱形炉的情况 [96,97]。Kim和Baek [98] 研究了具有灰色、恒定属性、吸收/发射和各向同性散射介质的完全发展的非轴对称管道流。Kaplan及其同事 [99] 模拟了一个非定常乙烯扩散火焰，将烟尘和燃烧气体视为灰色和非散射，但具有空间变化。Ramamurthy及其同事 [100] 研究了辐射管中的反应、辐射流，使用了更复杂的模型来描述燃烧气体的光谱行为，Song及其同事 [101] 研究了熔融玻璃射流。所有这些都在二维、圆柱形几何中使用了 ${S}_{4}$ -法，尽管Jendoubi等人 [102] 使用了 ${S}_{14}$ -方案来评估不同的散射行为。


<!-- Meanless: The Method of Discrete Ordinates (S ${}_{N}$ -Approximation) Chapter $\mid$-->

在迄今为止提出的公式中，由于前面陈述的原因，方程(16.4)的空间离散化是使用有限体积法进行的。原则上，DOM也可以与有限差分法一起用于空间离散化。对于这两种方法，如果不调用曲线坐标变换，很难处理复杂的三维几何形状。尽管如此，Howell和Beckner [103]尝试了这样做，他们使用“嵌入边界”来模拟不规则表面，而Adams和Smith [104]则模拟了一个复杂的炉膛。他们的结果清楚地显示了射线效应：使用粗糙的坐标网格（高达${S}_{8}$）和非常精细的空间网格，他们计算的辐射通量经历了非常强烈的非物理振荡。Sakami及其同事[105,106]展示了如何在非结构化、三角形、二维网格上进行空间微分。他们通过每个单元格追溯每条射线，使用伽辽金有限元类型的方案对整个单元格进行积分。Cheong和Song [107-110]提出了一个有些类似的方法。通过仔细的空间差分，他们展示了如何将标准的离散坐标法应用于非结构化网格和不规则几何形状。该方法也由Seo和Kim [111]进一步完善。

另一种被广泛用于求解方程(16.4)的方法是有限元法(FEM)。有限元法因其对复杂几何形状和网格拓扑的兼容性以及易于扩展到更高阶方案的能力而具有吸引力。已经使用了许多不同变体的FEM。这些包括传统的伽辽金法[112-115]、最小二乘FEM [116-118]、间断伽辽金FEM [119-121]、控制体积FEM [122]以及基于偶宇称公式的FEM公式[123]，在第16.8节中描述。近年来，使用无网格方法[124-128]的趋势也日益增长，特别是所谓的自然元法[126-128]，这是一种采用有限元类型公式的无网格方法。标准的离散坐标法也已经被公式化并用于具有空间变化折射率或所谓渐变介质的介质中[129,130]。该领域最早的发展可以追溯到Lemonnier和Le Dez [129]，他们为渐变一维平板公式化并应用了DOM。

### 16.6 有限角法 (FAM)

离散坐标法，以其标准形式，存在一些严重的缺点，例如虚假散射和射线效应。虚假散射是空间离散化的表现，而射线效应是由角度离散化引起的。此外，必须满足半程矩，方程(16.13)，才能准确评估表面通量，这使得将该方法应用于不规则几何形状非常困难。也许，该方法最严重的缺点是它不能保证辐射能量的守恒。这是因为标准离散坐标法使用简单的求积进行角度离散化，尽管通常如前几节所述，对空间离散化使用有限体积法。因此，该方法演变中的一个逻辑步骤是转向完全积分或“有限”的方法，无论是在空间还是方向上。这最初是由Briggs及其同事[131]在中子输运领域提出的。Raithby及其同事[132-135]给出了辐射传热的第一个公式。Chai及其同事[136-138]提出了略有不同的方案。Raithby [139]给出了一个很好的综述。

所谓的（辐射的）有限体积法使用对有限立体角的精确积分，这类似于在空间有限体积公式中对有限体积进行的积分[见方程(16.42)]。正是由于这两个公式的这种类似性质，该方法才被称为有限体积法，而实际上，这个术语指的是一种类似于有限体积法的方法，但是在角度空间中。为了消除这种混淆——特别是因为（空间中的）有限体积法也可以与原始的离散坐标法结合使用——此后，我们将把这种离散坐标法的变体称为有限角法（FAM）。有限角法在角度空间中是完全守恒的：对于任意几何形状，可以精确满足所有全矩和半矩，并且没有辐射能量的损失。角度网格可以适应每种特殊情况，例如准直辐照[136]。


<!-- Meanless: 594 Radiative Heat Transfer-->

<!-- Media -->

<!-- figureText: $\frac{\Delta {\theta }_{i}}{2}$ ${\Omega }_{i}$ -->

<img src="https://cdn.noedgeai.com/bo_d141kgv7aajc73bva480_31.jpg?x=586&y=229&w=661&h=565&r=0"/>

图 16.9 FAM中使用的有限控制角（或立体角）的图示。

<!-- Media -->

## 通用多维公式

在有限角法中，RTE不仅在空间上对控制体积（或单元格）进行积分，如从方程(16.42)推导方程(16.48)时所示，而且还对控制角进行积分。控制角可以通过将总立体角 ${4\pi }$ 分解为任意不重叠的有限大小的立体角，或通过沿极角和方位角使用相等或不相等的分区来选择。前者将产生非结构化角度网格，而后者将产生结构化角度网格。在任何一种情况下，人们都可以想象追踪一束射线锥穿过空间体积，而不是像原始离散坐标法中那样追踪单条射线。

在第16.5节中介绍的公式中，RTE首先被分解为 $n$ 个方向RTE。然后是在有限体积或单元格上进行空间积分。这两个步骤的最终结果是方程(16.48)。这里，我们不是在方向上分解RTE，而是首先进行有限体积（在空间中）积分。遵循从方程(16.42)推导方程(16.48)的完全相同的过程，我们得到

$$
\mathop{\sum }\limits_{{f = 1}}^{{N}_{f}}\left( {\widehat{\mathbf{s}} \cdot  {\widehat{\mathbf{n}}}_{f}}\right) {I}_{f}{A}_{f} =  - {\beta }_{P}{I}_{P}{V}_{P} + {\beta }_{P}{S}_{P}{V}_{P}. \tag{16.74}
$$

我们现在将角空间分解为一组 $n$ 个不重叠的有限立体角 ${\Omega }_{i}$，使得 $\mathop{\sum }\limits_{{i = 1}}^{n}{\Omega }_{i} = {4\pi }$，其中 $n$ 现在表示立体角（或方向）的总数。每个立体角被假定为围绕一个“中心”方向 ${\widehat{\mathbf{s}}}_{i}$（类似于控制体积或单元格的单元中心），如图16.9所示。接下来，我们将方程(16.74)在有限立体角 ${\Omega }_{i}$ 上积分，得到

$$
\mathop{\sum }\limits_{{f = 1}}^{{N}_{f}}\left\lbrack  {{\int }_{{\Omega }_{i}}\widehat{\mathbf{s}}{I}_{f}{d\Omega }}\right\rbrack   \cdot  {\widehat{\mathbf{n}}}_{f}{A}_{f} =  - {\beta }_{P}{V}_{P}{\int }_{{\Omega }_{i}}{I}_{P}{d\Omega } + {\beta }_{P}{V}_{P}{\int }_{{\Omega }_{i}}{S}_{P}{d\Omega },\;i = 1,2,\ldots ,n. \tag{16.75}
$$

注意到 ${I}_{f}$ 是面 $f$ 上的面积平均强度，我们接下来假设这个面积平均强度等于方向 ${\widehat{\mathbf{s}}}_{i}$ 上的强度，我们将其表示为 ${I}_{i,f}$。使用这个假设，${\int }_{{\Omega }_{i}}\widehat{\mathbf{s}}{I}_{f}{d\Omega } = {I}_{i,f}{\int }_{{\Omega }_{i}}\widehat{\mathbf{s}}{d\Omega } = {I}_{i,f}{\mathbf{S}}_{i}$，其中

$$
{\mathbf{S}}_{i} = {\int }_{{\Omega }_{i}}\widehat{\mathbf{s}}{d\Omega } = {\int }_{\Delta {\theta }_{i}}{\int }_{\Delta {\psi }_{i}}\widehat{\mathbf{s}}\sin {\theta d\theta d\psi }. \tag{16.76}
$$




<!-- Meanless: The Method of Discrete Ordinates ( ${S}_{N}$ -Approximation) Chapter-->

在方程(16.76)的最后一部分，使用了极坐标中立体角的定义[方程(1.29)]。将方程(16.76)代入方程(16.75)得到

$$
\mathop{\sum }\limits_{{f = 1}}^{{N}_{f}}{I}_{i,f}\left( {{\mathbf{S}}_{i} \cdot  {\widehat{\mathbf{n}}}_{f}}\right) {A}_{f} =  - {\beta }_{P}{V}_{P}{\int }_{{\Omega }_{i}}{I}_{P}{d\Omega } + {\beta }_{P}{V}_{P}{\int }_{{\Omega }_{i}}{S}_{P}{d\Omega },\;i = 1,2,\ldots ,n. \tag{16.77}
$$

最后，将平均的定义扩展到角度平均，我们得到

$$
\mathop{\sum }\limits_{{f = 1}}^{{N}_{f}}{I}_{i,f}\left( {{\mathbf{S}}_{i} \cdot  {\widehat{\mathbf{n}}}_{f}}\right) {A}_{f} =  - {\beta }_{P}{V}_{P}{I}_{i,P}{\Omega }_{i} + {\beta }_{P}{V}_{P}{S}_{i,P}{\Omega }_{i},\;i = 1,2,\ldots ,n, \tag{16.78}
$$

其中 ${I}_{i,P}$ 是沿方向 ${\widehat{\mathbf{s}}}_{i}$ 和单元格 $P$ 的平均（空间和角度）强度。同样，${S}_{i,P}$ 是平均源项，表示为

$$
{S}_{i,P} = \left( {1 - {\omega }_{P}}\right) {I}_{bP} + \frac{{\omega }_{P}}{4\pi }\mathop{\sum }\limits_{{j = 1}}^{n}{I}_{j,P}{\bar{\Phi }}_{ij}, \tag{16.79a}
$$

$$
{\bar{\Phi }}_{ij} = \frac{1}{{\Omega }_{i}}{\int }_{{\Omega }_{i}}{\int }_{{\Omega }_{j}}\bar{\Phi }\left( {{\widehat{\mathbf{s}}}^{\prime },\widehat{\mathbf{s}}}\right) d{\Omega }^{\prime }{d\Omega }. \tag{16.79b}
$$

在FAM中，方程(16.76)中的积分是解析进行的。角度空间通常由极角和方位角描述，方程(1.27)描述了方向向量 ${\widehat{\mathbf{s}}}_{i}$。此外，$\theta$ 和 $\psi$ 上的积分限通常取为 $\left\lbrack  {{\theta }_{i} - \Delta {\theta }_{i}/2,{\theta }_{i} + \Delta {\theta }_{i}/2}\right\rbrack$ 和 $\left\lbrack  {{\psi }_{i} - \Delta {\psi }_{i}/2,{\psi }_{i} + \Delta {\psi }_{i}/2}\right\rbrack$，其中 ${\theta }_{i}$ 和 ${\psi }_{i}$ 是围绕 ${\widehat{\mathbf{s}}}_{i}$ 所张立体角的“中心”极角和方位角，如图16.9所示。在这些规定下，将方程(1.27)代入方程(16.76)得到

$$
{\mathbf{S}}_{i} = {\int }_{{\Omega }_{i}}\widehat{\mathbf{s}}{d\Omega } = {\int }_{\Delta {\theta }_{i} - \Delta {\theta }_{i}/2}^{\Delta {\theta }_{i} + \Delta {\theta }_{i}/2}{\int }_{\Delta {\psi }_{i} - \Delta {\psi }_{i}/2}^{{\psi }_{i} + \Delta {\psi }_{i}/2}\left( {\sin \theta \cos \psi \widehat{\mathbf{r}} + \sin \theta \sin \psi \widehat{\mathbf{j}} + \cos \theta \widehat{\mathbf{k}}}\right) \sin {\theta d\theta d\psi }, \tag{16.80}
$$

解析积分后，得到

$$
{\mathbf{S}}_{i} = \left\lbrack  {\left( {\Delta {\theta }_{i} - \cos 2{\theta }_{i}\sin \Delta {\theta }_{i}}\right) \cos {\psi }_{i}\sin \left( \frac{\Delta {\psi }_{i}}{2}\right) }\right\rbrack  \widehat{\mathbf{i}}
$$

$$
 + \left\lbrack  {\left( {\Delta {\theta }_{i} - \cos 2{\theta }_{i}\sin \Delta {\theta }_{i}}\right) \sin {\psi }_{i}\sin \left( \frac{\Delta {\psi }_{i}}{2}\right) }\right\rbrack  \widehat{\mathbf{j}} + \sin 2{\theta }_{i}\sin \Delta {\theta }_{i}\left( \frac{\Delta {\psi }_{i}}{2}\right) \widehat{\mathbf{k}}. \tag{16.81}
$$

剩下要做的就是将面的强度 ${I}_{i,f}$ 与单元中心的强度 ${I}_{i,P}$ 联系起来。这属于空间离散格式的范畴，其过程和问题与前面为标准离散坐标法描述的相同。Raithby及其同事[132]在FAM的背景下开发了高阶空间精度的格式。然而，这样的高阶格式只能用于规则的结构化网格。即使对于这样的网格，它们也需要大量的解析和计算开销。鉴于Chai及其同事[46]讨论的稳定性考虑，通常更倾向于使用简单的阶梯格式。使用由方程(16.49)给出的阶梯格式，方程(16.78)可以写成与方程(16.51)完全相同的形式：

$$
\left\lbrack  {\mathop{\sum }\limits_{{f = 1}}^{{N}_{f,P}}\left( {\frac{\left| \left( {{\mathbf{S}}_{i} \cdot  {\widehat{\mathbf{n}}}_{f}}\right) \right|  + \left( {{\mathbf{S}}_{i} \cdot  {\widehat{\mathbf{n}}}_{f}}\right) }{2}{A}_{f}}\right)  + {\beta }_{P}{V}_{P}{\Omega }_{i}}\right\rbrack  {I}_{i,P} - \left\lbrack  {\mathop{\sum }\limits_{{f = 1}}^{{N}_{f,P}}\left( {\frac{\left| \left( {{\mathbf{S}}_{i} \cdot  {\widehat{\mathbf{n}}}_{f}}\right) \right|  - \left( {{\mathbf{S}}_{i} \cdot  {\widehat{\mathbf{n}}}_{f}}\right) }{2}{A}_{f}}\right) }\right\rbrack  {I}_{i,N\left( f\right) } = {\beta }_{P}{S}_{i,P}{V}_{P}{\Omega }_{i},
$$

$$
i = 1,2,\ldots ,n\text{.} \tag{16.82}
$$

与方程(16.49)的情况一样，方程(16.82)对任意形状和大小的单元格都有效，只要单元格是凸的并且由平面界面界定。它可以用于结构化和非结构化网格计算。对于结构化的二维正交笛卡尔网格，遵循方程(16.73)的推导，方程(16.82)简化为


<!-- Meanless: 596 Radiative Heat Transfer-->

$$
{I}_{i,P} = \frac{{R}_{i,P}}{{D}_{i,P}},\;i = 1,2,\ldots ,n;
$$

$$
{R}_{i,P} = {\beta }_{P}{\Delta x\Delta y}{\Omega }_{i}{S}_{i,P}
$$

$$
 + \frac{\left| {{\mathbf{S}}_{i} \cdot  {\widehat{\mathbf{n}}}_{w}}\right|  - \left| {{\mathbf{S}}_{i} \cdot  {\widehat{\mathbf{n}}}_{w}}\right| }{2}{I}_{i,W}{\Delta y} + \frac{\left| {{\mathbf{S}}_{i} \cdot  {\widehat{\mathbf{n}}}_{e}}\right|  - \left| {{\mathbf{S}}_{i} \cdot  {\widehat{\mathbf{n}}}_{e}}\right| }{2}{I}_{i,E}{\Delta y} + \frac{\left| {{\mathbf{S}}_{i} \cdot  {\widehat{\mathbf{n}}}_{n}}\right|  - \left| {{\mathbf{S}}_{i} \cdot  {\widehat{\mathbf{n}}}_{n}}\right| }{2}{I}_{i,N}{\Delta x} + \frac{\left| {{\mathbf{S}}_{i} \cdot  {\widehat{\mathbf{n}}}_{s}}\right|  - \left| {{\mathbf{S}}_{i} \cdot  {\widehat{\mathbf{n}}}_{s}}\right| }{2}{I}_{i,S}{\Delta x};
$$

$$
{D}_{i,P} = {\beta }_{P}{\Delta x\Delta y}{\Omega }_{i}
$$

$$
 + \frac{\left| {{\mathbf{S}}_{i} \cdot  {\widehat{\mathbf{n}}}_{w}}\right|  + \left| {{\mathbf{S}}_{i} \cdot  {\widehat{\mathbf{n}}}_{w}}\right| }{2}{\Delta y} + \frac{\left| {{\mathbf{S}}_{i} \cdot  {\widehat{\mathbf{n}}}_{e}}\right|  + \left| {{\mathbf{S}}_{i} \cdot  {\widehat{\mathbf{n}}}_{e}}\right| }{2}{\Delta y} + \frac{\left| {{\mathbf{S}}_{i} \cdot  {\widehat{\mathbf{n}}}_{n}}\right|  + \left| {{\mathbf{S}}_{i} \cdot  {\widehat{\mathbf{n}}}_{n}}\right| }{2}{\Delta x} + \frac{\left| {{\mathbf{S}}_{i} \cdot  {\widehat{\mathbf{n}}}_{s}}\right|  + \left| {{\mathbf{S}}_{i} \cdot  {\widehat{\mathbf{n}}}_{s}}\right| }{2}{\Delta x}, \tag{16.83}
$$

其中 ${\widehat{\mathbf{n}}}_{w} =  - \widehat{\mathbf{1}},{\widehat{\mathbf{n}}}_{e} = \widehat{\mathbf{1}},{\widehat{\mathbf{n}}}_{n} = \widehat{\mathbf{J}}$ 和 ${\widehat{\mathbf{n}}}_{s} =  - \widehat{\mathbf{J}}$。再次，根据传播射线的方向向量，相邻单元格对分子和分母的四个贡献项中的两个将消失。

边界条件的推导方式与标准离散坐标法类似，但对于漫射发射和反射表面，为了确保未与立体角 ${\Omega }_{i}$ 对齐的表面的辐射能量守恒，进行能量平衡是有利的。将方程(16.2)乘以 $\widehat{\mathbf{n}} \cdot  \widehat{\mathbf{s}}$ 并对所有出射方向积分，得到表面辐射出射度的表达式

$$
{J}_{w} = {\int }_{\widehat{\mathbf{n}} \cdot  \widehat{\mathbf{s}} > 0}I\widehat{\mathbf{n}} \cdot  \widehat{\mathbf{s}}{d\Omega } = {\int }_{\widehat{\mathbf{n}} \cdot  \widehat{\mathbf{s}} > 0}{\epsilon }_{w}{I}_{bw}\widehat{\mathbf{n}} \cdot  \widehat{\mathbf{s}}{d\Omega } + \left( {1 - {\epsilon }_{w}}\right) {\int }_{\widehat{\mathbf{n}} \cdot  \widehat{\mathbf{s}} < 0}I\left| {\widehat{\mathbf{n}} \cdot  \widehat{\mathbf{s}}}\right| {d\Omega }. \tag{16.84}
$$

我们现在考虑一个边界面，$f = w$。注意到边界面处的表面法线与方向无关，方程(16.84)可以重写为

$$
\left( {{\int }_{{\widehat{\mathbf{n}}}_{w} \cdot  \widehat{\mathbf{s}} > 0}I\widehat{\mathbf{s}}{d\Omega }}\right)  \cdot  {\widehat{\mathbf{n}}}_{w} = {\epsilon }_{w}{I}_{iw}\left( {{\int }_{{\widehat{\mathbf{n}}}_{w} \cdot  \widehat{\mathbf{s}} > 0}\widehat{\mathbf{s}}{d\Omega }}\right)  \cdot  {\widehat{\mathbf{n}}}_{w} - \left( {1 - {\epsilon }_{w}}\right) \left( {{\int }_{{\widehat{\mathbf{n}}}_{w} \cdot  \widehat{\mathbf{s}} < 0}I\widehat{\mathbf{s}}{d\Omega }}\right)  \cdot  {\widehat{\mathbf{n}}}_{w}. \tag{16.85}
$$

遵循前面为FAM描述的程序，方程(16.85)简化为

$$
{I}_{i,w,\text{ out }}\mathop{\sum }\limits_{{i,\text{ out }}}{\mathbf{S}}_{i} \cdot  {\widehat{\mathbf{n}}}_{w} = {\epsilon }_{w}{I}_{bw}\mathop{\sum }\limits_{{i,\text{ out }}}{\mathbf{S}}_{i} \cdot  {\widehat{\mathbf{n}}}_{w} - \left( {1 - {\epsilon }_{w}}\right) {I}_{i,w}\mathop{\sum }\limits_{{i,\text{ in }}}{\mathbf{S}}_{i} \cdot  {\widehat{\mathbf{n}}}_{w} \tag{16.86}
$$

或

$$
{I}_{i,w,\text{ out }} = {\epsilon }_{w}{I}_{bw} + \left( {1 - {\epsilon }_{w}}\right) \mathop{\sum }\limits_{{i,\text{ in }}}{I}_{i,w}\left| {{\mathbf{S}}_{i} \cdot  {\widehat{\mathbf{n}}}_{w}}\right| \sqrt{\mathop{\sum }\limits_{{i,\text{ out }}}\left( {{\mathbf{S}}_{i} \cdot  {\widehat{\mathbf{n}}}_{w}}\right) }, \tag{16.87}
$$

其中 ${I}_{i,w,{out}}$ 是离开边界面 $w$ 的漫射强度（对于所有出射方向 ${\Omega }_{i}$ 且 ${\mathbf{S}}_{i} \cdot  {\widehat{\mathbf{n}}}_{w} > 0$ 都相同），${\widehat{\mathbf{n}}}_{w}$ 是在 $w$ 处的单位表面法线，指向边界外（但进入相邻的体积元 $P$）。${I}_{i,w}$ 是离开相邻体积元 $P$ 进入边界面 $w$ 的强度。使用阶梯格式，我们可以设置 ${I}_{i,w} = {I}_{i,P}$ 对于 ${\mathbf{S}}_{i} \cdot  {\widehat{\mathbf{n}}}_{w} < 0$。

一旦所有内部强度 ${I}_{i,P}$ 和边界强度 ${I}_{i,w}$ 都被确定，入射辐射和辐射通量的内部值可以从下式找到

$$
{G}_{P} = \mathop{\sum }\limits_{i}{I}_{i,P}{\Omega }_{i},\;{\mathbf{q}}_{w} = \mathop{\sum }\limits_{i}{I}_{i,P}{\mathbf{S}}_{i}, \tag{16.88}
$$

而壁面通量由下式给出

$$
{q}_{w} = {\epsilon }_{w}\left( {{E}_{bw} - {H}_{w}}\right)  = {\epsilon }_{w}\left( {{I}_{bw}\mathop{\sum }\limits_{{i,\text{ out }}}{\mathbf{S}}_{i} \cdot  {\widehat{\mathbf{n}}}_{w} - \mathop{\sum }\limits_{{i,\text{ in }}}{I}_{i,w}\left| {{\mathbf{S}}_{i} \cdot  {\widehat{\mathbf{n}}}_{w}}\right| }\right) . \tag{16.89}
$$

注意，对于任意方向的表面，$\left| {{\mathbf{S}}_{i} \cdot  {\widehat{\mathbf{n}}}_{w}}\right|$ 的和可能不等于 $\pi$（对于入射或出射方向）；因此，为了保持一致性，方程(16.89)最右侧给出的 ${E}_{bro}$ 的有限体积表示是首选。在圆柱坐标中，FAM的公式化需要与DOM类似的修改，细节可以在[133,140,141]中找到。


<!-- Meanless: The Method of Discrete Ordinates (S ${}_{N}$ -Approximation) Chapter $\mid$-->

## 角网格的选择

在FAM中，立体角可以以多种替代方式选择，并且没有限制。对于二维平面计算，网格生成器通常在 $x - y$ 平面上生成网格。由于方位角 $\psi$ 是在 $x - y$ 平面中测量的（见图16.9），因此很自然地，角度离散化必须在 $\psi$ 中进行——无论是均匀还是非均匀的——以捕捉来自 $x - y$ 平面中所有方向的信息。另一方面，极角不需要离散化，因为 $z$ 方向延伸到 $\pm  \infty$。因此，对于二维平面几何的计算，使用一个跨越从0到 $\pi$ 的单一 $\theta$ 是可取的。任何其他选择都是多余的且计算上是浪费的。在三维中，${N}_{\theta }$ 和 ${N}_{\psi }$ 的确切选择取决于问题——特别是计算域的纵横比和边界条件。通常最简单的方法是均匀地划分 $\theta$ 和 $\psi$，这样我们就有一个均匀的 ${N}_{\theta } \times  {N}_{\psi }$ 角度网格。当然，这并不会导致相等的立体角。另一种划分角度空间的常用方法是选择相等的立体角，这基本上相当于选择导致其在单位球上的交点所包围的面积相等的 $\theta$ 和 $\psi$ 值。在二维轴对称计算中，使用 $r - z$ 平面中的二维网格，极角是从 $r - z$ 平面上的 $z$ 轴测量的，而方位角是从 $r$ 轴向纸面内测量的。由于在这种情况下没有前后对称性（如二维平面中），两个角度都需要解析。即使存在方位对称性，射线也可能从圆柱曲面上的一个点传播到曲面上的另一点，但具有不同的边界条件，而无需在 $r - z$ 平面上行进。这就需要解析方位角。通常的做法是在 $\psi$ 中至少使用4个细分。早期的研究，如[141]，表明除非介质光学非常厚，否则仅在 $\psi$ 中使用2个细分的结果会差很多。在许多商业程序中，二维轴对称情况通过将计算域视为一个薄的三维楔形来处理，该楔形在方位角上跨越整个楔形，并使用镜面反射边界条件或反射方向在楔形两个侧面上的旋转不变映射，如Cai等人[142]所讨论。

例16.7. 使用有限角法重复例16.1，将上半球和下半球作为立体角范围。

## 解

控制方程如前所述，

$$
\mu \frac{dI}{d\tau } + I = \left( {1 - \omega }\right) {I}_{b} + \frac{\omega }{4\pi }\left( {G + {A}_{1}{q\mu }}\right) .
$$

如果我们想以类似于例16.1的方式应用有限角法，即为每个立体角范围获得一个微分方程，那么我们只需要将控制方程在这些立体角上积分，而不是在体积上积分。假设在上半球有一个恒定的强度 ${I}^{ + }$，在下半球有一个恒定的强度 ${I}^{ - }$，我们得到 ${I}_{b} = G/{4\pi }$

$$
\text{上半球：}{\int }_{0}^{2\pi }{\int }_{0}^{1}\left\lbrack  {\mu \frac{d{I}^{ + }}{d\tau } + {I}^{ + } = \frac{1}{4\pi }\left( {G + {A}_{1}{\omega \mu q}}\right) }\right\rbrack  {d\mu d\psi }\text{,}
$$

$$
\text{下半球：}{\int }_{0}^{2\pi }{\int }_{-1}^{0}\left\lbrack  {\mu \frac{d{I}^{ - }}{d\tau } + {I}^{ - } = \frac{1}{4\pi }\left( {G + {A}_{1}{\omega \mu q}}\right) }\right\rbrack  {d\mu d\psi }\text{,}
$$

或

$$
\pi \frac{d{I}^{ + }}{d\tau } + {2\pi }{I}^{ + } = \frac{1}{2}G + \frac{1}{4}{A}_{1}{\omega q},
$$

$$
 - \pi \frac{d{I}^{ - }}{d\tau } + {2\pi }{I}^{ - } = \frac{1}{2}G - \frac{1}{4}{A}_{1}{\omega q}.
$$

从热流和入射辐射的定义，我们再次得到

$$
\left\{  \begin{array}{l} G \\  q \end{array}\right\}   = {\int }_{0}^{2\pi }{\int }_{-1}^{+1}\left\{  \begin{array}{l} 1 \\  \mu  \end{array}\right\}  {Id\mu d\psi } = \left\{  \begin{matrix} {2\pi }\left( {{I}^{ + } + {I}^{ - }}\right) \\  \pi \left( {{I}^{ + } - {I}^{ - }}\right)  \end{matrix}\right\}  ,
$$




<!-- Meanless: 598 Radiative Heat Transfer-->

如例16.1所示。因此，将上半球和下半球的方程相加和相减，我们得到

$$
\frac{dq}{d\tau } + G = G\;\text{ 或 }\;\frac{dq}{d\tau } = 0,
$$

$$
\frac{1}{2}\frac{dG}{d\tau } + {2q} = \frac{1}{2}{A}_{1}{\omega q}\;\text{ 或 }\;\frac{dG}{d\tau } =  - \left( {4 - {A}_{1}\omega }\right) q.
$$

对于边界条件，方程(16.87)，我们首先需要计算 ${\mathbf{S}}_{i}$：

$$
{\mathbf{S}}_{1} = {\int }_{0}^{2\pi }{\int }_{0}^{\pi /2}\left( {\sin \theta \cos \psi \widehat{\mathbf{r}} + \sin \theta \sin \psi \widehat{\mathbf{j}} + \cos \theta \widehat{\mathbf{k}}}\right) \sin {\theta d\theta d\psi } = \pi \widehat{\mathbf{k}},
$$

$$
{\mathbf{S}}_{2} = {\int }_{0}^{2\pi }{\int }_{\pi /2}^{\pi }\left( {\sin \theta \cos \psi \widehat{\mathbf{r}} + \sin \theta \sin \psi \widehat{\mathbf{j}} + \cos \theta \widehat{\mathbf{k}}}\right) \sin {\theta d\theta d\psi } =  - \pi \widehat{\mathbf{k}}.
$$

对于底部边界，我们有 $\widehat{\mathbf{n}} = \widehat{\mathbf{k}}$ 和 ${\mathbf{S}}_{1} \cdot  \widehat{\mathbf{n}} =  - {\mathbf{S}}_{2} \cdot  \widehat{\mathbf{n}} = \pi$，所以在

$$
\tau  = 0 : \;{I}^{ + } = {\epsilon }_{1}{I}_{b1} + \left( {1 - {\epsilon }_{1}}\right) {I}^{ - }.
$$

减去 $\left( {1 - {\epsilon }_{1}}\right) {I}^{ + }$，使用 $q$ 的定义，并除以 ${\epsilon }_{1}$ 得到

$$
\tau  = 0 : \;{I}^{ + } = {J}_{1}/\pi  = {I}_{b1} - \frac{1 - {\epsilon }_{1}}{{\epsilon }_{1}\pi }q,
$$

这当然与例16.1（以及任何漫射表面）相同。类似地，在顶部墙壁

$$
\tau  = {\tau }_{L} : \;{I}^{ - } = {J}_{2}/\pi .
$$

然后立即从例16.1中找到解（设 $1/{\mu }_{1}^{2} = 4$）

$$
\Psi  = \frac{q}{{J}_{1} - {J}_{2}} = \frac{1}{1 + \left( {1 + \frac{{A}_{1}\omega }{4}}\right) {\tau }_{L}},
$$

这与非对称 ${S}_{2}$ -近似的答案相同。更重要的是，本例中的分析表明，Schuster-Schwarzschild（或双通量）近似只是最低级的有限角法。

最后，在本例中，$\theta$ 方向被离散化，而 $\psi$ 则从0到 ${2\pi }$ 跨越整个范围。这一选择是为了使用与例16.1中相同的框架。另一种方法是考虑板到板的方向为 $x$ 或 $y$ 方向（而不是 $z$），并对 $\psi$ 进行离散化（在 $\theta$ 中只有一个角度跨越0到 $\pi$），如下一个例子所示。

例16.8. 使用有限角法重复例16.6，将总立体角分成四个相等的范围。

## 解

与例16.6一样，并遵循先前为二维平面计算提供的指导方针，我们将 $z$ 轴垂直于图16.7的纸面，极角 $\theta$ 从中测量。由于是二维的，最好将每个 ${\Omega }_{i}$ 分配给整个极角范围，以及方位角范围的四分之一。因此，按象限划分我们选择

$$
{\Omega }_{1} : \;0 \leq  \psi  < \frac{\pi }{2},\;0 \leq  \theta  \leq  \pi ,
$$

$$
{\Omega }_{2} : \;\frac{\pi }{2} \leq  \psi  < \pi ,\;0 \leq  \theta  \leq  \pi ,
$$

$$
{\Omega }_{3} : \;\pi  \leq  \psi  < \frac{3\pi }{2},\;0 \leq  \theta  \leq  \pi ,
$$

$$
{\Omega }_{4} : \;\frac{3\pi }{2} \leq  \psi  < {2\pi },\;0 \leq  \theta  \leq  \pi .
$$

立体角向量 ${\mathbf{S}}_{i}$ 通过 $\widehat{\mathbf{s}} = \sin \theta \cos \psi \widehat{\mathbf{t}} + \sin \theta \sin \psi \widehat{\mathbf{j}} + \cos \theta \widehat{\mathbf{k}}$ 得到

$$
{\mathbf{S}}_{i} = {\int }_{\Delta {\psi }_{i}}{\int }_{0}^{\pi }\widehat{\mathbf{s}}\sin {\theta d\theta d\psi } = {\int }_{\Delta {\psi }_{i}}\left\lbrack  {\frac{\pi }{2}\cos \psi \widehat{\mathbf{t}} + \frac{\pi }{2}\sin \psi \widehat{\mathbf{y}} + 0\widehat{\mathbf{k}}}\right\rbrack  {d\psi } = {\left. \left( \frac{\pi }{2}\sin \psi \widehat{\mathbf{t}} - \frac{\pi }{2}\cos \psi \widehat{\mathbf{y}}\right) \right| }_{\Delta {\psi }_{i}},
$$




<!-- Meanless: The Method of Discrete Ordinates ( ${S}_{N}$ -Approximation) Chapter $\mid$-->

或

$$
{\mathbf{S}}_{1} = \frac{\pi }{2}\left( {\widehat{\mathbf{i}} + \widehat{\mathbf{j}}}\right) ,\;{\mathbf{S}}_{2} =  - \frac{\pi }{2}\left( {\widehat{\mathbf{i}} - \widehat{\mathbf{j}}}\right) ,\;{\mathbf{S}}_{3} =  - \frac{\pi }{2}\left( {\widehat{\mathbf{i}} + \widehat{\mathbf{j}}}\right) ,\;{\mathbf{S}}_{4} = \frac{\pi }{2}\left( {\widehat{\mathbf{i}} - \widehat{\mathbf{j}}}\right) ,
$$

这与例16.6中的方向相同，除了因子 $\pi$（它给出了 ${\Omega }_{i}$ 的总立体角），以及例16.6中的 ${\widehat{\mathbf{s}}}_{i}$ 是 $x-y$ 平面上的投影，而本例中的 ${\mathbf{S}}_{i}$ 实际上位于 $x-y$ 平面上。推导 ${\mathbf{S}}_{i}$ 的另一种方法是在方程(16.81)中代入 $\Delta {\theta }_{i} = \pi, {\theta }_{i} = \pi /2$ 和 $\Delta {\psi }_{i} = \pi /2$，得到

$$
{\mathbf{S}}_{i} = \pi \frac{1}{\sqrt{2}}\left( {\cos {\psi }_{i}\widehat{\mathbf{r}} + \sin {\psi }_{i}\widehat{\mathbf{j}}}\right) .
$$

使用 ${\psi }_{1} = \pi /4, {\psi }_{2} = {3\pi }/4, {\psi }_{3} = {5\pi }/4$ 和 ${\psi }_{4} = {7\pi }/4$ 会得到与前面相同的结果。现在，设 ${A}_{k} = \frac{1}{2}L, {\beta V} =$ $\kappa {\left( \frac{1}{2}L\right) }^{2} = \frac{1}{4}L$ 和 ${\Omega }_{i} = \pi$，我们从方程(16.83)得到

$$
{\Omega }_{1}\left( {i = 1}\right)  : \;{I}_{1,P} = \frac{\frac{1}{4}{L\pi }{S}_{1,P} + {I}_{1,w,P}\left| {{\mathbf{S}}_{1} \cdot  \left( {-\widehat{\mathbf{I}}}\right) }\right| \frac{1}{2}L + {I}_{1,s,P}\left| {{\mathbf{S}}_{1} \cdot  \left( {-\widehat{\mathbf{J}}}\right) }\right| \frac{1}{2}L}{\frac{1}{4}{L\pi } + {\mathbf{S}}_{1} \cdot  \widehat{\mathbf{I}}\frac{1}{2}L + {\mathbf{S}}_{1} \cdot  \widehat{\mathbf{J}}\frac{1}{2}L},P = 1,2,3,4,
$$

其他三个方向也类似。计算点积并简化后，可以写出对所有节点和所有方向的通用形式

$$
{I}_{i,P} = \frac{1}{3}\left( {{S}_{i,P} + {I}_{i,{xf},P} + {I}_{i,{yf},P}}\right) ,\;P,i = 1,2,3,4,
$$

其中 ${I}_{i,{xf},P}$ 是穿过 $x =$ 常数面进入体积 $P$ 的 ${\mathbf{S}}_{i}$ 方向上的强度，${I}_{i,{yf},P}$ 也类似。因此

$$
{I}_{1,1} = \frac{1}{3}\left( {{S}_{1,1} + {I}_{1,w,1} + {I}_{1,s,1}}\right)  = \frac{1}{3}\left( {{S}_{1,1} + 0 + {I}_{bw}}\right) ,
$$

$$
{I}_{1,2} = \frac{1}{3}\left( {{S}_{1,2} + {I}_{1,w,2} + {I}_{1,s,2}}\right)  = \frac{1}{3}\left( {{S}_{1,2} + {I}_{1,P} + {I}_{bw}}\right) ,
$$

...,

这与 ${S}_{2}$ -近似结合阶梯格式的结果完全相同。这是可以预期的，因为（i）这里应用的有限体积法（在空间上）使用了阶梯格式，并且（ii）对于这个非常简单的情况，${S}_{2}$ -近似满足所有半矩。有限角法的优势在于它能守恒辐射能。与例16.6一样，由于空间网格粗糙，射线效应不明显。

附录F包含了由Chai及其同事[136-138]开发的程序FVM2D.f，该程序求解任意矩形外壳的方程(16.83)、(16.87)和(16.89)。

例16.9. 使用程序FVM2D.f重复例16.8，以允许精细的网格和坐标分辨率。特别地，考虑只有底部条带的一部分，$- {0.1} < x <  + {0.1}$，被加热的情况。

## 解

与例16.6一样，我们将 $z$ 轴垂直于纸面。我们将使用 $N \times  N$ 个单元格在 $x$ -y平面上，并将总立体角离散化为 $1 \times  M$ 个子角，由于问题的二维性质，我们限制自己只有一个极角。FVM2D.f对于光学中等情况 ${\tau }_{L} = {\kappa L} = 1$ 时，顶部表面的辐照结果如图16.10所示。答案与蒙特卡罗法（见第20章）的精确结果进行了比较，也与上一章描述的 ${P}_{1}$ -和 ${P}_{3}$ -近似的结果进行了比较。注意，对于具有各向同性散射的灰色介质的辐射平衡，解仅取决于消光系数（无论散射或吸收/再发射的比例如何），本例与图15.12中呈现的情况相同（除了这里考虑了部分加热的条带）。图16.10a显示了对于固定的 ${10} \times  {10}$ 空间网格，方向离散化的影响。观察到，当我们从 $1 \times  4$ 方向增加到 $1 \times  8$ 方向时，精度提高，但对于更精细的方向分辨率，精度会恶化。射线效应在图16.10b中对于固定的 $1 \times  8$ 方向离散化更为明显：随着空间网格变得更精细，越来越多的远离加热壁的单元格接收不到沿所选坐标的直接辐射。显然，对于一个（相对粗糙的）$1 \times  8$ 方向离散化，一个（粗糙的）${10} \times  {10}$ 空间网格给出了最准确的答案，对于更精细的网格，精度迅速下降。注意，在有限角法中，射线效应不如标准离散坐标法（其离散方向与平均立体角相对）那么明显，因此预计后者在当前情况下表现更差。

相比之下，${P}_{N}$ 的答案相对准确，因为它们不受射线效应的影响，并且它们的答案与空间离散化无关（超出了为避免明显截断误差所需的最小 $N \times  N$ 离散化）。


<!-- Meanless: 600 Radiative Heat Transfer-->

<!-- Media -->

<!-- figureText: 0.14 0.14 ${10} \times  {10}$ FAM $1 \times  8$ ${20} \times  {20}$ FAM ${\tau }_{L} = 1$ 0.12 ${40} \times  {40}$ FAM 1×8 ${80} \times  {80}$ FAM 1×8 ${P}_{3}$ 0.10 ${P}_{1}$ 蒙特卡罗 0.08 0.06 0.04 0.02 -0.4 -0.2 0 0.2 0.4 (b) ${10} \times  {10}$ FAM 1X4 ${10} \times  {10}$ FAM 1×8 ${\tau }_{L} = 1$ 0.12 10X10 FAM 1×16 ${10} \times  {10}$ FAM 1X32 0.10 ${P}_{3}$ 蒙特卡罗 $H/\sigma {T}_{w}^{ * }$ 0.08 0.06 0.04 0.02 -0.4 -0.2 0.2 0.4 (a) -->

<img src="https://cdn.noedgeai.com/bo_d141kgv7aajc73bva480_37.jpg?x=212&y=223&w=1411&h=668&r=0"/>

图 16.10 方形外壳中底部有加热带的顶部表面辐照：(a)方向离散化的影响，(b)空间离散化的影响。

<!-- Media -->

在这个特定的例子中，${P}_{1}$ 实际上比 ${P}_{3}$ 表现得更好。这并不意味着 ${P}_{N}$ -近似优于有限角法（或者 ${P}_{1}$ 比 ${P}_{3}$ 更好），从图15.12中给出的附加结果以及[143]中对这个问题和其他几个问题给出的 ${P}_{N}$ 和FAM的直接比较可以看出。

Murthy和Mathur [69,144,145]指出，一般来说，方程(16.83)、(16.87)和(16.89)由于立体角悬垂而产生误差。例如，对于单元格$\bar{P}$的面$f$，部分具有${\mathbf{S}}_{i} \cdot  {\widehat{\mathbf{n}}}_{f} > 0$（指向体积元外）的立体角${\Omega }_{i}$实际上可能与单元格重叠。同样，立体角边界不太可能在各处都与实体边界完美对齐。他们通过像素化，即通过将${\Omega }_{i}$分解成更小的部分来确定重叠分数，从而提高了方法的准确性。此外，注意到标准的逐行迭代方法在光学厚的情况下导致收敛速度慢得无法接受，他们引入了一种新方案，该方案同时更新一个单元格内的所有方向强度，从而使收敛率基本上与光学厚度无关[69,70]。Hassanzadeh及其同事[146]也开发了一种方法，通过以平均强度$G/{4\pi }$而不是所有方向强度进行迭代来加速光学厚介质的收敛。还提出了对该方法的其他几项改进。Kim和Huh [147]注意到，大多数研究人员将总立体角${4\pi }$分解为$N \times  N$个相等的极角$\theta$和方位角$\psi$的段。这使得靠近极点$\left( {\theta  = 0,\pi }\right)$的${\Omega }_{i}$非常小，而靠近赤道$\left( {\theta  = \pi /2}\right)$的${\Omega }_{i}$很大。在他们的方法中，称为${\mathrm{{FT}}}_{n}$有限体积法，他们建议，对于$n$个不同的极角${\theta }_{i}$，应该在极点附近选择较少的方位角，即随着${\theta }_{i}$的增长，分布为$4,8,\ldots ,{2n} - 4,{2n},{2n},{2n} - 4,\ldots ,8,4$。这导致$n\left( {n + 2}\right)$个不同的立体角（等于标准${S}_{n}$方案中的坐标数），所有${\Omega }_{i}$的大小大致相等。${\mathrm{{FT}}}_{n}$角度离散化方案也已被其他研究人员探索[148,149]。Kamdem [76]建议，可以通过在端点附近使用更高的方位角细化来提高FAM的准确性。例如，如果在某个二维问题中$\psi$在0和$\pi$之间离散化，则建议在$\psi  = 0$和$\psi  = \pi$附近进行更高的角度细化。

有限角法已被用于非结构化网格，以模拟复杂的二维和三维几何形状[59,144,150]。近年来，该方法也已用于非灰体介质。Sun等人[151]展示了其在由分子气体组成的三维介质中的应用。他们使用了一种混合的蒙特卡罗/FAM方法，其中蒙特卡罗方法仅用于处理介质的非灰体性质（波数选择），而FAM则用于在非结构化网格上求解RTE。Liu及其同事[152]也展示了如何使用域分解将具有非结构化网格的有限角法并行化。该方法还被用于许多组合传热问题[153,154]，并被包含在几个重要的商业CFD代码中，例如FLUENT [155]，其中也使用了非结构化的有限体积法。近年来，有限元法已被证明可用于RTE的空间和角度离散化[156-158]。FEM用于角度离散化的一个主要优点是，类似于其用于空间离散化的优点，它可以优雅地扩展以支持角度自适应[157]。


<!-- Meanless: The Method of Discrete Ordinates ( ${S}_{N}$ -Approximation) Chapter $\mid$-->

<!-- Media -->

<!-- figureText: 0.8 0.8 0.6 $y/L$ 0.4 0.2 0 -0.5 -0.3 -0.1 0.1 0.3 0.5 $x/L$ (b) 0.6 $y/L$ 0.4 0.2 0 -0.5 -0.3 -0.1 0 0.1 0.3 0.5 $x/L$ (a) -->

<img src="https://cdn.noedgeai.com/bo_d141kgv7aajc73bva480_38.jpg?x=354&y=225&w=1131&h=563&r=0"/>

图 16.11 在与例16.9中考虑的几何形状类似的几何形状中，辐射平衡的温度分布：(a) 使用${S}_{8}$近似的DOM，(b) 使用与${S}_{8}$近似完全相同的方向向量的FAM。深色阴影表示热区域，而浅色阴影表示冷区域。结果转载自[162]。

<!-- Media -->

与标准DOM一样，FAM也已经被公式化并证明适用于具有空间变化折射率（渐变介质）的介质[159,160]。特别是，Liu [159]使用FAM进行角度离散化和有限体积法进行三维空间离散化，为渐变介质制定了RTE的离散形式。使用了结构化正交网格和笛卡尔坐标。该方法针对一维平板的已发表结果进行了验证，并为二维矩形进行了演示。最近，Asllanaj和Fumeron [160]将FAM应用于使用非结构化网格离散化的复杂二维渐变介质。

## DOM与FAM的比较

历史上，有限角法的开发初衷是为了守恒辐射能，就像空间有限体积法（作为有限差分法的竞争者）被开发出来以保证通量的局部和全局守恒一样。自其诞生以来，少数研究[139,150,161,162]直接比较了FAM和标准DOM，即在两种计算中使用相同数量的方向和单元格。这些研究表明，FAM的一个附带好处是，在大多数DOM结果中明显存在的射线效应在一定程度上得到了缓解。从概念上讲，对立体角进行积分具有平均效应，因此有望减轻沿坐标方向的能量“条纹化”。为了研究这两种方法中的射线效应，Mazumder及其同事[150,162]在与例16.9中考虑的类似的二维几何形状中进行了计算，即一个方形空腔，除了底部壁上$- {0.05} \leq  x \leq  {0.05}$之间有一个加热斑块外，所有侧面都是冷黑壁。一个类似的测试问题——加热带放置在$- {0.1} \leq  x \leq  {0.1}$之间是唯一的区别——现在被辐射界视为一个基准问题，并已被其他几位研究人员[74,143,163,164]以及例16.9中考虑。虽然概念上简单，但已发现它对大多数RTE求解方法都非常具有挑战性。计算是在一系列光学厚度下，在纯各向同性散射介质以及仅有吸收和发射的介质中进行的。如果辐射平衡存在，各向同性散射的灰色介质与纯吸收-发射（非散射）的灰色介质表现相同，因为局部吸收后接着再发射等同于各向同性散射。这里，为了突出射线效应，图16.11显示了在光学厚度${\kappa L} = {0.01}$的纯吸收-发射介质中进行辐射平衡计算的结果。这些结果是使用精细的（${160} \times  {160}$）网格获得的。在DOM结果中，沿坐标方向的能量条纹化（射线效应）清晰可见。在FAM结果中也存在一些条纹化的迹象，特别是当远离热斑，即更靠近顶壁时。然而，平均而言，射线效应似乎被FAM有所缓解。Sankar和Mazumder [150]解决了例16.9中考虑的相同问题，但使用了由9,492个三角形单元格离散化的精细非结构化网格。DOM和FAM都使用阶梯格式进行空间离散化，而使用${S}_{8}$格式（或FAM中等效数量的角度，即$1 \times  {40}$）进行角度离散化。图16.12显示，如果精细的空间网格与精细的角度网格耦合，FAM会产生准确的结果。尽管预测的热流中仍有轻微的振荡，但它们不像使用粗糙的角度网格和精细的空间网格时那么明显（参见例16.9）。另一方面，DOM仍然表现出强烈的射线效应，表现为两个壁上热流的剧烈振荡。


<!-- Meanless: 502 Radiative Heat Transfer-->

<!-- Media -->

<!-- figureText: 0.08 0.12 无量纲热流, $q/\sigma {T}_{h}^{4}$ 0.1 0.08 0.06 0.04 $\square$ 蒙特卡罗 0.02 $\operatorname{DOM}\left( {S}_{8}\right)$ FAM (1 x 40 angles) ${}^{0}$ 0.2 0.4 0.6 0.8 距离底壁的距离, $y/L$ (b) 无量纲热流, $q/\sigma {T}_{h}^{4}$ 0.06 0.04 0.02 $\square$ 蒙特卡罗 $\operatorname{DOM}\left( {S}_{8}\right)$ FAM (1 x 40 angles) 0 0.2 0.4 0.6 0.8 距离左壁的距离, $x/L$ (a) -->

<img src="https://cdn.noedgeai.com/bo_d141kgv7aajc73bva480_39.jpg?x=320&y=225&w=1193&h=607&r=0"/>

图 16.12 例16.9中考虑的问题的壁面热流，光学厚度为 ${\tau }_{L} = 1$：(a) 顶壁，(b) 侧壁。结果转载自[150]。${T}_{h}$ 是底部热斑的温度。

<!-- Media -->

在其他研究中，Liu及其同事[165]用通用贴体（曲线）坐标[166]表达了RTE，并将标准离散坐标法和有限角法应用于许多二维和三维问题。他们发现两种方法所需的CPU时间相似，而有限角法总是稍微更准确。Fiveland和Jessee [167]以及Kim和Huh [168]也得出了类似的结论，注意到有限角法在光学薄介质中尤其优于标准离散坐标法，因为它对射线效应不那么敏感。Coelho及其同事[169]比较了有限角法与离散传递法[88]的性能，并且像Selçuk和Kayakol [87]一样，发现有限角法经济得多。有限角法的主要优点是选择坐标的自由度更大，以及有限角法守恒辐射能。此外，处理复杂外壳对于有限角法来说更自然。例如，Baek及其同事[170-172]使用贴体坐标来研究具有灰色、恒定属性介质的几个三维外壳中的辐射。

### 16.7 修正离散坐标法

在第16.5节中已经注意到，离散坐标法（以其标准或有限角形式）可能会受到射线效应的影响，如果方向离散化相对于空间离散化较粗，并且如果介质包含强发射的小源（来自壁面或介质内部）。这促使Ramankutty和Crosbie [173,174]将边界发射与介质发射分开，就像在第15.8节的修正微分近似中所做的那样，即令

$$
I\left( {\mathbf{r},\widehat{\mathbf{s}}}\right)  = {I}_{w}\left( {\mathbf{r},\widehat{\mathbf{s}}}\right)  + {I}_{m}\left( {\mathbf{r},\widehat{\mathbf{s}}}\right) . \tag{16.90}
$$

与壁面相关的强度场可以通过第16.5节中概述的任何标准方法求解，而 ${I}_{m}$ 的RTE和边界条件变为

$$
\frac{d{I}_{m}}{ds} = \kappa {I}_{b} - \beta {I}_{m}\left( \widehat{\mathbf{s}}\right)  + \frac{{\sigma }_{s}}{4\pi }{\int }_{4\pi }{I}_{m}\left( {\widehat{\mathbf{s}}}^{\prime }\right) \Phi \left( {{\widehat{\mathbf{s}}}^{\prime },\widehat{\mathbf{s}}}\right) d{\Omega }^{\prime } + \frac{{\sigma }_{s}}{4\pi }{\int }_{4\pi }{I}_{w}\left( {\widehat{\mathbf{s}}}^{\prime }\right) \Phi \left( {{\widehat{\mathbf{s}}}^{\prime },\widehat{\mathbf{s}}}\right) d{\Omega }^{\prime }, \tag{16.91}
$$




<!-- Meanless: The Method of Discrete Ordinates (S ${}_{N}$ -Approximation) Chapter $\mid$-->

$$
{I}_{m}\left( {{\mathbf{r}}_{w},\widehat{\mathbf{s}}}\right)  = \frac{1 - \epsilon }{\pi }{\int }_{\widehat{\mathbf{n}} \cdot  {\widehat{\mathbf{s}}}^{\prime } < 0}{I}_{m}{\widehat{\mathbf{s}}}^{\prime })\left| {\widehat{\mathbf{n}} \cdot  {\widehat{\mathbf{s}}}^{\prime }}\right| d{\Omega }^{\prime }. \tag{16.92}
$$

再次，方程(16.91)和(16.92)对灰体介质有效，或在光谱基础上对非灰体介质有效，适用于具有不透明、漫射发射和漫反射壁面的外壳。与方程(16.4)和(16.5)比较，很明显，${I}_{m}$ 的修正方程可以通过将方程(16.4)中的发射项 $\kappa {I}_{b}$ 替换为

$$
\kappa {I}_{b} + \frac{{\sigma }_{s}}{4\pi }\mathop{\sum }\limits_{{j = 1}}^{n}{w}_{j}{I}_{w}\left( {\widehat{\mathbf{s}}}_{j}\right) \Phi \left( {{\widehat{\mathbf{s}}}_{j},{\widehat{\mathbf{s}}}_{i}}\right) .
$$

来以相同的方式求解。对于给定的温度和辐射属性场，与壁面相关的强度只需计算一次，在存在散射、壁面反射或高阶空间格式时，${I}_{m}$ 所需的迭代过程中无需重新评估。或者，方程(16.91)和(16.92)也可以用有限角法求解，方法是将与 ${I}_{w}$ 相关的散射源添加到方程(16.79a)的 ${S}_{i,P}$ 中。Ramankutty和Crosbie将修正的离散坐标法应用于具有黑体壁的边界发射主导的矩形外壳，使用了 ${I}_{w}$ 的解析解。对于这类问题，当使用低阶 ${S}_{N}$ -技术和精细空间网格时，标准离散坐标法表现出强烈的射线效应，而修正的离散坐标法完全减轻了这种效应。Sakami及其同事[175,176]和Baek及其同事[172,177]研究了具有这种修正离散坐标法的不规则几何形状，后者对壁面发射部分采用了直接的蒙特卡罗方案，即使在极端的边界源项下也给出了极好的精度。Sankar和Mazumder [150]更进一步，对 ${I}_{m}$ 调用了 ${P}_{1}$ 近似，而 ${I}_{w}$ 则使用有限角法确定。这在第15.8节中被描述为混合微分近似，并且与修正微分近似基本相同，只是壁面发射的强度分量是使用FAM而不是精确的角系数公式求解的。这使得在复杂几何形状中，特别是有障碍物的情况下，求解相对容易，在这些情况下，计算角系数矩阵本身就是一项艰巨的任务。Sankar和Mazumder [150]报告的结果显示出与前述关于修正DOM的工作相似的趋势：该方法被发现在壁面发射主导的情况下非常准确。

然而，与修正微分近似一样，修正离散坐标法或其变体在介质内部存在强烈的、孤立的发射源时，无法提高基础标准方法的准确性。这一点也得到了Sankar和Mazumder [150]的证实。Coelho [73]和Li和Werther [178]扩展了该方法，将所有发射（壁面和介质）与 ${I}_{w}$ 关联起来，只有散射辐射通过离散坐标计算。虽然这减轻了所有射线效应，但尚不清楚这种方法有什么增益（因为与发射相关的强度的求解基本上与整个问题一样复杂）。

### 16.8 偶宇称公式

方向性RTEs，即标准离散坐标法求解的方程(16.4)或有限角法求解的方程(16.79)，是具有双曲特性的一阶偏微分方程，可能导致不物理的结果（负强度）并可能显示出强烈的射线效应。将它们转换为二阶边值问题可以消除不物理的结果，并可能减轻射线效应。在传热领域，这首先由Song和Park [179]尝试。通过定义偶宇称和奇宇称强度，可以获得二阶公式，这涉及到沿一条线的前向$\widehat{\mathbf{s}}$和后向-$\widehat{\mathbf{s}}$两个方向：

$$
\psi \left( {\mathbf{r},\widehat{\mathbf{s}}}\right)  = I\left( {\mathbf{r},\widehat{\mathbf{s}}}\right)  + I\left( {\mathbf{r}, - \widehat{\mathbf{s}}}\right) , \tag{16.93}
$$

$$
\phi \left( {\mathbf{r},\widehat{\mathbf{s}}}\right)  = I\left( {\mathbf{r},\widehat{\mathbf{s}}}\right)  - I\left( {\mathbf{r}, - \widehat{\mathbf{s}}}\right) . \tag{16.94}
$$

在接下来的讨论中，我们将自己限制在各向同性散射的更简单情况。在这些条件下，将方向 $\widehat{\mathbf{s}}$ 和 $- \widehat{\mathbf{s}}$ 的RTE (16.1) 相加和相减，得到

$$
\widehat{\mathbf{s}} \cdot  \nabla \phi \left( {\mathbf{r},\widehat{\mathbf{s}}}\right)  =  - \beta \left( \mathbf{r}\right) \psi \left( {\mathbf{r},\widehat{\mathbf{s}}}\right)  + {2\kappa }\left( \mathbf{r}\right) {I}_{b}\left( \mathbf{r}\right)  + \frac{2{\sigma }_{s}\left( \mathbf{r}\right) }{4\pi }{\int }_{2\pi }\psi \left( {\mathbf{r},{\widehat{\mathbf{s}}}^{\prime }}\right) d{\Omega }^{\prime }, \tag{16.95}
$$

$$
\widehat{\mathbf{s}} \cdot  \nabla \psi \left( {\mathbf{r},\widehat{\mathbf{s}}}\right)  =  - \beta \left( \mathbf{r}\right) \phi \left( {\mathbf{r},\widehat{\mathbf{s}}}\right) . \tag{16.96}
$$




<!-- Meanless: 604 Radiative Heat Transfer-->

注意，现在对所有方向的积分只覆盖 ${2\pi }$（半个单位球）。使用方程(16.96)，奇宇称项很容易被消除，从而得到一个抛物型的二阶偏微分方程，

$$
\widehat{\mathbf{s}} \cdot  \nabla \left( {\frac{1}{\beta \left( \mathbf{r}\right) }\widehat{\mathbf{s}} \cdot  \nabla \psi \left( {\mathbf{r},\widehat{\mathbf{s}}}\right) }\right)  = \beta \left( \mathbf{r}\right) \psi \left( {\mathbf{r},\widehat{\mathbf{s}}}\right)  - {2\kappa }\left( \mathbf{r}\right) {I}_{b}\left( \mathbf{r}\right)  - \frac{2{\sigma }_{s}\left( \mathbf{r}\right) }{4\pi }{\int }_{2\pi }\psi \left( {\mathbf{r},{\widehat{\mathbf{s}}}^{\prime }}\right) d{\Omega }^{\prime }, \tag{16.97}
$$

这需要在双向射线的每一端（边界交点）都有一个边界条件。将方程(16.93)和(16.94)代入(16.2)得到

$$
\psi \left( {{\mathbf{r}}_{w},\widehat{\mathbf{s}}}\right)  + \phi \left( {{\mathbf{r}}_{w},\widehat{\mathbf{s}}}\right)  = {2\epsilon }\left( {\mathbf{r}}_{w}\right) {I}_{b}\left( {\mathbf{r}}_{w}\right)  + \frac{\rho \left( {\mathbf{r}}_{w}\right) }{\pi }{\int }_{\widehat{\mathbf{n}} \cdot  {\widehat{\mathbf{s}}}^{\prime } < 0}\left\lbrack  {\psi \left( {{\mathbf{r}}_{w},{\widehat{\mathbf{s}}}^{\prime }}\right)  - \phi \left( {{\mathbf{r}}_{w},{\widehat{\mathbf{s}}}^{\prime }}\right) }\right\rbrack  \left| {\widehat{\mathbf{n}} \cdot  {\widehat{\mathbf{s}}}^{\prime }}\right| d{\Omega }^{\prime },\;\widehat{\mathbf{n}} \cdot  \widehat{\mathbf{s}} > 0. \tag{16.98}
$$

由于 $\psi$ 和 $\phi$ 的双向性，边界条件也必须应用于 $- \widehat{\mathbf{s}}$ 方向，或

$$
\psi \left( {{\mathbf{r}}_{w},\widehat{\mathbf{s}}}\right)  - \phi \left( {{\mathbf{r}}_{w},\widehat{\mathbf{s}}}\right)  = {2\epsilon }\left( {\mathbf{r}}_{w}\right) {I}_{b}\left( {\mathbf{r}}_{w}\right)  + \frac{\rho \left( {\mathbf{r}}_{w}\right) }{\pi }{\int }_{\widehat{\mathbf{n}} \cdot  {\widehat{\mathbf{s}}}^{\prime } < 0}\left\lbrack  {\psi \left( {{\mathbf{r}}_{w},{\widehat{\mathbf{s}}}^{\prime }}\right)  + \phi \left( {{\mathbf{r}}_{w},{\widehat{\mathbf{s}}}^{\prime }}\right) }\right\rbrack  \left| {\widehat{\mathbf{n}} \cdot  {\widehat{\mathbf{s}}}^{\prime }}\right| d{\Omega }^{\prime },\;\widehat{\mathbf{n}} \cdot  \widehat{\mathbf{s}} < 0. \tag{16.99}
$$

这两个只在 $\widehat{\mathbf{n}} \cdot  \widehat{\mathbf{s}}$ 的符号上有所不同；它们可以很容易地组合起来，并且 $\phi$ 被消除，从而给出了由方程(16.97)描述的射线两端的适当边界条件

$$
\psi \left( {{\mathbf{r}}_{w},\widehat{\mathbf{s}}}\right)  - \operatorname{sign}\left( {\widehat{\mathbf{n}} \cdot  \widehat{\mathbf{s}}}\right) \frac{1}{\beta \left( \mathbf{r}\right) }\widehat{\mathbf{s}} \cdot  \nabla \psi \left( {\mathbf{r},\widehat{\mathbf{s}}}\right) 
$$

$$
 = {2\epsilon }\left( {\mathbf{r}}_{w}\right) {I}_{b}\left( {\mathbf{r}}_{w}\right)  + \frac{\rho \left( {\mathbf{r}}_{w}\right) }{\pi }{\int }_{\widehat{\mathbf{n}} \cdot  {\widehat{\mathbf{s}}}^{\prime } < 0}\left\lbrack  {\psi \left( {{\mathbf{r}}_{w},{\widehat{\mathbf{s}}}^{\prime }}\right)  + \operatorname{sign}\left( {\widehat{\mathbf{n}} \cdot  \widehat{\mathbf{s}}}\right) \frac{1}{\beta \left( \mathbf{r}\right) }\widehat{\mathbf{s}} \cdot  \nabla \psi \left( {\mathbf{r},{\widehat{\mathbf{s}}}^{\prime }}\right) }\right\rbrack  \left| {\widehat{\mathbf{n}} \cdot  {\widehat{\mathbf{s}}}^{\prime }}\right| d{\Omega }^{\prime }. \tag{16.100}
$$

在偶宇称公式中，方向离散化仅在 ${2\pi }$ 或半球上进行，即只需要一半的坐标，但需要求解二阶偏微分方程。空间离散化可以通过中心差分[179]、有限体积[180,181]或有限元[112, 167,182,183]来完成；方向离散化通常是根据标准离散坐标法进行的，除了Becker及其同事[183]也采用有限元进行方向离散化。Liu等人[181]、Fiveland和Jessee [167]以及Kang和Song [182]将标准离散坐标法的实现与偶宇称公式进行了比较，发现总的来说，精度相似，但偶宇称模型需要更多的CPU时间，特别是在光学薄的情况下。

### 16.9 其他相关方法

辐射的方向性所带来的困难促使许多研究人员开发各种近似方案来离散化方向（类似于标准离散坐标法）或对立体角范围进行平均（类似于有限角法）。在这里，我们将简要讨论这些相关模型中最重要的几个。

通量法。我们在例16.7中已经看到，双通量法或Schuster-Schwarzschild近似只不过是有限角法的一种粗略实现。双通量法已被许多研究者应用于一维问题，通常是在辐射通量的精确确定不那么关键的情况下[184-189]。注意到该方法的准确性或多或少仅限于没有各向异性散射的一维问题，Chin和Churchill [190,191]开发了六通量法，主要用于处理存在准直辐照时的强各向异性散射。在这种方法中，强度被分解为前向、后向和四个（通常相等）侧向分量。这些分量如何分解或在总立体角 ${4\pi }$ 上平均是任意的。该方法已被各种研究人员应用于多个二维和三维问题[191-194]。

离散传递法（DTM）。离散传递法（DTM）由Lockwood和Shah [88]开发，比离散坐标法（DOM）的认真发展早了几年。离散传递法与离散坐标法相似，因为它也选择离散的方向。它也与蒙特卡罗法（见第20章）有关，因为强度的射线是从一个表面追踪到另一个表面的。本质上，在壳体的边界上建立节点，从这些节点向预定方向发出射线。然后追踪射线穿过内部有限体积时的相互作用（通过在体积中沉积能量而减弱射线，通过发射和入散射而增强），直到它碰到另一个表面。该方法与标准DOM非常相似，并带有其缺点（不守恒，易受射线效应影响，对于非黑体壁和/或散射需要迭代），但使用了更繁琐的光束追踪（等同于指数格式）。在发射点建立一个立体角“铅笔”，当光束穿过内部体积时会产生不准确性。射线追踪与蒙特卡罗法非常相似。虽然DTM避免了蒙特卡罗法的统计散射，但效率较低，因为每条射线携带单一的统计信息（给定位置、方向、波长），而蒙特卡罗模拟中统计选择的束携带多种统计信息（在非灰色、反射和/或散射环境中尤为重要）。尽管如此，由于其早期出现、早期处理不规则几何形状的能力，以及已被纳入几个重要的商业CFD代码中，该方法仍然享有一定的普及度。Cumber [195,196]为该方法提供了几项改进，Coelho和Carvalho [197]提出了DTM的守恒公式，Versteeg及其同事[198,199]量化了该方法固有的误差。多位研究人员进行了二维和三维计算，主要用于炉膛和其他燃烧应用[200-207]。还对DTM与其他方法进行了一些比较。Selçuk和Kayakol [87]比较了二维矩形几何的DTM和 ${S}_{4}$ 解，发现 ${S}_{4}$ 方法在精度相当的情况下，计算机时间少了三个数量级。同样，Coelho及其同事[169]解决了具有障碍物的二维外壳中的辐射问题，使用了DTM、DOM、FAM、蒙特卡罗和区域法（见第17章）；他们发现DOM和FAM最经济。最后，Keramida及其同事使用DTM和六通量法模拟了天然气燃烧炉[194]，发现六通量法更优越。在这些传统应用领域之外，DTM最近也被用于求解多维渐变介质的RTE [208]。


<!-- Meanless: The Method of Discrete Ordinates (S ${}_{N}$ -Approximation) Chapter $\mid$-->

YIX法。该方法最初由Tan和Howell [209]开发，辐射传输方程以积分形式表示[例如，参见方程(9.28)]，即入射辐射 $G$ 和辐射通量 $q$ 可以在任何点（介质内部或边界上）作为三重积分（距离该点的距离和立体角）进行评估。从每个点（沿离散坐标）的积分涉及某些可以预先确定的几何函数。虽然该方法有潜力比DOM和FAM更有效，但每个问题的设置工作量很大，并且难以推广。该方法似乎仅由开发它的团队使用[210-214]。

DRESOR法。在许多应用中，需要高角度分辨率地确定强度的角度分布。例如，在许多测量技术中，强度仅沿某些视线测量，这些信息用于确定辐射特性。虽然DOM和FAM都计算沿不同方向的强度，但选择的方向通常相当有限。蒙特卡罗法（见第20章）能够计算任意方向的强度。然而，这需要巨大的计算时间和资源，并且结果常常充满统计误差。DRESOR（散射或反射能量比分布）法由Zhou, Cheng及其同事[215-217]开发，旨在以高角度分辨率计算强度。在DRESOR法中，RTE的源项，即方程(16.1)的右侧，首先被分成两部分：一个表示无自吸收的等温单位体积的发射项，另一个表示介质在远点周围单位体积中散射并进入方向$\widehat{\mathbf{s}}$周围单位立体角的能量项。同样，RTE的边界条件，即方程(16.2)的右侧，被分成两项：一个表示无自吸收的等温单位面积的发射项，另一个表示辐照到边界上一点然后反射到方向$\widehat{\mathbf{s}}$周围单位立体角的能量项。在DRESOR法中，上述散射和反射项用四个比率表示。这些比率使用类似于蒙特卡罗法的射线追踪程序计算。正是这种计算使得DRESOR法在计算上非常昂贵。近年来，为了减少该方法的计算开销，所谓的方程求解DRESOR法[218,219]已被提出并用于三维问题。自诞生以来，DRESOR法已扩展到一维瞬态RTE的解[220,221]、具有渐变折射率的一维（平面平行和圆柱形）介质[221,222]以及平面平行非灰介质[223]。

除了上面讨论的具体方法外，还提出了对DOM和FAM的其他变体和改进。最近，一种称为广义源有限体积法（GS-FVM）[224]的方法被开发出来。该方法从根本上基于FAM，并且也能够以高角度分辨率计算强度。发现该方法生成的强度与使用反向蒙特卡罗法生成的结果吻合得很好。为了减轻射线效应，提出了一种基于所谓入射能量传输方程的解的方法[225]。在这种方法中，不是求解标准的RTE（用于强度），而是为入射辐射开发了一个新方程。该方法仅对平面平行介质进行了演示。


<!-- Meanless: 606 Radiative Heat Transfer-->

### 16.10 结论

有限角法，由于其普遍适用性以及能够通过精细的角度网格生成高精度解的能力，如今可能是解决复杂情况下RTE的最流行方法。当然，该方法也有其弱点：它会受到射线效应的影响（尽管比DOM轻），并且当精细的空间网格必须伴随精细的角度网格才能获得精确解时，计算成本可能非常高。此外，DOM和FAM都受到虚假散射的影响。两种方法在光学厚介质和/或强散射介质中收敛性差。尽管有这些缺点，DOM和FAM是唯一（除了统计蒙特卡罗法）可以通过细化空间和角度网格，在不增加其数学复杂性的情况下，达到任意精度水平的方法。

离散坐标法与第17章的区域法有关，后者是直到20世纪60年代首选的RTE求解方法。区域法在光学厚介质中也表现不佳，但它不受射线效应和虚假散射的影响。然而，它不能在存在各向异性散射的情况下使用，并且难以应用于不规则几何。另一方面，${P}_{1}$ -近似更容易应用和求解，并且对于光学厚介质（DOM和FAM难以应用的地方）以及在冷环境中的热辐射介质，给出了非常准确的结果。但该方法在由辐射壁包围和/或受方向性辐照的光学薄介质中会失效。最后，也许最重要的是，更高阶的 ${P}_{N}$ -近似在数学上难以公式化。

## 问题

16.1 考虑一个被大的、等温的、灰色的板所包含的灰色、等温且各向同性散射的介质，板的温度分别为 ${T}_{1}$ 和 ${T}_{2}$，发射率分别为 ${\epsilon }_{1}$ 和 ${\epsilon }_{2}$。使用 ${S}_{2}$ -近似确定板之间的辐射通量。

16.2 考虑一个大的、等温（温度 ${T}_{w}$）、灰色且漫射（发射率 $\epsilon$）的壁面，毗邻一个半无限的灰色吸收/发射和线性各向异性散射介质。该介质是等温的（温度 ${T}_{m}$）。使用 ${S}_{2}$ -近似确定作为离板距离函数的辐射通量。

16.3 考虑相距 $1\mathrm{\;m}$ 的平行黑板，温度恒定为 ${T}_{1}$ 和 ${T}_{2}$。由于压力变化，（灰色）吸收系数等于

$$
\kappa  = {\kappa }_{1} + {\kappa }_{1}x;\;{\kappa }_{0} = {0.01}{\mathrm{\;{cm}}}^{-1};\;{\kappa }_{1} = {0.0002}{\mathrm{\;{cm}}}^{-2},
$$

其中 $x$ 从板1测量。介质不散射辐射。通过精确方法和 ${S}_{2}$ -近似确定辐射平衡下的无量纲热流 $\Psi  = q/\sigma \left( {{T}_{1}^{4} - {T}_{2}^{4}}\right)$。

16.4 半径为 ${100\mu }\mathrm{m}$ 的黑色球形颗粒悬浮在相距 $1\mathrm{\;m}$ 的两块冷黑平行板之间。颗粒以 $\pi /{10}\mathrm{\;W}/$ 颗粒的速率产生热量，这些热量必须通过热辐射去除。板间的颗粒数量由下式给出

$$
{N}_{T}\left( z\right)  = {N}_{0} + {\Delta Nz}/L,\;0 < z < L;\;{N}_{0} = {\Delta N} = {212}\text{ 颗粒 }/{\mathrm{{cm}}}^{3}.
$$

(a) 确定局部吸收系数和局部产热率；引入光学坐标并确定整个间隙的光学厚度。

(b) 如果要采用 ${S}_{2}$ -近似，控制传热的相关方程和边界条件是什么？

(c) 顶部和底部表面的热流率是多少？颗粒释放的总能量是多少？最大颗粒温度是多少？


<!-- Meanless: The Method of Discrete Ordinates ( ${S}_{N}$ -Approximation) Chapter $\mid$-->

16.5 两个无限长的同心圆柱体，半径为 ${R}_{1}$ 和 ${R}_{2}$，发射率为 ${\epsilon }_{1}$ 和 ${\epsilon }_{2}$，具有相同的恒定表面温度 ${T}_{w}$。圆柱体之间的介质具有恒定的吸收系数 $\kappa$ 并且不散射；介质内部发生均匀的热生成 ${Q}^{\prime \prime \prime }$。如果辐射是唯一的传热方式，使用 ${S}_{2}$ -近似确定介质中的温度分布和壁面上的热流。

16.6 一个无限的、黑色的、等温的板界定了一个充满黑球的半无限空间。在离板任意给定距离 $z$ 处，颗粒数密度相同，即 ${N}_{T} = {6.3662} \times  {10}^{8}{\mathrm{\;m}}^{-3}$。然而，悬浮球体的半径随着离表面的距离单调减小

$$
a = {a}_{o}{e}^{-z/L};\;{a}_{o} = {10}^{-4}\mathrm{\;m},\;L = 1\mathrm{\;m}.
$$

(a) 确定作为 $z$ 的函数的吸收系数（你可以做大颗粒假设）。

(b) 确定作为 $z$ 的函数的光学坐标。半无限空间的总光学厚度是多少？

(c) 假设辐射平衡存在，并使用 ${S}_{2}$ -近似，建立边界条件并求解热流和温度分布（作为 $z$ 的函数）。

16.7 考虑两块平行的黑板，都在 ${1000}\mathrm{\;K}$，相距 $2\mathrm{\;m}$。板间的介质发射和吸收（但不散射），吸收系数为