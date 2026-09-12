# A 题药材烘干模型参考文献使用指南

> 本文档用于指导论文与各分题方法文档的文献引用。它回答三个问题：哪些论断需要外部文献、应把引用放在哪里、哪些结论不能借文献过度背书。题面参数、附件数据和本项目计算结果应分别引用官方材料或项目输出，不能用外部论文代替。

## 1. 引用原则

### 1.1 三类来源必须分开

1. **官方材料**：题目给出的几何尺寸、初始条件、物性公式、边界数据、收缩半径和结果模板，应引用 `data/official/A题/` 下的题面与附件。
2. **外部学术文献**：用于说明模型形式、常见假设、适用范围、数值方法和干燥机理具有学术依据。
3. **本项目结果**：57.5 h、温湿度剖面、网格收敛、守恒残差和敏感性结果，应引用程序输出、图表或项目内验证文档，不能引用外部论文。

### 1.2 外部文献不能替题面参数背书

本项目中的 \(h=25\ \mathrm{W/(m^2K)}\)、\(h_m=8\times10^{-7}\ \mathrm{m/s}\)、\(D(C,T)\) 经验式以及 \(R(t)\) 数据均来自题面或附件。外部论文只能说明类似参数和模型在食品或药材干燥中常见，不能把其他材料的测量值写成本文参数来源。

### 1.3 优先引用顺序

同一论断存在多篇候选文献时，优先级如下：

1. 与药材、芳香植物或长圆柱几何直接相关的研究；
2. `Journal of Food Engineering`、`Drying Technology` 等领域期刊中的方法或综述；
3. 数值分析期刊中的有限体积、退化扩散和时间推进文献；
4. 一般性综述只用于研究背景，不用来支撑非常具体的离散公式。

正文通常每个关键论断引用 1–2 篇即可。不要在同一句后堆叠五六篇内容重复的文献。

## 2. 快速引用地图

| 模块 | 需要文献支持的论断 | 推荐文献键 | 建议插入位置 |
| --- | --- | --- | --- |
| 研究背景 | 热风干燥是耦合传热传质过程；建模可用于预测内部温湿分布 | R01、R02、R03 | 论文引言；`PROJECT.md` §2 的论文化版本 |
| 药材背景 | 药用或芳香植物的干燥速率和品质受温度、有效扩散率影响 | R04、R05 | 论文引言末段 |
| Q1 热方程 | 傅里叶导热与能量守恒构成瞬态热模型 | R01、R03 | `docs/q1/heat_conduction.md` §5 |
| Q1 水分方程 | 食品/植物材料常用有效 Fick 扩散描述内部水分迁移 | R02、R03、R06 | `docs/q1/moisture_diffusion.md` §6–7 |
| 长圆柱降维 | 细长枝条或圆柱体可在条件合适时采用一维径向扩散模型 | R06、R07 | `docs/q1/moisture_diffusion.md` §3、§15；Q3 方法 1 |
| Robin 边界 | 表面值应由内部扩散与外部对流传质共同决定 | R08、R09 | `docs/q1/moisture_diffusion.md` §10；Q2 §11 |
| Biot 数 | Biot 数用于比较内外传递阻力，且可随干燥状态变化 | R09 | `docs/q1/moisture_diffusion.md` §15 |
| 界面平衡 | 空气湿度与固体干基含水率通常需通过平衡含水率/吸附等温线关联 | R10 | `docs/q1/moisture_diffusion.md` §12；Q2、Q4 模型局限 |
| Q2 热湿耦合 | 温度、含水率和变物性使干燥问题呈非线性耦合 | R01、R02、R03 | `docs/q2/coupled_heat_moisture_derivation.md` §4、§12 |
| 潜热局限 | 较完整模型可包含蒸发冷却或吸附热；本文省略潜热属于题设闭合简化 | R02、R11 | Q1 §1；Q2 §11.2、§19 |
| 有限体积法 | 有限体积法适合守恒型扩散方程，并可处理 Robin 边界 | R07、R12、R13 | Q1 数值方法；Q2 §15；Q4 §16 |
| Rannacher 启动 | 初始数据或边界不光滑时，以后向欧拉启动可改善 Crank–Nicolson 的早期行为 | R14 | `docs/q1/heat_conduction.md` §9 |
| Q3 退化扩散 | 浓度相关扩散可用 Kirchhoff 变换处理；退化扩散需使用稳健离散 | R15、R16 | Q3 独立复核 §2.2；Q4 §16.2–16.3 |
| 真实表面硬化 | 食品干燥中确有表面硬化和结构变化，但需力学/多相机理支持 | R17 | Q3 模型评价，不能用于证明旧离散停滞 |
| Q4 收缩 | 收缩会改变几何尺度、孔隙结构和热质传递过程 | R18、R19、R20 | `docs/q4/moving_boundary_heat_moisture_derivation.md` §1、§5、§15.2 |
| 移动边界 | 收缩干燥可采用移动边界或随材料运动的坐标方法 | R19、R21 | Q4 §8–10、§17 |

## 3. 分题引用指南

### 3.1 问题 1：常物性热传导与水分扩散

#### 3.1.1 研究对象与一维径向简化

建议在 `docs/q1/moisture_diffusion.md` §3 或论文“模型假设”中，在说明忽略轴向梯度之后加入 R06。该文研究圆柱形巴拉圭冬青枝条，明确使用长圆柱的一维 Fick 扩散模型，和本项目的几何最接近。

R07 可作为数值建模补充：其研究圆柱状食品中的一维瞬态扩散，采用全隐式有限体积方法，并比较固定/变化体积和扩散系数模型。

建议表述：

> 对于长径比较大且主要传递阻力沿径向分布的枝条状植物材料，一维长圆柱扩散模型已有干燥研究采用。因此，本文在轴对称条件下忽略轴向梯度，仅保留径向传热传质 [R06–R07]。

注意：文献只能说明该简化具有先例。本项目仍应使用 \(L/R=12.5\)、端面面积比例或二维敏感性分析说明它对当前几何是否足够准确。

#### 3.1.2 傅里叶导热与 Fick 扩散

在 `docs/q1/heat_conduction.md` §5 的控制方程之后引用 R01 或 R03；在 `docs/q1/moisture_diffusion.md` §6–7 引用 R02、R03。R02 特别适合说明：简单扩散模型有明确适用边界，扩散系数随含水率变化时不能盲目套用常系数解析式。

建议表述：

> 食品与植物材料热风干燥通常由内部热传导、有效水分扩散以及表面对流交换共同控制；在连续介质近似下，可分别由傅里叶定律和有效 Fick 扩散建立守恒方程 [R01–R03]。

#### 3.1.3 第三类传质边界和 Biot 数

在 `docs/q1/moisture_diffusion.md` §10 的 Robin 条件后引用 R08。该研究直接比较了恒定表面浓度和对流传质两类边界，适合支撑“表面含水率不是预先给定，而由内外阻力共同决定”。

在 §15 的质量传递 Biot 数后引用 R09。避免把 \(\mathrm{Bi}_m>0.1\) 写成绝对定律；更稳妥的写法是“提示内部梯度不可忽略”，因为 Biot 数及控制机制可能随干燥状态变化。

#### 3.1.4 空气含湿量与药材含水率的基准问题

在 `docs/q1/moisture_diffusion.md` §12 必须引用 R10。该文说明药用和芳香植物的空气平衡湿度与固体平衡含水率通常通过 EMC/ERH 关系连接。

建议表述：

> 严格的气固界面传质通常需要吸附等温线或平衡含水率关系，将空气相湿度映射为固体表面的平衡含水状态 [R10]。由于题目未提供该关系，本文把附件中的烘房水分浓度视为可进入题给 Robin 关系的等效外部状态量；该处理是题设下的闭合假设，而非两种 kg/kg 定义之间的恒等关系。

### 3.2 问题 2：变物性热湿耦合

#### 3.2.1 耦合框架

在 `docs/q2/coupled_heat_moisture_derivation.md` §4 或 §12 引用 R01–R03。R02 对本项目尤其重要，因为它指出：干燥中的热质耦合和含水率相关扩散系数往往不可忽略。

建议表述：

> 干燥模型中的温度与水分迁移通常相互影响，且有效扩散率可依赖温度和含水率；因此变物性条件下应采用非线性耦合求解，而不能继续把两个场完全独立处理 [R02–R03]。

这里应紧接着说明本项目的特殊耦合路径是题面公式决定的：\(C\to\rho,c_p,k\to T\)，以及 \((C,T)\to D\to C\)。

#### 3.2.2 不考虑潜热的正确写法

R02 和 R11 不支持“潜热不重要”。相反，它们说明更完整的食品干燥模型可能包含蒸发冷却、相变或吸附热。因此，引用它们时应把省略潜热写成模型层级选择。

建议表述：

> 完整多相干燥模型可进一步考虑蒸发冷却、吸附热及多相水分输运 [R02, R11]。本题未提供把干基含水率通量转换为真实蒸发质量通量所需的界面平衡与体积质量尺度，因此本文采用题面参数可闭合的表观热容模型，不在能量方程中加入不可识别的潜热项。

不要写：

> 药材干燥过程中潜热可以忽略。

#### 3.2.3 有限体积与非线性迭代

在 Q2 §15 引用 R12–R13，说明有限体积离散的核心优势是以界面通量差保持局部守恒，并能够一致处理 Neumann/Robin 边界。

Picard 交替迭代属于常见的工程非线性求解策略，本项目更重要的证据是迭代收敛记录和与独立求解器的交叉验证，不必为了每个算法步骤额外堆积引用。

### 3.3 问题 3：长期干燥时间与退化扩散

#### 3.3.1 终止时间结果不引用外部论文

57.5 h 是由题面参数和本项目程序计算得到的结果，应引用：

- `src/q3_drying_time.py` 的正式输出；
- `outputs/q3/data/` 的未舍入数据；
- `docs/q3/drying_time_independent_review.md` 的独立离散复核。

外部论文只能支撑干燥后期可能由内部扩散控制、扩散率随含水率降低而下降等一般机理，不能证明本项目的具体终止时间。

#### 3.3.2 Kirchhoff 变换与界面平均

在 `docs/q3/drying_time_independent_review.md` §2.2 及 Q4 §16.2–16.3 引用 R15–R16。

- R15 直接讨论浓度相关扩散系数下 Fick 第二定律的 Kirchhoff 变换；
- R16 讨论退化非线性抛物问题的有限体积离散，用于说明强退化扩散需要专门的数值处理。

建议表述：

> 对 \(D=D(C)\) 的非线性扩散，可引入 Kirchhoff 势 \(\Phi(C)=\int^C D(s)\,\mathrm ds\)，使连续扩散通量写成势梯度。相应的界面割线平均能够保留跨单元的积分传输能力，避免点值平均被局部极小扩散系数非物理支配 [R15–R16]。

R15 可以支撑 Kirchhoff 变换思想，但“当前割线平均使本算例在 \(N=200\) 即网格无关”仍是本项目自己的数值结果，应由压力测试和网格收敛表证明。

#### 3.3.3 表面硬化文献的使用边界

R17 是高相关、高被引的真实表面硬化机理文献，可放在模型评价或讨论部分，用于说明食品材料在干燥中可能因结构和力学变化产生表面硬化。

但不得把 R17 放在 `docs/q3/drying_time_method.md` 旧方法 22 后，声称它证明“当前模型预测空气越干越无法完成干燥”。独立复核已经证明该停滞来自顶点型网格与点值调和平均的组合。

正确表述：

> 真实食品干燥中的表面硬化涉及组织收缩、玻璃化、孔隙演变和力学响应等机制 [R17]；当前单一有效扩散模型未显式包含这些机制，因此不能把离散格式导致的通量停滞解释为真实表面硬化。

#### 3.3.4 Rannacher、Richardson 与守恒残差

- 在 Q1 的 Crank–Nicolson 启动说明中引用 R14；其应用领域虽不是食品干燥，但对非光滑初始数据下 Rannacher 启动的数值作用给出了分析依据。
- Richardson 外推是通用数值分析方法。论文若篇幅有限，可说明公式与网格比，不必单列专门文献；若学校模板要求所有数值方法均有来源，可补充标准数值分析教材。
- 有限体积法的机器精度守恒残差只验证通量组装与线性求解闭合，不能单独作为空间精度证据。此处引用 R12 说明守恒结构即可，精度必须依靠网格收敛和独立离散验证。

### 3.4 问题 4：收缩与移动边界

#### 3.4.1 收缩为何必须进入模型

在 `docs/q4/moving_boundary_heat_moisture_derivation.md` §1、§5 或 §15.2 引用 R18、R20。

- R18 是食品对流干燥收缩建模的经典综述；
- R20 直接研究收缩和孔隙率对传热传质的影响。

建议表述：

> 生物材料在脱水过程中常发生显著收缩，几何尺度和孔隙结构的变化会进一步影响内部传热传质，因此固定尺寸模型可能系统性偏离真实干燥过程 [R18, R20]。

#### 3.4.2 移动边界和材料坐标

在 Q4 §8–10 引用 R19；若论文需要说明 ALE 是可选的更一般方法，可在数值方法比较中补充 R21。

建议表述：

> 收缩干燥可表述为随时间变化区域上的移动边界问题，并通过随材料运动的坐标或 ALE 方法映射到计算域 [R19, R21]。本文依据附件仅给出整体半径的条件，进一步采用均匀径向相似收缩，将物理区域映射到固定材料坐标 \(\xi=r/R(t)\)。

需要注意：\(u=r\dot R/R\) 与坐标移动项严格抵消，是本项目假设下的推导结果。文献只能支持移动边界建模路线，不能代替本文自己的坐标变换推导。

#### 3.4.3 识别收缩效应的对照设计

“附录 4 固定半径”与“附录 4 移动半径”的对照属于本项目实验设计，不需要外部引用。R18–R20 可用于解释为什么必须控制物性差异后再归因于收缩。

## 4. 四问代码中的数值方法与引用要求

本节以当前代码为准，区分“代码现在实际做了什么”和“论文最终应如何表述”。引用只能说明方法依据，不能掩盖代码与文档之间的不一致。

### 4.1 总览

| 问题 | 当前代码的主要数值方法 | 必须或建议引用 | 不需要学术引用、靠代码验证即可 |
| --- | --- | --- | --- |
| Q1 | 顶点型径向有限体积；Crank–Nicolson；Rannacher 启动；水分 Picard 迭代；界面调和平均；Thomas 三对角求解；Bessel 解析基准 | R07、R12–R14、R22–R26 | Excel 导出、线性插值、参数校验、Thomas 的具体循环实现 |
| Q2 | 顶点型有限体积；后向欧拉；时间中点冻结物性；分块 Picard 交替；带状矩阵直接求解 | R02、R12–R13、R24、R26 | `solve_banded` 调用、输出定点匹配、工作簿写入 |
| Q3 | 复用 Q2 离散；人工分段变步长；全域阈值检测；网格加密与 Richardson 外推 | R12、R24–R27 | 60 s 输出对齐、平台均值的具体取法、断言语句 |
| Q4 | 材料坐标固定域；顶点型有限体积；Kirchhoff 割线平均；后向欧拉/Picard；移动域守恒；空间和时间收敛 | R15–R16、R19、R21、R24–R25、R28 | 半径文件插值、体外位置留空、四位小数导出规则 |

### 4.2 Q1 代码：`src/q1_fvm.py`

#### 4.2.1 径向有限体积与中心奇点处理

实际代码位置：

- `make_grid`（约第 192 行）：节点位于 \(r_i=i\Delta r\)，中心和表面使用半控制体；
- `assemble_operator`（约第 208 行）：由界面通量散射形成对称三对角算子；
- `check_operator`（约第 222 行）：检查内部行和为零以及 Robin 表面行和。

这属于**顶点型或节点中心有限体积**，不是格心型有限体积。论文中可引用 R12–R13 说明守恒有限体积及一般边界条件的理论基础，并引用 R07 作为圆柱干燥中的相近应用。

建议表述：

> 在径向节点上构造控制体，中心与表面节点分别对应半控制体。对圆柱守恒方程逐控制体积分，将散度项转化为界面通量差，从而自然消除 \(r=0\) 的坐标奇点并保持局部守恒 [R07, R12–R13]。

代码审查提示：`docs/q1/heat_conduction.md` §9 当前写成“单元中心有限体积法”，与代码不一致。正式论文应改为“顶点型/节点中心有限体积法”，除非后续确实重构为格心网格。

#### 4.2.2 Crank–Nicolson 与 Rannacher 启动

实际代码位置：`solve_heat` 和 `solve_mass` 中默认 `theta=0.5`；最初两个整步分别拆成两个半步后向欧拉。

- Crank–Nicolson 原始方法引用 R22；
- Rannacher 启动的收敛和抑制非光滑初始数据影响引用 R14。

建议表述：

> 主时间推进采用 \(\theta=1/2\) 的 Crank–Nicolson 格式 [R22]。由于初始均匀场与时变 Robin 边界可能在起始时刻产生较陡响应，前两个时间步各拆为两个后向欧拉半步，即采用 Rannacher 启动，以减弱早期数值振荡 [R14]。

不要只写“采用无条件稳定的 Crank–Nicolson 方法”就结束。稳定性不等于无振荡，也不代替时间步收敛检验。

#### 4.2.3 水分非线性 Picard 迭代

实际代码位置：`solve_mass` 约第 433–465 行。每次迭代用 \((C^n+C^{(m)})/2\) 更新扩散系数，再解冻结系数的线性三对角系统。

R24 可支撑固定点/Picard 迭代的一般数值基础。论文仍须报告本项目的收敛容差、最大迭代次数和实际最大迭代数；这些运行参数不能由文献替代。

建议表述：

> 对含水率相关扩散系数采用 Picard 固定点迭代：在每次迭代中冻结当前估计对应的界面扩散系数，求解线性化系统，直至相邻迭代的无穷范数变化低于给定容差 [R24]。

#### 4.2.4 调和平均的依据和适用边界

实际代码位置：`solve_mass.build` 中使用相邻节点扩散系数的调和平均。

R26 可用于说明：在分段常数或强非均匀扩散介质中，调和平均对应串联阻力，是经典界面处理。但必须同时写明本项目的限制：当 \(D(C)\) 在一个单元内随未知量指数退化时，点值调和平均可能被端点极小值支配。Q1 的短时高含水率区可将其作为基线；Q2–Q4 的低含水率正式结果不能仅凭 R26 继续使用调和平均。

#### 4.2.5 三对角求解器

Q1 自行实现 Thomas 消去，Q2–Q4 使用 `scipy.linalg.solve_banded`。这只是离散线性系统的求解实现，竞赛论文通常写“利用三对角追赶法/带状矩阵直接求解”即可，不必单独引用论文。

如果正文详细展开追赶法递推公式，可选用标准数值线性代数教材作为来源；没有必要引用 SciPy 软件文档作为论文的核心学术文献。代码层面应继续保留与稠密求解器对拍的单元测试。

#### 4.2.6 圆柱 Bessel 解析解基准

实际代码位置：`run_bessel_benchmark` 约第 732 行。代码在恒定环境温度和 Robin 换热边界下，求解特征方程

\[
\lambda J_1(\lambda)-\mathrm{Bi}\,J_0(\lambda)=0,
\]

并用 Bessel 级数与有限体积温度场对照。

R23 可作为圆柱瞬态导热解析解和 Bessel 展开的标准来源。论文应说明该解析解只验证常物性热方程、圆柱几何和 Robin 边界的离散，不验证水分方程、变物性耦合或长期边界外推。

### 4.3 Q2 代码：`src/q2_coupled.py`

#### 4.3.1 后向欧拉与时间中点物性

实际代码使用后向欧拉储存项和新时间层边界值；\(D\)、\(\rho c_p\) 与 \(k\) 则在旧状态和当前 Picard 状态的中点估计。论文应将它描述为：

> 后向欧拉时间推进下的半隐式/冻结系数 Picard 线性化。

不能因为系数取了时间中点，就把整体格式直接称为二阶 Crank–Nicolson。主时间离散仍具有后向欧拉的一阶特征，精度应以减半步长实验判断。

后向欧拉是标准方法，通常无需单独引用；若与 Q1 对比时间格式，可同时引用 R22 和通用数值分析教材。非线性固定点迭代引 R24。

#### 4.3.2 分块交替耦合

代码在每个 Picard 迭代内先更新水分，再根据新含水率更新热物性并求温度，是一种分块顺序耦合，而不是一次性求解完整单体 Jacobian。

论文建议写：

> 每个时间层采用分块 Picard 迭代，依次求解水分子问题和温度子问题，直至两个场的最大更新量同时满足容差 [R24]。

R01–R03 支撑物理耦合，R24 支撑非线性迭代；不要引用干燥综述来证明具体迭代一定收敛。收敛性应由最大迭代数、残差和减小步长结果证明。

#### 4.3.3 当前调和平均尚不是最终离散

`q2_coupled.py` 文件头已经注明待改为 Kirchhoff 积分平均。Q2 前 3 h 尚未进入最强退化区，因此旧结果可能变化很小，但正式论文的方法章节应在修复并重算后再引用 R15–R16。

在修复前，不应把 R26 写成“调和平均对本问题始终正确”的依据。

#### 4.3.4 回归、极值原理与守恒检查

代码包含：退化到 Q1 常物性的回归测试、网格/时间步比较、极值范围检查、中心—表面径向序检查和逐步守恒残差。

这些属于**软件验证证据**，不需要分别配参考文献。可引用 R12 解释有限体积的守恒结构，引用 R25 规范离散误差报告；最终可信度仍来自本项目实际测试结果。

### 4.4 Q3 代码：`src/q3_drying_time.py`

#### 4.4.1 人工分段变步长

当前步长依次为 1、5、20、60 s，并使所有切换点和步长与 60 s 输出网格对齐。该调度是依据边界变化和扩散时间尺度设计的工程策略，不需要为四个具体数值寻找文献。

论文需要给出固定 10 s、20 s 或统一缩小步长的对照，证明终止时间稳定。若讨论刚性扩散问题的一般隐式时间积分，可引用 R27；但当前正式 Q3 代码并没有使用 BDF，不能写成“Q3 使用 BDF”。

#### 4.4.2 阈值检测不是连续事件求根

当前代码在每个已接受时间步后检查 `C.max() < 0.15`，并在 60 s 网格上首次越过时停止。它属于离散步点阈值检测，没有在两个时间层之间进行插值、二分或事件求根。

因此论文应报告“60 s 输出网格上的首次达标时刻”，或在后续实现连续事件定位后再称为连续事件时间。该逻辑不需要外部引用，但需要前后两个输出时刻的夹逼证据。

#### 4.4.3 Richardson 外推和网格不确定度

R25 可作为网格收敛、观测阶和离散不确定度报告的规范来源。使用 Richardson 外推前必须满足：

1. 三组网格求解的是同一个离散模型；
2. 网格加密比明确；
3. 解进入渐近收敛区；
4. 时间误差相对空间误差足够小，或二者分别检验；
5. 被比较量应尽可能使用连续事件定位，而非被 60 s 量化后的值。

当前 Q3 的 `N=200/400/800` 外推建立在点值调和平均上，而该格式在强退化区存在结构性缺陷。因此旧 Richardson 数字只能作为历史诊断，不能因引用 R25 就升级为最终精度证据。应在 Kirchhoff 修正后重新计算。

#### 4.4.4 独立 BDF 复核

`docs/q3/drying_time_independent_review.md` 记录了格心有限体积加 SciPy BDF 的独立复核，但该 BDF 求解器当前不在 `src/` 或 `tests/` 中。R27 可说明 BDF 适用于刚性常微分方程组，SciPy 官方文档可作为软件实现说明；论文中更重要的是强调它与正式代码在网格、时间推进和事件定位上相互独立。

若要把该复核写成可重复的正式验证，应先把求解器整理进 `tests/`，固定容差和 Jacobian 稀疏结构，再引用 R27。没有可运行代码时，不要展开过多实现细节。

#### 4.4.5 代码与项目记忆不一致

当前 Q3 仍存在三项已知过渡状态：

1. 仍调用 `harm(Dm)`，尚未采用项目已冻结的 Kirchhoff 正式方案；
2. 文件头仍把薄层称为真实“硬壳层”，并声称 \(N=800\) 必需；该解释已被独立复核否定；
3. Picard 达到最大迭代次数后没有像 Q2、Q4 那样显式告警。

这些是实现/文档一致性问题，不是增加引用能解决的问题。正式论文引用前必须先完成修复与重算。

### 4.5 Q4 代码：`src/q4_moving_boundary.py`

#### 4.5.1 材料坐标和移动域离散

代码以 \(\xi=r/R(t)\) 把移动区域映射到固定区间，在新时间层使用 \(R(t_{n+1})^2W_i\) 作为储存权重，并以 \(hR(t_{n+1})\)、\(h_mR(t_{n+1})\) 组装表面项。

R19 可支撑食品收缩干燥的移动边界建模，R21 可支撑 ALE/随动网格的一般路线。本文的材料坐标抵消关系仍应由自己的链式法则推导，不要用文献标题代替推导。

建议表述：

> 将物理域 \(0\le r\le R(t)\) 映射到固定材料坐标 \(\xi=r/R(t)\)。在均匀相似收缩假设下，材料运动项与坐标移动项抵消，收缩通过储存尺度 \(R^2\) 和表面尺度 \(R\) 进入有限体积方程 [R19, R21]。

#### 4.5.2 Kirchhoff 势和解析 `exp1` 实现

代码不使用查找表，而是计算

\[
\Psi(C)=C\exp(-a/C)-aE_1(a/C),
\]

再以割线斜率构造界面扩散系数。R15–R16 支撑 Kirchhoff 变换和退化扩散离散；R28 可作为指数积分 \(E_1\) 定义与恒等式的权威数学来源。

建议表述：

> 对含水率指数退化的扩散系数，引入 Kirchhoff 势并以相邻节点势差直接构造界面通量 [R15–R16]。针对题给指数形式，势函数可由指数积分 \(E_1\) 写成解析表达式，其定义和微分关系见 R28。

解析式正确性还应由数值微分检查 \(\Psi'(C)=\exp(-a/C)\) 验证；引用不能替代这一单元测试。

代码审查还显示：`docs/q4/moving_boundary_heat_moisture_derivation.md` §16.3 主要推导了格心网格下“最后一个单元中心—真实表面”的半网格非线性 Robin 重构，而当前代码采用顶点型网格，\(\xi=1\) 本身就是表面节点，并直接把 Robin 通量装入该节点的半控制体。这是两种不同的边界离散。论文若以当前代码为准，应描述后者，不能声称程序实际求解了 §16.3 的半网格非线性表面根；若决定采用 §16.3，则需要先修改代码并重新验证。

#### 4.5.3 Q4 的时间推进和事件时间

Q4 采用后向欧拉、Picard 交替和分段变步长；当 \(C_{\max}<1.01C_\mathrm{target}\) 时将步长缩小到 1 s。当前 `t_f` 是**秒级离散越过时刻**，不是严格连续求根得到的事件时间。

论文更稳妥的写法是“阈值附近以 1 s 步长定位，事件时间误差受最后跨越区间约束”。若确实需要称为连续事件时间，应在相邻时间层间增加二分推进或其他事件定位，并做 1、0.5、0.25 s 收敛比较。

#### 4.5.4 移动域守恒、极限测试和效应分解

代码用材料权重 \(2\sum_iW_iC_i\) 检查水分库存，并设置：

- \(h_m=0\) 时平均干基含水率不变；
- \(h=h_m=0\) 时均匀初态保持不变；
- 三档空间网格与减半时间步；
- 平台均值/附件末值的边界敏感性；
- 附录 3 固定半径、附录 4 固定半径、附录 4 收缩半径三组对照。

这些设计本身不需要逐项引用。R12 支撑守恒离散的一般结构，R25 支撑网格/时间离散误差报告，R18–R20 支撑为何要控制物性后单独识别收缩效应。

### 4.6 哪些代码细节通常不放学术引用

以下内容写清楚即可，不建议占用参考文献名额：

- NumPy/Pandas/OpenPyXL 的数组、表格和工作簿读写；
- 分段线性插值本身；
- 输出文件命名、目录结构和 macOS 大小写冲突检查；
- 四位小数保守导出阈值；
- 三对角矩阵的具体存储格式；
- `assert`、回读验收和数据完整性检查；
- 输出时刻精确落在时间步上的整数设计。

如果论文附录专门列软件环境，可以引用 NumPy、SciPy 的软件论文或官方文档，但它们不应替代数值方法的理论来源。

## 5. 论文各章节的最小引用配置

若比赛论文篇幅紧张，正文至少保留以下引用组合：

### 5.1 引言

- 药材/植物干燥背景：R04 或 R05；
- 热质耦合建模综述：R01、R02。

### 5.2 模型建立

- 有效 Fick 扩散和圆柱模型：R03、R06；
- 对流传质边界：R08；
- 气固平衡含水率限制：R10。

### 5.3 数值方法

- 有限体积法：R12；
- Kirchhoff 变换与退化扩散：R15、R16；
- Rannacher 启动仅在正文确实写到时引用 R14。

### 5.4 收缩模型

- 收缩综述：R18；
- 移动边界模型：R19；
- 收缩影响实例：R20。

### 5.5 模型评价

- 潜热和多相机理限制：R11；
- 真实表面硬化与当前模型能力边界：R17。

按这一最小配置，正文约使用 12–14 篇外部文献即可建立完整证据链。其余文献可放入附录、方法文档或作为备选。

## 6. 不应采用的引用方式

1. **不要用外部论文证明题面经验式正确。** 附录 3、4 的 \(\rho,c_p,k,D\) 公式属于比赛给定条件。
2. **不要引用其他材料的扩散率数值作为本药材参数。** 不同植物组织、温度和含水率基准下的 \(D_\mathrm{eff}\) 不可直接移植。
3. **不要用表面硬化文献证明旧 Q3 停滞结果。** R17 只证明真实现象存在，不证明当前简化模型已捕捉该现象。
4. **不要用有限体积文献代替网格收敛。** R12–R13 说明方法基础，不保证当前实现自动准确。
5. **不要声称文献证明潜热可忽略。** 当前无潜热口径来自题设闭合与参数可识别性，而不是普遍物理规律。
6. **不要把空气含湿量和固体干基含水率写成同一物理量。** 即使单位都显示为 kg/kg，分母和热力学含义仍不同。
7. **不要给项目计算结果配无关外部论文。** 终止时间、误差和敏感性应回指本项目表格、代码与复核。

## 7. 推荐参考文献表

文献键用于本指南内部映射。正式论文应按学校或竞赛模板重新编号，并统一采用 GB/T 7714 或指定格式。

### 7.1 干燥建模与药材背景

**R01** Akter F, Muhury R, Sultana A, Deb U K. A comprehensive review of mathematical modeling for drying processes of fruits and vegetables[J]. *International Journal of Food Science*, 2022, 2022. DOI: [10.1155/2022/6195257](https://doi.org/10.1155/2022/6195257). [Consensus record](https://consensus.app/papers/a-comprehensive-review-of-mathematical-modeling-for-akter-muhury/bb5376e787a55ea19fea3ec05e376794/?utm_source=chatgpt).

**R02** Perré P. The proper use of mass diffusion equations in drying modeling: Introducing the Drying Intensity Number[J]. *Drying Technology*, 2015, 33: 1949–1962. DOI: [10.1080/07373937.2015.1076836](https://doi.org/10.1080/07373937.2015.1076836). [Consensus record](https://consensus.app/papers/the-proper-use-of-mass-diffusion-equations-in-drying-perré/61a3ee2bbdaa5ecc894d1c921c1c8d95/?utm_source=chatgpt).

**R03** Maroulis Z, Kiranoudis C, Marinos-Kouris D. Heat and mass transfer modeling in air drying of foods[J]. *Journal of Food Engineering*, 1995, 26: 113–130. DOI: [10.1016/0260-8774(94)00040-G](https://doi.org/10.1016/0260-8774(94)00040-G). [Consensus record](https://consensus.app/papers/heat-and-mass-transfer-modeling-in-air-drying-of-foods-maroulis-kiranoudis/0ce825f9876e5965adb4f9e1c345ca43/?utm_source=chatgpt).

**R04** Guo H, Chen Y, Xu W, et al. Assessment of drying kinetics, textural and aroma attributes of *Mentha haplocalyx* leaves during the hot air thin-layer drying process[J]. *Foods*, 2022, 11: 784. DOI: [10.3390/foods11060784](https://doi.org/10.3390/foods11060784). [Consensus record](https://consensus.app/papers/assessment-of-drying-kinetics-textural-and-aroma-guo-chen/c3fe5406c843552ab483553b6514199b/?utm_source=chatgpt).

**R05** Thamkaew G, Galindo F. A review of drying methods for improving the quality of dried herbs[J]. *Critical Reviews in Food Science and Nutrition*, 2020, 61: 1763–1786. DOI: [10.1080/10408398.2020.1765309](https://doi.org/10.1080/10408398.2020.1765309). [Consensus record](https://consensus.app/papers/a-review-of-drying-methods-for-improving-the-quality-of-thamkaew-galindo/59b07237f24b5e1596b6c6cc3f7489c5/?utm_source=chatgpt).

### 7.2 圆柱扩散、边界条件与平衡含水率

**R06** Faggion H, Tussolini L, Freire F, Freire J, Zanoelo E. Mechanisms of heat and mass transfer during drying of mate (*Ilex paraguariensis*) twigs[J]. *Drying Technology*, 2016, 34: 474–482. DOI: [10.1080/07373937.2015.1060498](https://doi.org/10.1080/07373937.2015.1060498). [Consensus record](https://consensus.app/papers/mechanisms-of-heat-and-mass-transfer-during-drying-of-mate-faggion-tussolini/41798dc0034c54b59b0c5f32b6b6da50/?utm_source=chatgpt).

**R07** da Silva W D, e Silva C M D P S, Farias V S O, Gomes J. Diffusion models to describe the drying process of peeled bananas: Optimization and simulation[J]. *Drying Technology*, 2012, 30: 164–174. DOI: [10.1080/07373937.2011.628554](https://doi.org/10.1080/07373937.2011.628554). [Consensus record](https://consensus.app/papers/diffusion-models-to-describe-the-drying-process-of-peeled-silva-silva/5e801fb2c02a577090f150b8ade27b44/?utm_source=chatgpt).

**R08** Kaya A, Aydin O, Demirtaş C. Concentration boundary conditions in the theoretical analysis of convective drying process[J]. *Journal of Food Process Engineering*, 2007, 30: 564–577. DOI: [10.1111/j.1745-4530.2007.00131.x](https://doi.org/10.1111/j.1745-4530.2007.00131.x). [Consensus record](https://consensus.app/papers/concentration-boundary-conditions-in-the-theoretical-kaya-aydin/9334b8f474ab54ee82f4c9a1138262d0/?utm_source=chatgpt).

**R09** Giner S A, Irigoyen R M T, Cicuttín S, Fiorentini C. The variable nature of Biot numbers in food drying[J]. *Journal of Food Engineering*, 2010, 101: 214–222. DOI: [10.1016/j.jfoodeng.2010.07.005](https://doi.org/10.1016/j.jfoodeng.2010.07.005). [Consensus record](https://consensus.app/papers/the-variable-nature-of-biot-numbers-in-food-drying-giner-irigoyen/9045a4d378805408b4e37e36c69b870a/?utm_source=chatgpt).

**R10** Soysal Y, Öztekin S. Equilibrium moisture content equations for some medicinal and aromatic plants[J]. *Journal of Agricultural Engineering Research*, 1999, 74: 317–324. DOI: [10.1006/jaer.1999.0463](https://doi.org/10.1006/jaer.1999.0463). [Consensus record](https://consensus.app/papers/equilibrium-moisture-content-equations-for-some-soysal-öztekin/3d5d17915dbd535cb93fb6f609d32d26/?utm_source=chatgpt).

**R11** Purlis E. Modelling convective drying of foods: A multiphase porous media model considering heat of sorption[J]. *Journal of Food Engineering*, 2019, 263: 132–146. DOI: [10.1016/j.jfoodeng.2019.05.028](https://doi.org/10.1016/j.jfoodeng.2019.05.028). [Consensus record](https://consensus.app/papers/modelling-convective-drying-of-foods-a-multiphase-porous-purlis/d9d810029c385e9d8ec5a84816f4e63c/?utm_source=chatgpt).

### 7.3 数值方法

**R12** Droniou J. Finite volume schemes for diffusion equations: Introduction to and review of modern methods[J]. *Mathematical Models and Methods in Applied Sciences*, 2014, 24: 1575–1619. DOI: [10.1142/S0218202514400041](https://doi.org/10.1142/S0218202514400041). [Consensus record](https://consensus.app/papers/finite-volume-schemes-for-diffusion-equations-droniou/55fefda1cc01578abb7c2bc0bd753bb1/?utm_source=chatgpt).

**R13** Herbin R, Vignal M. Error estimates on the approximate finite volume solution of convection diffusion equations with general boundary conditions[J]. *SIAM Journal on Numerical Analysis*, 2000, 37: 1935–1972. DOI: [10.1137/S0036142999351388](https://doi.org/10.1137/S0036142999351388). [Consensus record](https://consensus.app/papers/error-estimates-on-the-approximate-finite-volume-solution-herbin-vignal/e702706a62e252c1a94db531db577cff/?utm_source=chatgpt).

**R14** Giles M, Carter R. Convergence analysis of Crank–Nicolson and Rannacher time-marching[J]. *Journal of Computational Finance*, 2006: 89–112. DOI: [10.21314/JCF.2006.152](https://doi.org/10.21314/JCF.2006.152). [Consensus record](https://consensus.app/papers/convergence-analysis-of-cranknicolson-and-rannacher-giles-carter/a92b6a784db454129f721617142ce4c5/?utm_source=chatgpt).

**R15** Gama R. The Kirchhoff transformation and Fick’s second law with concentration-dependent diffusion coefficient[J]. *WSEAS Transactions on Heat and Mass Transfer*, 2021, 16. DOI: [10.37394/232012.2021.16.9](https://doi.org/10.37394/232012.2021.16.9). [Consensus record](https://consensus.app/papers/the-kirchhoff-transformation-and-the-fick’s-second-law-gama-gama/eb82de69940057feb563717aa4e53928/?utm_source=chatgpt).

**R16** Arbogast T, Huang C S, Zhao X. Finite volume WENO schemes for nonlinear parabolic problems with degenerate diffusion on non-uniform meshes[J]. *Journal of Computational Physics*, 2019, 399. DOI: [10.1016/j.jcp.2019.108921](https://doi.org/10.1016/j.jcp.2019.108921). [Consensus record](https://consensus.app/papers/finite-volume-weno-schemes-for-nonlinear-parabolic-arbogast-huang/df9387dace8057e4a585735fc3db4038/?utm_source=chatgpt).

### 7.4 表面硬化、收缩与移动边界

**R17** Gulati T, Datta A. Mechanistic understanding of case-hardening and texture development during drying of food materials[J]. *Journal of Food Engineering*, 2015, 166: 119–138. DOI: [10.1016/j.jfoodeng.2015.05.031](https://doi.org/10.1016/j.jfoodeng.2015.05.031). [Consensus record](https://consensus.app/papers/mechanistic-understanding-of-casehardening-and-texture-gulati-datta/d689d8a969575bf6a69de752ced6eee1/?utm_source=chatgpt).

**R18** Mayor L, Sereno A. Modelling shrinkage during convective drying of food materials: A review[J]. *Journal of Food Engineering*, 2004, 61: 373–386. DOI: [10.1016/S0260-8774(03)00144-4](https://doi.org/10.1016/S0260-8774(03)00144-4). [Consensus record](https://consensus.app/papers/modelling-shrinkage-during-convective-drying-of-food-mayor-sereno/0db6b7c884d05802ac5abbda910fa565/?utm_source=chatgpt).

**R19** Adrover A, Brasiello A, Ponso G. A moving boundary model for food isothermal drying and shrinkage: General setting[J]. *Journal of Food Engineering*, 2019. DOI: [10.1016/j.jfoodeng.2018.09.018](https://doi.org/10.1016/j.jfoodeng.2018.09.018). [Consensus record](https://consensus.app/papers/a-moving-boundary-model-for-food-isothermal-drying-and-adrover-brasiello/b3d191916ab155f3ae2aca95f881cbcc/?utm_source=chatgpt).

**R20** Aprajeeta J, Gopirajah R, Anandharamakrishnan C. Shrinkage and porosity effects on heat and mass transfer during potato drying[J]. *Journal of Food Engineering*, 2014, 144: 119–128. DOI: [10.1016/j.jfoodeng.2014.08.004](https://doi.org/10.1016/j.jfoodeng.2014.08.004). [Consensus record](https://consensus.app/papers/shrinkage-and-porosity-effects-on-heat-and-mass-transfer-aprajeeta-gopirajah/a59e0aab613153afa44302e8a086d04e/?utm_source=chatgpt).

**R21** Das R, Prasad K. Finite element modeling in heat and mass transfer of potato slice dehydration, nonisotropic shrinkage kinetics using arbitrary Lagrangian–Eulerian algorithm and artificial neural network[J]. *Journal of Food Process Engineering*, 2024. DOI: [10.1111/jfpe.14545](https://doi.org/10.1111/jfpe.14545). [Consensus record](https://consensus.app/papers/finite-element-modeling-in-heat-and-mass-transfer-of-potato-das-prasad/1ca14c16762658ffabc1ad1ce7fce927/?utm_source=chatgpt).

### 7.5 代码数值方法补充文献

**R22** Crank J, Nicolson P. A practical method for numerical evaluation of solutions of partial differential equations of the heat-conduction type[J]. *Mathematical Proceedings of the Cambridge Philosophical Society*, 1947, 43(1): 50–67. DOI: [10.1017/S0305004100023197](https://doi.org/10.1017/S0305004100023197).

**R23** Özisik M N. *Heat Conduction* [M]. 2nd ed. New York: John Wiley & Sons, 1993. 用于圆柱瞬态导热、第三类边界、Bessel 特征函数展开和解析基准。

**R24** Kelley C T. *Iterative Methods for Linear and Nonlinear Equations* [M]. Philadelphia: Society for Industrial and Applied Mathematics, 1995. DOI: [10.1137/1.9781611970944](https://doi.org/10.1137/1.9781611970944).

**R25** Celik I B, Ghia U, Roache P J, Freitas C J, Coleman H, Raad P E. Procedure for estimation and reporting of uncertainty due to discretization in CFD applications[J]. *Journal of Fluids Engineering*, 2008, 130(7): 078001. DOI: [10.1115/1.2960953](https://doi.org/10.1115/1.2960953).

**R26** Shashkov M, Steinberg S. Solving diffusion equations with rough coefficients in rough grids[J]. *Journal of Computational Physics*, 1996, 129(2): 383–405. DOI: [10.1006/jcph.1996.0257](https://doi.org/10.1006/jcph.1996.0257).

**R27** Hairer E, Wanner G. *Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems* [M]. 2nd ed. Berlin: Springer, 1996. DOI: [10.1007/978-3-642-05221-7](https://doi.org/10.1007/978-3-642-05221-7). 仅用于 Q3 独立 BDF 复核或刚性积分讨论；当前正式 Q3/Q4 代码使用后向欧拉，不能把 BDF 写成主求解器方法。

**R28** NIST Digital Library of Mathematical Functions. Exponential and logarithmic integrals: definitions and integral representations[EB/OL]. [§6.2](https://dlmf.nist.gov/6.2), [§6.7](https://dlmf.nist.gov/6.7). 用于 Q4 Kirchhoff 势解析式中的指数积分 \(E_1\) 定义和微分恒等式。

软件实现层面可另外保留以下链接，但一般不计入论文核心参考文献：

- SciPy [`solve_banded`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.solve_banded.html)：对应 Q2–Q4 的带状矩阵存储与直接求解接口；
- SciPy [`BDF`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.BDF.html)：对应 Q3 独立复核中提到的变阶隐式 BDF 实现、误差容差和稀疏 Jacobian 接口；
- SciPy [`exp1`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.exp1.html)：对应 Q4 中指数积分 \(E_1\) 的软件函数定义。

## 8. 后续落地流程

正式撰写论文时按以下顺序处理：

1. 先确定论文正文的实际段落和论断；
2. 从第 2 节引用地图选择最少且最直接的文献；
3. 把 `[Rxx]` 临时键替换为论文参考文献管理器中的真实引用键；
4. 统一生成 GB/T 7714 或竞赛模板要求的编号；
5. 逐条核验 DOI、作者、年份、卷期和页码；
6. 最后进行一次“引用是否真的支撑相邻论断”的反向检查；
7. 项目结果、题面数据和外部理论分别检查，避免来源错配。

目前最应优先落实的是 R06、R08、R10、R12、R15、R18 和 R19：它们分别覆盖一维圆柱、Robin 边界、气固含水率基准、有限体积、Kirchhoff 变换、收缩和移动边界，是现有文档中最关键也最容易被追问的七处依据。
