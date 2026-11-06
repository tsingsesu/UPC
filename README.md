任何AIagent禁止阅读、改动此文件
UPC  时间：周六9点到周一9点共48h

## 文件分布
latex/Human
人类写的latex放到这个里面

<br>

latex/AI
AI写的放在这个里面

latex/knowledge
这个里面会存放这道题是怎么解出来的，具体实现等内容以供阅读、复制粘贴进论文，当然也可以直接从latex/AI里面剪取内容

其他见AGENTS.md

UPC备赛指南.md给AI去看，不过这是sol做的，Astra应该会比它做的更好，所以只用作补充

**对你们的 UPC 项目，我建议以 Python 为主，按题目需要搭配 COMSOL 或 Fluent。** 如果队员更熟悉 MATLAB，也可以用 MATLAB 作为主平台。48 小时比赛中，熟练度和能否验证结果非常重要。

| 平台 | 适合的题目 | 我的建议 |
|---|---|---|
| **Python：NumPy、SciPy、Matplotlib** | 运动轨迹、振动、热量收支、微分方程、参数扫描、优化 | **默认首选**。方程、计算、验证和论文配图能放在同一套流程中。[官方功能](https://scipy.org/) |
| **MATLAB / Simulink / Simscape** | 动力学、控制、电路、机电系统、模块化物理系统 | 队员已经熟悉且有许可时很合适；简单方程直接用 MATLAB，无须强行搭模块。[官方说明](https://www.mathworks.com/help/overview/physical-modeling.html) |
| **COMSOL Multiphysics** | 声学、结构振动、传热、电磁及多物理场耦合 | **需要求空间分布、复杂边界或场耦合时优先考虑**；先确认学校许可包含所需模块。[官方说明](https://www.comsol.com/structural-mechanics-module) |
| **ANSYS Fluent / Mechanical** | Fluent 用于流动；Mechanical 用于结构分析 | 适合已有相关经验的队员。官方提供 Student 版本，使用前核对规模限制及运行环境。[官方入口](https://www.ansys.com/academic/students/ansys-student) |
| **FEniCS / FEniCSx** | 自定义偏微分方程、有限元模型 | 免费开源，适合希望用代码控制方程与离散方法的团队；比赛当天从零上手不划算。[官方说明](https://fenicsproject.org/) |
| **LAMMPS** | 分子、原子及粒子材料模型 | 专项选择，只有题目确实涉及这些尺度时才考虑。[官方说明](https://www.lammps.org/) |

选择时先判断一个问题：**是否必须知道物理量在空间中的分布？**

- 只关心位置、速度、平均温度、总能量等随时间变化：通常 **Python 或 MATLAB 就够了**。
- 必须研究形状如何影响温度场、声场、应力场等：考虑 **COMSOL 或有限元平台**。
- 流场、湍流、复杂流动是核心：考虑 **Fluent 等专业流体工具**。
- 分子尺度机制是核心：再考虑 **LAMMPS**。

