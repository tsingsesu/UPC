# AI使用记录（尚不等同于完整赛事披露）

实际后端：DeepSeek deepseek-v4-flash，Claude Code 2.1.207；代理请求记录返回HTTP 200。Codex负责决策、审查、修正、独立验证和英文稿。

下列记录保留已知的委派输入与反馈。原始完整对话及代码初次输出尚未归集，不能宣称已经符合官方完整输出要求。

## 首次任务输入

读取src/Mai实施计划.md，严格按合同仅实现src/Mai_penalty.m与src/Mai执行反馈.md。不运行命令，不使用子代理。不读取README.md或latex/Human。

## 首次任务附带计划（原文）

# B题数值实现合同

Codex已选择B题。你必须使用现有Claude Code所接入的DeepSeek，不能转交子代理，不运行命令。只允许创建/修改 `src/Mai_penalty.m` 与 `src/Mai执行反馈.md`；禁止读README.md、latex/Human、提交文件夹、任何个人模板。读题目/题目.md和本计划即可。不要替Codex写论文，不编造运行结果。

实现一个MATLAB R2024a主函数 `result=Mai_penalty(mode)`，mode默认'full'，'quick'样本少；源码使用英文合法函数名（MATLAB入口限制），中文注释。Optimization及Statistics工具箱已验证存在。项目根由mfilename路径计算，无个人路径。创建data/派生数据及img（正式结果），只输出中文名的JSON/CSV/PNG/MAT文件。全输出进入同一个结构result，JSON便于论文引用。MATLAB主力，Python独立复核另由Codex负责。

物理模型：r'=v，v'=[0;0;-9.80665]-K*Cd*norm(u)*u+cross(b,u)，u=v-w，K=rho*pi*R^2/(2*m)，rho=1.225,m=.43,R=.11。b=[bx;by;bz]单位1/s，是有效Magnus系数向量，不可谎称可唯一恢复真实自旋。已知r0=[0;0;.11]，tobs=[.05,.1]，观测两行[1.5,.45,.44;2.96,.86,.74]。拟合theta=[vx vy vz Cd by bz]，bx=0默认；使用lsqnonlin有界加权最小二乘，bounds=[10,-20,-5,.05,-1.5,-1.5]至[50,20,25,.6,1.5,1.5]。这些全是建模情景范围，不是文献已测区间。可用多个初值核查局部解；记录退出标志、残差、边界命中。sigma=.01 m为独立各坐标高斯噪声假设，不能由6点数据估计；另做.002,.005,.02 m情景。r0暂精确，明确限制。

积分拟合阶段使用固定步RK4 h=.005 s，精确到观测时刻；预测阶段ode45+事件定位x=11（direction +1）或下降触地z=R（direction -1），不能把发射初始z=R当触地终止：从t=.1的状态开始事件积分。最大飞行时间2s，未到球门/先触地记录分类和NaN，不静默忽略。预测点是约束最大似然参数的plug-in点，不要称无条件唯一概率众数。另做球心x=11.11的全越线几何敏感性，中心可通行条件abs(y)<3.66-R,z<2.44-R,z>=R。

基线：1.重力弹道，v0由两次观测线性最小二乘求得，给拟合残差和预测；2.常加速度恰好插值v0=(4*r1-r2-3*r0)/(.1)，a=(r2-2*r1+r0)/(.05^2)，解析求x=11最早正实根；3.仅阻力Cd与v0拟合；4.完整有效Magnus模型。

不确定性：参数雅可比J是观测位置/参数（不要混淆已除sigma的J），covTheta=sigma^2*pinv(J'*J)，SVD报告奇异值/条件数；预测映射g(theta)=[yGoal,zGoal]的数值雅可比G，C=G*covTheta*G'。联合95%椭圆半轴sqrt(chi2inv(.95,2)*eig(C))，边际区间±1.959964 sqrt(diag(C))，明确只是局部高斯线性化且不处理边界。参数bootstrap：从最优模型生成观测+sigma randn，重新拟合、记录失败/边界/未入门，full=600 quick=40；给有效样本分位数、均值、中位数、标准差、越过球门比例（bootstrap稳定性频率，不是实际进球概率）。不能把筛选进球样本后的区间冒充无条件区间。输出所有样本CSV。

验证：1.RK4在Cd=b=0时与解析重力弹道比较；2.无重力、无Magnus时阻力一维解析v=v0/(1+K*Cd*v0*t),x=log(1+K*Cd*v0*t)/(K*Cd)比较（积分器支持g参数）；3.同一最优theta下h=.01,.005,.0025观测误差和ode45预测maxstep=.02,.01,.005比较；4.无风Magnus功率v dot cross(b,v)=0及阻力耗能率非正检查；5.独立合成轨迹测试：随机已知参数vx[24,34],vy[3,10],vz[4,9],Cd[.15,.4],by[-.15,.15],bz[-.6,.1]，高精度ode45生成真值，加入sigma=.01，反演并预测，full=150 quick=10；输出入门/触地分类、RMSE和95%椭圆覆盖（记录分母，边界和失败，不能称实验确认）。另对球门可达样本使用较低/较高自旋等覆盖。

模型偏差情景：分别固定bx=-1,+1；风w=[0,-3,0]和[0,3,0]；变Cd=Cd0*(speed/30)^p的p=-.5,+.5；b(t)=b0*exp(-t/tau)的tau=.5,1秒。每个情景必须重新拟合相同早期观测再外推，给残差、预测点、相对基准偏差、分类。这些是压力测试，非概率区间。逐个改变因素。

正式图由MATLAB画英文轴/图例、SI单位、无标题：轨迹与观测、球门预测椭圆与bootstrap散点、噪声水平与边际区间半宽、模型情景比较。图片名中文PNG，避免过多绘图。所有结果来自result和样本。随机种子20261107。

写清软件版本、算法步长/容差、种子、数值检查指标和耗时到result。每批计算fprintf进度。不允许只写骨架/TODO；若工具写文件失败说明。反馈只能说已实现未运行，运行由Codex做。尽量约500行以内清晰实现，科学正确优先。


## 执行方反馈文件（原文）

# Mai执行反馈

合同方：Codex → Claude Code（DeepSeek）。日期：2026-11-07。B题数值实现。

## 交付范围

按 `src/Mai实施计划.md` 合同，仅创建：

- `src/Mai_penalty.m`（主实现，MATLAB R2024a）
- `src/Mai执行反馈.md`（本文件）

未运行任何命令，未使用子代理，未读取 README.md、latex/Human、提交文件夹或任何个人模板。除 `题目/题目.md` 与实施计划外无其他输入。

## 本次实际完成的内容

`Mai_penalty(mode)` 已按合同实现为完整单文件函数，`mode` 默认 `'full'`，`'quick'` 降低 Bootstrap 与合成反演样本数。**尚未运行**，以下“实现”均指代码已写入，不含任何实测数值结果。

1. **路径**：由 `mfilename('fullpath')` 计算项目根，自动创建 `data/`、`img/`，不含个人绝对路径。
2. **物理模型**：`r'=v`，`v'=[0;0;-9.80665]-K*Cd*||u||*u+cross(b,u)`，`u=v-w`，`K=rho*pi*R^2/(2*m)`，`rho=1.225`、`m=0.43`、`R=0.11`；`b=[bx;by;bz]` 有效 Magnus 系数（1/s），`bx=0` 固定，代码与结论文本均写明不能唯一恢复真实自旋。
3. **拟合**：`theta=[vx vy vz Cd by bz]`，有界加权 `lsqnonlin`（trust-region-reflective，加权残差 `(pred-obs)/sigma`，sigma=0.01 m），上下界严格取合同值；多初值逐个记录退出标志、SSE、RMSE、残差、边界命中，最终取 SSE 最小解。拟合阶段固定 RK4 `h=0.005 s`，按步数精确记录到 0.05/0.10 s。
4. **预测**：从 `t=0.10` 的状态开始单独 `ode45` 事件积分，事件为 `x=11`（direction +1）与 `z=R`（direction -1），上限 2 s，按先发生事件分类 `goal/ground/cross_out/timeout`，未到球门或先触地明确记录，NaN 不静默忽略。另做球心全越线 `x=11.11` 几何敏感性。中心可通行条件 `abs(y)<3.66-R, z<2.44-R, z>=R`。预测点标注为约束最大似然参数的 plug-in 点。
5. **不确定性**：观测位置雅可比 `J`（显式说明未除 sigma），`covTheta=sigma^2*pinv(J'*J)`，SVD 奇异值与条件数；预测映射数值雅可比 `G`，`C=G*covTheta*G'`；联合 95% 椭圆半轴 `sqrt(chi2inv(.95,2)*eig(C))`，边际区间 `±1.959964*sqrt(diag(C))`，并写明仅为局部高斯线性化、不处理边界。噪声情景 `.002/.005/.010/.020 m`。
6. **Bootstrap**：从最优模型预测位置生成观测 + sigma·randn 重新拟合，full=600 / quick=40；记录失败、边界命中、入门状态，输出有效样本的分位数（2.5/50/97.5）、均值、中位数、标准差、进球稳定性频率（含成功分母与全部抽样分母两种比例），所有样本写入 CSV；文本明确这不是实际进球概率，未把筛选后的进球样本区间冒充无条件区间。
7. **基线**：重力弹道（v0 线性最小二乘）、常加速度精确插值（解析求 x=11 最早正实根）、仅阻力（Cd 与 v0）、完整有效 Magnus，四者均给出观测残差与预测。
8. **验证**：(a) Cd=b=0 时 RK4 与解析重力弹道比较；(b) 无重力无 Magnus 的阻力一维解析解比较（积分器支持 `g` 参数）；(c) `h=.01/.005/.0025` 与 `ode45 MaxStep=.02/.01/.005` 比较；(d) 无风时 `v·(b×v)` 与阻力耗能率非正检查；(e) 随机已知参数（范围内）高精度 ode45 合成轨迹，加 sigma=0.01 噪声后反演并预测，full=150 / quick=10，输出分类、参数 RMSE、预测 RMSE、95% 椭圆覆盖（记录分母、边界、失败），并对可达样本按自旋较低/较高分组统计覆盖。
9. **模型偏差情景**：固定 `bx=-1/+1`、风 `[0,±3,0]`、变 Cd（p=-0.5/+0.5）、`b(t)=b0*exp(-t/tau)`（tau=0.5/1.0 s），每个情景重新拟合同一早期观测后外推，给出残差、预测点、相对基准偏差、分类；逐个因素单独改变。
10. **图**：4 张正式 PNG，中文文件名，英文轴标/图例、SI 单位、无标题：轨迹与观测、球门预测椭圆与 Bootstrap 散点、噪声水平与边际区间半宽、模型情景比较。
11. **输出**：统一 `result` 结构（软件版本、算法步长/容差、种子 20261107、数值检查指标、耗时、全部模块结果与说明），写 `data/Mai罚球结果.json` 与 `.mat`；另有拟合残差、自助样本、合成反演样本、模型偏差情景、噪声情景 CSV。全部为中文文件名。
12. **其他**：各批计算 `fprintf` 进度；约 950 行（含中文注释与结构体字段说明），未留 TODO/骨架。

## 我明确未做 / 不能声称的

- **未运行**：无任何拟合值、RMSE、椭圆、覆盖率的实测数字；本反馈不含伪造结果。
- **未替 Codex 写论文**：不包含论文正文、结论表述或结果解释。
- **模型限制（代码内已注释）**：`bx=0` 不可识别真实自旋；`sigma` 为建模假设不可由 6 个观测数估计；`r0` 暂按精确已知处理；线性化区间不覆盖边界；Bootstrap 频率非进球概率；情景为压力测试。
- **工具箱依赖**：`lsqnonlin`（Optimization）、`chi2inv`/`prctile`（Statistics and Machine Learning），合同称已验证存在；`exportgraphics` 需 R2020a+（R2024a 满足）。

## 交给 Codex 的运行建议

在项目根下：

```matlab
cd C:\Users\lhl87\Desktop\UPC
result = Mai_penalty('quick');   % 先快速检查流程与文件写出
result = Mai_penalty('full');    % 正式结果，输出 data/ 与 img/
```

若首次运行报错，优先检查：`lsqnonlin` 是否可达、`ode45` 事件是否在 2 s 内触发、`exportgraphics` 图形句柄是否可离屏渲染。所有数值结论以 Codex 实际运行结果为准。


## 第二次任务输入

读取src/Mai实施计划.md，严格按新任务将英文论文完整译成中文LaTeX，只写指定两个文件，不运行命令，不使用子代理。

# 中文论文翻译任务
只允许写latex/AI/B题中文论文.tex和src/Mai执行反馈.md。不要修改其它文件，不运行命令，不用子代理，不读README.md、latex/Human。读取latex/AI/B题英文论文.tex以及data/派生数据/Mai罚球结果.json。
将英文论文完整翻译成中文论文，用单个ctexart文档（UTF8,fontset=windows），公式、数值、图表、段落结论、引用与英文保持一致，不缩写成摘要。保留队号TEAM NUMBER REQUIRED占位，说明中文稿仅供队友阅读。中文字体使用SimSun，数学保持LaTeX。采用xelatex，页眉队号，图片相对../../img/，图片文件名不变。正文和参考文献不能瞎编新结果。特别保留：600次bootstrap53次边界165次框外，不能把72.5%称进球概率；150合成测试146覆盖，147有效平面对，全部分母146/150，3次预测触地，RMSE约.727m。解释1cm噪声下drag-only不能拒绝，不能仅凭6参拟合6点称Magnus已确证。AI披露尚需完整输入输出归集，不宣称已可提交。
无需重新生成图片或数据。若数值有不一致以已运行JSON为准并反馈。完成后写反馈，说明只实现未编译。


## 中文论文首次输出（原文）

```latex
\documentclass[11pt,a4paper,UTF8,fontset=windows]{ctexart}
\usepackage[margin=25mm]{geometry}
\usepackage{amsmath,amssymb,booktabs,graphicx,fancyhdr,hyperref,float}
\usepackage{xeCJK}
\setCJKmainfont{SimSun}
\hypersetup{colorlinks=true,linkcolor=blue,urlcolor=blue,citecolor=blue}
\pagestyle{fancy}\fancyhf{}\fancyhead[C]{队伍：TEAM NUMBER REQUIRED}\fancyfoot[C]{\thepage}
\setlength{\headheight}{15pt}\setlength{\parskip}{5pt}\setlength{\parindent}{0pt}
\graphicspath{{../../img/}}
\newcommand{\R}{\mathbb R}
\begin{document}
\begin{center}
{\large 队伍：TEAM NUMBER REQUIRED}\par\vspace{12mm}
{\Large B 题}\par\vspace{10mm}
{\LARGE 由 0.10 秒运动预测点球}\par\vspace{8mm}
{\large 一种受物理约束的估计及其不确定性}
\end{center}
\vspace{6mm}
\begin{center}
{\small 中文译本仅供队友阅读；正式提交以英文稿为准，队号待定。}
\end{center}
\vspace{4mm}
\begin{abstract}
两次短时位置测量可以约束足球的速度与曲率，但在没有关于空气动力作用的假设时，无法唯一确定其未来运动。我们构造了一个有界最大似然估计器，采用重力、二次阻力以及一个有效马格努斯（Magnus）项。将拟合出的轨迹外推至其首次穿越球门平面，并把触地视为一个竞争事件。对于题给观测，名义模型预测在踢出后 $0.401$ s 时于约 $(11,2.60,1.97)$ m 处穿越。若改用常加速度插值，则给出 $(11,2.49,1.92)$ m，这说明了外推假设的影响。测量不确定性会被强烈放大：假设各坐标误差相互独立且为 1 cm，则两个球门平面坐标的局部 95\% 边际半宽均约为 0.97 m。我们把这种条件抽样不确定性，与由风、未观测到的自旋分量、速度相关阻力以及自旋衰减引起的结构不确定性区分开来。数值验证采用了解析极限情形、步长细化、独立积分以及合成预测检验。由于没有独立的实验轨迹可用，这些检验确立的是数值一致性而非实验有效性。因此，一个有实用价值守门员预测应当是一个点、一个不确定性区域以及一种失效分类；仅凭题给信息，无法支撑一个高度精确的无条件落点。
\end{abstract}
\newpage
\section{问题理解与交付内容}
我们采用 B 题~\cite{upc} 所规定的坐标与观测。球的初始质心为 $\mathbf r_0=(0,0,0.11)$ m；随后的两次观测为
\begin{equation}
\mathbf r_1=(1.50,0.45,0.44)\ \mathrm m,\qquad
\mathbf r_2=(2.96,0.86,0.74)\ \mathrm m,
\end{equation}
分别对应 $t_1=0.05$ s 与 $t_2=0.10$ s。我们的主要预测是球心首次向前穿过 $x=11$ m 时的位置，遵循题目对平面的约定。另有一个单独考虑有限半径的检验，用来评估球是否从球门口通过。本模型不预测守门员的反应，也不预测球与门柱的碰撞。

所要求的一般方法必须能接受实测位置、估计轨迹、传播不确定性，并拒绝那些先触地或未能到达球门平面的轨迹。我们逐步显式实现了每一环节。所有报告的不确定性都以所述噪声与受力假设为条件；而这两者题目均未完全给定。

\section{信息极限与解析基线}
\subsection{为什么插值不等于验证}
在初始位置固定之后，数据只包含六个标量约束。一个足够灵活的六参数模型可以精确插值这些观测，因而不留下任何可用来估计噪声的残差自由度。包括三维自旋、可变阻力、风以及相机偏置在内的其它未知量，都无法全部自由地反演出来。

为做一个简单的构造性说明，向一条光滑轨迹加入 $\mathbf c\,t(t-t_1)(t-t_2)$ 可以保持全部三个观测位置不变，却改变此后的位置。物理限制会约束允许的 $\mathbf c$；而仅凭观测无法做到这一点。因此我们采用一个降阶的物理模型，并另行检验偏离该模型的情形。

\subsection{常加速度}
取 $h=0.05$ s 与 $\mathbf r(t)=\mathbf r_0+\mathbf v_0t+\tfrac12\mathbf at^2$，精确插值得
\begin{equation}
\mathbf v_0=\frac{4\mathbf r_1-\mathbf r_2-3\mathbf r_0}{2h}
=(30.4,9.4,6.9)\ \mathrm{m\,s^{-1}},\qquad
\mathbf a=\frac{\mathbf r_2-2\mathbf r_1+\mathbf r_0}{h^2}
=(-16,-16,-12)\ \mathrm{m\,s^{-2}}.
\end{equation}
方程 $30.4t-8t^2=11$ 的最小正解为 $0.405008$ s，给出 $(y,z)=(2.494824,1.920367)$ m。保留这些多余位数仅是为了可复现性检查。该加速度是推算出的短时平均值，而不是一个保证始终不变的作用力定律。

仅重力基线在扣除已知的重力位移后，用线性最小二乘估计初速度。仅阻力基线则估计三个速度分量与一个阻力系数。其残差与预测将在数值结果中与完整模型作比较。

\subsection{短观测窗口放大误差}
假设 $\mathbf r_0$ 精确，且 $\mathbf r_1,\mathbf r_2$ 的每个坐标均带有标准差为 $\sigma$ 的独立误差。加速度的不确定性为
\begin{equation}
\operatorname{sd}(a_j)=\frac{\sqrt5\,\sigma}{h^2}.
\end{equation}
当 $\sigma=0.01$ m 时，该值为 $8.94\ \mathrm{m\,s^{-2}}$。在固定时刻 $t=8h=0.40$ s，二次插值的权重为 $(21,-48,28)$，于是其坐标不确定性为 $\sqrt{48^2+28^2}\,\sigma=55.57\sigma$。球门平面预测还需要不确定的穿越时刻；仅在固定时刻传播不确定性是不充分的。

\section{物理模型与合理范围}
\subsection{方程与参数}
令 $\mathbf u=\mathbf v-\mathbf w$ 为相对空气的速度。牛顿第二定律近似为
\begin{align}
\dot{\mathbf r}&=\mathbf v,\\
\dot{\mathbf v}&=-g\mathbf e_z-K C_D|\mathbf u|\mathbf u+\mathbf b\times\mathbf u,
&K&=\frac{\rho\pi R^2}{2m}.
\label{eq:model}
\end{align}
阻力项遵循常规的二次作用力定律。风洞研究已经确立，足球的阻力与侧向力取决于表面结构、流态与旋转~\cite{asai}。这些研究支持把这些机理纳入模型，但并不能标定本题中的具体足球。

有效马格努斯向量 $\mathbf b$ 的量纲为每秒。若升力系数近似正比于自旋参数 $R\omega/|\mathbf u|$，则横向加速度具有这种叉积形式。因此 $\mathbf b$ 综合了自旋与空气动力效率；它并不是对角速度的测量。常数 $\mathbf b$ 与 $C_D$ 是短飞行近似，后文会对其施加压力测试，而不是把它们当作普遍定律。

\begin{table}[H]\centering
\caption{输入与假设。搜索边界是建模选择，而不是实测总体区间。}
\begin{tabular}{lll}\toprule
量 & 取值 & 性质\\\midrule
初始质心与观测 & 如题所给 & 题目数据\\
球半径 $R$ & $0.11$ m & 触地假设\\
球质量 $m$ & $0.43$ kg & 代表性选择\\
空气密度 $\rho$ & $1.225$ kg m$^{-3}$ & 环境假设\\
重力 $g$ & $9.80665$ m s$^{-2}$ & 标准代表值\\
阻力系数 $C_D$ & $[0.05,0.60]$ & 搜索域\\
$b_y,b_z$ & $[-1.5,1.5]$ s$^{-1}$ & 搜索域\\
坐标误差 $\sigma$ & $2,5,10,20$ mm & 相互独立的噪声情景\\
\bottomrule\end{tabular}\end{table}
所选质量与半径与 IFAB 关于球的规定相符，即 410--450 g、周长 68--70 cm~\cite{ifab}。密度与重力是假设而非现场测量。初始速度边界为 $v_x\in[10,50]$、$v_y\in[-20,20]$、$v_z\in[-5,25]$ m s$^{-1}$。

在参考模型中 $\mathbf w=\mathbf0$ 且 $b_x=0$，于是
\begin{equation}
\boldsymbol\theta=(v_{0x},v_{0y},v_{0z},C_D,b_y,b_z).
\end{equation}
固定 $b_x$ 是为了削弱一个弱可辨识方向；这并不能证明真实自旋的纵向分量为零。若三个自旋分量都自由，则由六个观测将得出至少七个未知量。同样，轨迹数据主要识别的是乘积 $KC_D$，因此独立地改变质量、密度与拟合阻力，会对同一种不确定性重复计算。

\subsection{尺度与守恒检验}
拟合得到的初速度约为 $32.6$ m s$^{-1}$。取假设的粘度为 $1.81\times10^{-5}$ Pa s，可得 $\mathrm{Re}=2\rho R|\mathbf v|/\mu\simeq4.9\times10^5$。这是对采用空气动力模型的合理性检验，而不是对其系数的标定。在 $T\simeq0.40$ s 内，无量纲变化量 $KC_DUT\simeq0.26$、$|\mathbf b|T\simeq0.15$、$gT/U\simeq0.12$ 表明，阻力、曲率与重力都可能实质性地影响外推。

在静止空气中，机械能满足
\begin{equation}
\frac{d}{dt}\left(\frac12m|\mathbf v|^2+mgz\right)
=-mKC_D|\mathbf v|^3\le0,
\end{equation}
因为 $\mathbf v\cdot(\mathbf b\times\mathbf v)=0$。这才是恰当的耗散平衡。当流动空气做功时，地面参考系中的机械能不必单调减少。

\section{推断、事件检测与不确定性}
\subsection{有界最大似然}
对选定的 $\sigma$，我们在所述边界内最小化
\begin{equation}
\chi^2(\boldsymbol\theta)=\sum_{i=1}^2
\frac{\|\mathbf r(t_i;\boldsymbol\theta)-\mathbf r_i\|^2}{\sigma^2},
\end{equation}
所用工具为 \texttt{lsqnonlin}~\cite{lsq}。各坐标方差相等意味着单独改变 $\sigma$ 只会缩放目标函数，而不改变内点最优解。多个起始点用于检查局部一致性；它们并不证明全局最优。

点预测是受约束的最大似然参数估计在球门平面上的像。正是在这个精确意义上，它是我们模型条件下的“最可能”预测。若要求边际落点概率密度的众数，还需要一个概率化的参数模型并对冗余参数积分；我们不声称给出这种无条件众数。

观测值使用四阶 Runge--Kutta 法、步长 $0.005$ s 进行求值。后续飞行使用带方向事件的 \texttt{ode45}~\cite{ode}：$x-11=0$ 递增、$z-R=0$ 递减。事件积分从 $0.10$ s 处的拟合状态开始，以免把起脚触地误判为立即落地。未发生穿越、数值失败以及更早触地都会被显式记录。适用性前提是球在观测窗口内处于空中。

题目给出的中心平面判据为 $|y|<3.66$、$0<z<2.44$。更严格的球门口净空判据为 $|y|<3.66-R$、$R\le z<2.44-R$。我们还计算了 $x=11+R$ 处的球心位置，作为有限尺寸敏感性分析；这只是一个理想薄平面近似，而不是详细的碰撞后模型。

\subsection{局部置信域}
令 $\mathbf f(\boldsymbol\theta)\in\R^6$ 收集预测观测值，$J=\partial\mathbf f/\partial\boldsymbol\theta$，$G=\partial(y_g,z_g)/\partial\boldsymbol\theta$，其中包含穿越时刻的变化。对于满秩内点解，
\begin{equation}
C_\theta\simeq\sigma^2(J^TJ)^{-1},\qquad C_g\simeq GC_\theta G^T.
\end{equation}
局部 95\% 联合椭圆为
\begin{equation}
(\mathbf q-\widehat{\mathbf q})^TC_g^{-1}(\mathbf q-\widehat{\mathbf q})\le5.991.
\end{equation}
边际坐标半宽为 $1.959964\sqrt{(C_g)_{jj}}$。它们不是联合 95\% 的矩形边界。在参数约束起作用时，高斯线性化只是近似，可能错误地表示覆盖率。

\subsection{重采样与模型不确定性}
参数化自助法由拟合模型生成带噪观测，对每个样本重新拟合，并传播每个成功拟合的结果。我们在报告坐标分位数的同时，也报告边界命中、失败、触地与平面穿越。这描述的是估计器的抽样稳定性。落在球门口内的自助预测比例既不是实验测定的得分概率，也不是贝叶斯后验概率。

结构检验每次只改变一个机理：固定 $b_x=\pm1$ s$^{-1}$、横向风 $w_y=\pm3$ m s$^{-1}$、$C_D=C_{D0}(|\mathbf u|/30)^p$ 且 $p=\pm0.5$、以及 $\mathbf b(t)=\mathbf b_0e^{-t/\tau}$ 且 $\tau=0.5,1$ s。每一种都针对同一组观测重新拟合。这些刻意选定的压力情形不赋予出现概率，也不得合并进统计置信区间。

\section{数值结果与验证}
% VERIFIED_RESULTS_BEGIN
\subsection{点预测与竞争基线}
四个独立起点收敛到同一个内点解，参数为
\begin{equation}
\widehat{\boldsymbol\theta}=(30.41546,9.41672,6.90801,0.370687,-0.060326,-0.357128).
\end{equation}
前三个分量单位为 m s$^{-1}$，第四个为无量纲量，后两个单位为 s$^{-1}$。预测的穿越点为 $(y_g,z_g)=(2.599441,1.967005)$ m，时刻 $0.4005035$ s。保留这种数值精度是为了可复现性；答案应报告为 $0.401$ s 时 $(11,2.60,1.97)$ m。
\begin{table}[H]\centering
\caption{观测拟合与外推。极小的插值残差不能证明物理有效性。}
\begin{tabular}{lrrrr}\toprule
模型 & RMSE (mm) & $t_g$ (s) & $y_g$ (m) & $z_g$ (m)\\\midrule
仅重力 & 10.376 & 0.3706 & 3.217 & 1.957 \\
常加速度 & 0.000 & 0.4050 & 2.495 & 1.920 \\
仅阻力 & 5.029 & 0.4058 & 3.217 & 1.879 \\
阻力与有效马格努斯 & 0.000 & 0.4005 & 2.599 & 1.967 \\
\bottomrule\end{tabular}\end{table}
仅阻力模型的目标函数值为 $\chi^2=1.517$，对应两个名义残差自由度；在 1 cm 高斯噪声假设下，这相当于近似拟合优度尾概率 $0.468$。因此，在该噪声水平下，这些观测本身不足以确立需要那两个额外横向力参数。在 2 mm 噪声下，同样的仅阻力偏差给出 $\chi^2=37.9$，说明解释结果强烈依赖测量质量。我们保留这个更灵活的模型，是把它作为一个条件性的、能感知曲率的预测，而不是因为它的精确拟合就宣称它唯一正确。

\begin{figure}[H]\centering
\includegraphics[width=.78\linewidth]{Mai罚球轨迹与观测.png}
\caption{从出手到首次穿越球门平面的名义拟合飞行轨迹。只有两个早期位置是观测值。}
\end{figure}
在中心平面 $x=11.11$ m 处，有限半径敏感性给出 $0.40493$ s 时 $(y,z)=(2.61689,1.97822)$ m。额外的行程使球心在球门平面内偏移约 2.1 cm。在两种判据下，名义射门都能通过理想球门口。

\subsection{测量不确定度与自助稳定性}
当独立坐标误差的标准差为 0.01 m 时，局部预测协方差为
\begin{equation}
C_g\simeq\begin{pmatrix}0.245132&0.000117\\0.000117&0.243555\end{pmatrix}\ \mathrm{m^2}.
\end{equation}
边际 95\% 半宽分别为 0.9704 m 与 0.9673 m。在 2、5、20 mm 噪声水平下，相应的半宽在每个坐标上分别缩放至约 0.194、0.485 与 1.94 m。这一缩放是局部结果；这并不是在断言有界非线性重采样在每个噪声水平上都精确按此缩放。
95\% 联合椭圆的半轴为 1.212 m 与 1.208 m。未缩放的观测雅可比条件数为 15.0；它依赖于参数单位，仅作为数值诊断报告，而不是不变的辨识度指标。

全部 600 次自助拟合均成功终止；其中 53 次至少命中一个参数边界。所有样本都到达球门平面，但有 165 次未满足有限半径的球门口净空判据。后者仍保留在坐标汇总与图中。重新估计的平面坐标的经验 2.5--97.5 百分位区间为：$y$ 向 $[1.703,3.617]$ m，$z$ 向 $[0.967,2.937]$ m。这些是估计器的重采样分位数，而不是真实射门的后验分布。72.5\% 的落入球门口频率仅是一个稳定性诊断量。边界命中以及延伸到画面之外的区域，都是结果的实质性特征。

\begin{figure}[H]\centering
\includegraphics[width=.88\linewidth]{Mai球门预测椭圆与自助散点.png}
\caption{在 1 cm 坐标噪声下由 600 次自助重拟合得到的球门平面预测。保留了全部平面穿越点，包括落在球门口之外的预测。椭圆为局部 95\% 联合置信近似。}
\end{figure}

\subsection{验证与合成预测检验}
主积分、优化与重采样由 MATLAB R2024a 完成，随机种子为 20261107。独立的 Python 实现使用 SciPy 1.17.1、DOP853 积分与有界非线性最小二乘。它与 MATLAB 的名义球门平面坐标一致到小于 $2\times10^{-11}$ m。这一极小差异确立的是实现一致性，而不是真实世界预测精度。
RK4 重力基准的最大位置误差为 $2.7\times10^{-15}$ m，一维二次阻力基准为 $5.2\times10^{-13}$ m。在拟合参数固定时，步长 0.01、0.005、0.0025 s 使预测穿越点的变化小于 $3\times10^{-10}$ m。独立地把事件积分器的最大步长改为 0.02、0.01、0.005 s，结果变化小于 $4\times10^{-12}$ m。事件积分的容差为 RelTol=$10^{-9}$ 与 AbsTol=$10^{-11}$。细化检验针对的是积分误差，而不是拟合物理中的不确定性。在舍入误差量级的误差下，我们不声称任何收敛阶。一项独立的积分能量检验把机械能损失与阻力耗散平衡到 $2\times10^{-13}$ J。

在合成检验中，150 个真值向量各自独立、均匀地抽取自 $v_x\in[24,34]$、$v_y\in[3,10]$、$v_z\in[4,9]$ m s$^{-1}$、$C_D\in[0.15,0.40]$、$b_y\in[-0.15,0.15]$、$b_z\in[-0.60,0.10]$ s$^{-1}$。这些范围定义的是一个设计好的检验集，而不是观测到的总体统计。高精度 ode45 轨迹独立于 RK4 拟合程序地生成早期观测，然后加入标准差为 1 cm 的噪声，每个案例从两个起点重新拟合。

全部 150 次反演均收敛，其中 31 次命中边界。全部真值轨迹都到达球门平面；但有 3 条拟合轨迹先触地。对 147 组双方都穿越平面的配对，双坐标 RMSE 为 0.727 m。局部 95\% 椭圆覆盖了这 147 个真实穿越点中的 146 个（99.32\%）；若把 3 个未穿越的预测计为未覆盖，则为 146/150（97.33\%）。这两个数字都不能证明精确的 95\% 标定，尤其是在只有 150 个案例且存在有效约束的情况下。球门口内/球门口外/触地的分类在 123/150 个案例中一致。合成数据服从参考力模型族，因此这些结果并未针对真实球飞行检验结构正确性。

\subsection{结构压力情景}
\begin{table}[H]\centering
\caption{单因子变化，每一种都针对同一组观测重新拟合。距离是相对名义穿越点的偏移；这些是情景，而不是置信边界。}
\begin{tabular}{lrrr}\toprule
情景 & $y_g$ (m) & $z_g$ (m) & 偏移 (cm)\\\midrule
$b_x=-1$ & 2.5516 & 2.0402 & 8.75 \\
$b_x=+1$ & 2.6601 & 1.9037 & 8.77 \\
风 $w_y=-3$ & 2.6004 & 1.9669 & 0.10 \\
风 $w_y=+3$ & 2.5983 & 1.9664 & 0.13 \\
$C_D$ 指数 $p=-0.5$ & 2.5983 & 1.9648 & 0.25 \\
$C_D$ 指数 $p=+0.5$ & 2.6005 & 1.9690 & 0.22 \\
$b$ 衰减 $\tau=0.5$ s & 2.7135 & 1.9503 & 11.53 \\
$b$ 衰减 $\tau=1.0$ s & 2.6611 & 1.9579 & 6.23 \\
\bottomrule\end{tabular}\end{table}
每个压力情形都把早期位置拟合到数值精度，但未来位置仍有差异。在所选情景中最大偏移约为 11.5 cm，出现在 0.5 s 衰减情形。定常横向风在重新拟合后的残余影响很小，因为自由的初速度与力参数吸收了它的大部分早期影响；这并不意味着此后突发的阵风是无害的。未观测到的纵向自旋造成约 9 cm 的偏移。压力范围仅作示例，并不能界定所有物理上可能的未来扰动。
\begin{figure}[H]\centering
\includegraphics[width=.9\linewidth]{Mai模型情景预测比较.png}
\caption{单因子重拟合后相对名义球门平面预测的位移。单位为厘米；所选情景不携带任何概率权重。}
\end{figure}

% VERIFIED_RESULTS_END

\section{实际含义、局限与改进}
名义球心位于球门的正 $y$ 一侧，且位置足够高，以至于垂直方向的不确定性有影响。在 0.10 s 时，距离预测穿越点仅剩约 0.30 s。输入小数位所暗示的精度不能与已知的相机精度混为一谈：仅按厘米取整，将是一个不同于标准差为 1 cm 的独立高斯噪声的误差模型。

首要局限是测量区间过短。增加帧数可以区分曲率变化并降低速度--加速度耦合；经过标定的相机误差协方差与同步时间戳会使置信估计更有依据。共同相机偏置、初始位置不确定、深度重建以及相关坐标误差，都不在独立误差情景的覆盖范围内。对于一般已知协方差 $\Sigma$，只需把 $J^TJ/\sigma^2$ 替换为 $J^T\Sigma^{-1}J$。

常数有效马格努斯力结构无法分辨非定常尾流、突发阵风或最后一次观测之后的受力变化。对于任意的未来扰动，不存在数据驱动的有限界。低自旋飘摆与反弹需要额外经过标定的模型。对于在到达球门之前触地的轨迹，我们的方法会报告另一种结果，而不是悄悄地向地面以下外推。

独立实现与合成验证检查的是方程与算法。它们并未确立这一特定足球、场地与相机服从所选模型。这里没有独立测量的完整轨迹；实验确认仍然有待完成。边界与合成抽样范围是明确的场景选择，而不是普遍的足球统计。

\section{结论}
题给观测支持一个条件性的名义预测，约为 $0.401$ s 时 $(11,2.60,1.97)$ m。一个六参数的重力--阻力--马格努斯模型给出从位置到事件定位预测的可复现路径，而解析基线揭示了外推假设的影响。不确定性必须成为答案的一部分：在 1 cm 坐标误差情景下，边际 95\% 半宽约为 1 m，远大于数值积分误差。结构压力情形必须单独报告。因此，这一通用方法返回的是一个点、一个条件不确定性区域、若干诊断量以及一个结果类别，而不是一个缺乏依据的精确普适预测。

\begin{thebibliography}{9}
\bibitem{upc} University Physics Competition. 2026 contest, Problem B. \url{https://www.uphysicsc.com/2026contest.html}. Accessed 7 November 2026.
\bibitem{asai} T. Asai, K. Seo, O. Kobayashi and R. Sakashita. Fundamental aerodynamics of the soccer ball. \emph{Sports Engineering} 10, 101--109 (2007). \url{https://doi.org/10.1007/BF02844207}.
\bibitem{ifab} The IFAB. Laws of the Game 2026/27, Law 2: The Ball. \url{https://www.theifab.com/laws/latest/the-ball/}. Accessed 7 November 2026.
\bibitem{ode} MathWorks. \texttt{ode45}: Solve nonstiff differential equations. \url{https://www.mathworks.com/help/matlab/ref/ode45.html}. Implementation run with MATLAB R2024a.
\bibitem{lsq} MathWorks. \texttt{lsqnonlin}: Solve nonlinear least-squares problems. \url{https://www.mathworks.com/help/optim/ug/lsqnonlin.html}.
\bibitem{rules} University Physics Competition. Contest rules, including generative-AI disclosure. \url{https://www.uphysicsc.com/contestrules.html}. Accessed 7 November 2026.
\end{thebibliography}

\section*{人工智能使用报告}
OpenAI Codex 协助了选题、建模决策、来源核查、独立数值验证与稿件准备。Claude Code 经由用户现有的服务商路由到 DeepSeek \texttt{deepseek-v4-flash}，协助了 MATLAB 实现。人类队员必须核验并对所有采纳内容承担责任。

\textbf{披露草稿：尚不可提交。}官方规则要求提供 AI 工具的精确输入与完整输出~\cite{rules}。本段只是角色概述，不能替代那些记录。采纳的完整对话与执行记录必须在提交前附上并接受审查；队号也有待确认。没有自动提交任何参赛材料。
\end{document}

```
