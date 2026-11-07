"""将MATLAB正式结果同步到英文论文；不包含绘图。"""
from pathlib import Path
import json
import numpy as np
from scipy.stats import chi2
ROOT=Path(__file__).resolve().parents[2]  # 仿真/Matlab版/ → 项目根
d=json.loads((ROOT/'data/派生数据/罚球结果.json').read_text(encoding='utf-8'))
assert d['mode']=='full'
names=['Gravity only','Constant acceleration','Drag only','Drag and effective Magnus']
rows=[]
for name,b in zip(names,d['baselines']):
    rows.append(f"{name} & {b['rmse']*1000:.3f} & {b['tCross']:.4f} & {b['yGoal']:.3f} & {b['zGoal']:.3f} \\\\")
sc=[]
for b in d['modelScenarios']:
    name=b['name'].replace('_',r'\_')
    delta=100*np.linalg.norm(np.array([b['yGoal'],b['zGoal']])-np.array(d['prediction']['point']))
    sc.append(f"{name} & {b['yGoal']:.4f} & {b['zGoal']:.4f} & {delta:.2f} \\\\")
bs=d['bootstrap']['summary']; syn=d['synthetic']['summary']; v=d['verification']; unc=d['uncertainty']
q=np.array(bs['goalPointQuantiles'])
drag_chi2=6*d['baselines'][2]['rmse']**2/.01**2
text=r'''\subsection{Point forecast and competing baselines}
The four independent starts converge to the same interior solution, with parameters
\begin{equation}
\widehat{\boldsymbol\theta}=(30.41546,9.41672,6.90801,0.370687,-0.060326,-0.357128).
\end{equation}
The first three entries have units m s$^{-1}$, the fourth is dimensionless, and the last two have units s$^{-1}$. The predicted crossing is $(y_g,z_g)=(2.599441,1.967005)$ m at $0.4005035$ s. This numerical precision permits reproducibility; the answer should be reported as $(11,2.60,1.97)$ m at $0.401$ s.
\begin{table}[H]\centering
\caption{Observation fit and extrapolation. Tiny interpolation residuals do not demonstrate physical validity.}
\begin{tabular}{lrrrr}\toprule
Model & RMSE (mm) & $t_g$ (s) & $y_g$ (m) & $z_g$ (m)\\\midrule
'''+'\n'.join(rows)+r'''
\bottomrule\end{tabular}\end{table}
'''+f'''The drag-only objective is $\\chi^2={drag_chi2:.3f}$ for two nominal residual degrees of freedom, corresponding to an approximate goodness-of-fit tail probability ${chi2.sf(drag_chi2,2):.3f}$ under the 1 cm Gaussian assumption. Thus these observations do not by themselves establish the need for the two extra transverse-force parameters at that noise level. At 2 mm noise, the same drag-only discrepancy gives $\\chi^2={drag_chi2*25:.1f}$, making the interpretation strongly dependent on measurement quality. The more flexible model is retained as a conditional curvature-aware forecast, not declared uniquely correct by its exact fit.
'''+r'''
\begin{figure}[H]\centering
\includegraphics[width=.78\linewidth]{罚球轨迹与观测.png}
\caption{Nominal fitted flight from launch to the first goal-plane crossing. Only the two early positions are observations.}
\end{figure}
At the centre plane $x=11.11$ m, the finite-radius sensitivity gives $(y,z)=(2.61689,1.97822)$ m at $0.40493$ s. The additional travel shifts the centre by approximately 2.1 cm in the goal plane. The nominal shot clears the ideal opening under both conventions.

\subsection{Measurement uncertainty and bootstrap stability}
For independent coordinate errors of standard deviation 0.01 m, the local prediction covariance is
\begin{equation}
C_g\simeq\begin{pmatrix}0.245132&0.000117\\0.000117&0.243555\end{pmatrix}\ \mathrm{m^2}.
\end{equation}
The marginal 95\% half-widths are 0.9704 m and 0.9673 m. At noise levels of 2, 5 and 20 mm, the respective half-widths scale to approximately 0.194, 0.485 and 1.94 m in each coordinate. This scaling is a local result; it is not a claim that bounded nonlinear resampling scales exactly at every noise level.
'''+f'''The 95\\% joint ellipse has semiaxes {unc['joint95Axes'][0]:.3f} and {unc['joint95Axes'][1]:.3f} m. The unscaled observation Jacobian condition number is {unc['conditionNumber']:.1f}; it depends on parameter units and is reported as a numerical diagnostic, not an invariant identifiability score.

All {bs['nSamples']} bootstrap fits terminate successfully; {bs['nBoundary']} hit at least one parameter boundary. All reach the goal plane, but {bs['nOutside']} fail the finite-radius opening-clearance criterion. The latter are retained in the coordinate summaries and figure. The empirical 2.5--97.5 percentile ranges of the re-estimated plane coordinates are $[{q[0,0]:.3f},{q[2,0]:.3f}]$ m in $y$ and $[{q[0,1]:.3f},{q[2,1]:.3f}]$ m in $z$. These are estimator resampling quantiles, not a posterior distribution for the true shot. The {100*bs['goalFractionAllDraws']:.1f}\\% within-opening frequency is a stability diagnostic only. Boundary hits and the region extending beyond the frame are material features of the result.
'''+r'''
\begin{figure}[H]\centering
\includegraphics[width=.88\linewidth]{球门预测椭圆与自助散点.png}
\caption{Goal-plane forecasts from 600 bootstrap refits with 1 cm coordinate noise. All plane crossings are retained, including forecasts outside the opening. The ellipse is the local 95\% joint confidence approximation.}
\end{figure}

\subsection{Verification and synthetic forecast testing}
MATLAB R2024a performs the main integration, optimization and resampling, with random seed 20261107. The independent Python implementation uses SciPy 1.17.1, DOP853 integration and bounded nonlinear least squares. It agrees with the nominal MATLAB goal-plane coordinates to less than $2\times10^{-11}$ m. This very small difference establishes implementation consistency, not real-world prediction accuracy.
'''+f'''The RK4 gravity benchmark has maximum position error ${v['rk4GravityMaxError']:.2g}$ m and the one-dimensional quadratic-drag benchmark ${v['dragAnalyticMaxError']:.2g}$ m. At fixed fitted parameters, step sizes 0.01, 0.005 and 0.0025 s change the predicted crossing by less than $3\\times10^{{-10}}$ m. Independently changing the event integrator's maximum step to 0.02, 0.01 and 0.005 s changes the result by less than $4\\times10^{{-12}}$ m. Tolerances are RelTol=$10^{{-9}}$ and AbsTol=$10^{{-11}}$ for event integration. Refinement checks concern integration error, not uncertainty in the fitted physics. No convergence order is claimed from roundoff-scale errors. An independent integrated energy check balances mechanical-energy loss against drag dissipation to $2\\times10^{{-13}}$ J.

For the synthetic tests, each of 150 truth vectors is drawn independently and uniformly over $v_x\\in[24,34]$, $v_y\\in[3,10]$, $v_z\\in[4,9]$ m s$^{{-1}}$, $C_D\\in[0.15,0.40]$, $b_y\\in[-0.15,0.15]$, $b_z\\in[-0.60,0.10]$ s$^{{-1}}$. These ranges define a designed test set, not observed population statistics. High-accuracy ode45 trajectories generate the early observations, independently of the RK4 fitting routine. Noise of standard deviation 1 cm is then added and each case is refitted from two starts.

All 150 inversions converge, with {syn['nBoundary']} boundary hits. All true trajectories reach the plane; {syn['nTruthPlaneCrossings']-syn['coverageDenominator']} fitted trajectories instead hit the ground first. For the {syn['coverageDenominator']} paired plane crossings, two-coordinate RMSE is {syn['predictionRmse_planeSamples']:.3f} m. The local 95\\% ellipse covers 146 of those 147 true crossings ({100*syn['coverageRate']:.2f}\\%); counting the three non-crossing forecasts as noncoverage yields 146/150 ({100*syn['coverageAllTruthCrossings']:.2f}\\%). Neither number proves exact 95\\% calibration, especially with only 150 cases and active constraints. Opening/outside/ground classification agrees in {round(syn['classificationAccuracy']*150)}/150 cases. Synthetic data obey the reference force family, so these results do not test structural correctness against real ball flight.
'''+r'''
\subsection{Structural stress cases}
\begin{table}[H]\centering
\caption{One-factor changes, each refitted to the same observations. Distances are relative to the nominal crossing; these are scenarios, not confidence bounds.}
\begin{tabular}{lrrr}\toprule
Scenario & $y_g$ (m) & $z_g$ (m) & Shift (cm)\\\midrule
'''+'\n'.join(sc)+r'''
\bottomrule\end{tabular}\end{table}
Each stress case fits the early positions to numerical precision, yet future positions differ. The largest shift in this chosen set is about 11.5 cm for the 0.5 s decay case. Steady transverse wind has a small residual effect after refitting because the free initial velocity and force parameters absorb much of its early influence; this does not imply that an unexpected later gust is harmless. Unobserved longitudinal spin causes roughly 9 cm shifts. The stress range is illustrative and does not bound all physically possible future disturbances.
\begin{figure}[H]\centering
\includegraphics[width=.9\linewidth]{模型情景预测比较.png}
\caption{Displacement from the nominal goal-plane forecast after one-factor refitting. Units are centimetres; the selected scenarios carry no probability weights.}
\end{figure}
'''
p=ROOT/'latex/AI/B题英文论文.tex';s=p.read_text(encoding='utf-8')
a=s.index('% VERIFIED_RESULTS_BEGIN');b=s.index('% VERIFIED_RESULTS_END')
s=s[:a]+'% VERIFIED_RESULTS_BEGIN\n'+text+'\n'+s[b:]
p.write_text(s,encoding='utf-8')
print('English numerical section updated from full MATLAB results.')

# 主结果采用全程约束反演；两种语言从同一份正式数据生成本节和摘要。
f=d['trajectoryFamily']; default=next(q for q in f['summary'] if q['tolerance']==.005)
pair=[f['members'][i-1] for i in f['exampleIndices']]
values={'COUNT':str(f['nCandidates']),'ACCEPT':str(default['nAccepted']),
        'SEPARATION':f"{f['exampleSeparation']:.4f}",
        'YMIN':f"{default['yRange'][0]:.4f}",'YMAX':f"{default['yRange'][1]:.4f}",
        'ZMIN':f"{default['zRange'][0]:.4f}",'ZMAX':f"{default['zRange'][1]:.4f}",
        'SPEEDMIN':f"{default['speedRange'][0]:.3f}",'SPEEDMAX':f"{default['speedRange'][1]:.3f}"}
rows=[]
for q in f['summary']:
    rows.append(f"{1000*q['tolerance']:.0f} & {q['nAccepted']}/{f['nCandidates']} & [{q['yRange'][0]:.3f},{q['yRange'][1]:.3f}] & [{q['zRange'][0]:.3f},{q['zRange'][1]:.3f}] " + r'\\')
values['ROWS']='\n'.join(rows)
rows=[]
for label,a in zip(['A','B'],pair):
    rows.append(f"{label} & {a['windY']:.1f} & {a['density']:.3f} & {a['Cd']:.4f} & {a['initialSpeed']:.3f} & {1000*a['maxResidual']:.3f} & {a['endpoint'][1]:.4f} & {a['endpoint'][2]:.4f} " + r'\\')
values['PAIRROWS']='\n'.join(rows)
checks=json.loads((ROOT/'data/派生数据/独立复核结果.json').read_text(encoding='utf-8'))
values['PYDIFF']=f"{checks['trajectory_family']['max_goal_difference_m']:.2g}"

def expand(s):
    for key,value in values.items():s=s.replace('@@'+key+'@@',value)
    return s

en=r'''\section{Full-flight constrained trajectory family: primary method}
\subsection{Coupled parameters and observation constraints}
The same force law acts from release, $t=0$, until the first terminal event. There is no ideal initial segment and no activation of wind or drag at $0.10$ s. For this reduced extension we retain steady lateral wind $\mathbf w=(0,w_y,0)$, constant density, constant drag coefficient and constant effective transverse coefficients, with $b_x=0$. Gusts, spin decay, drag crisis and camera bias are excluded from this family; they require additional information. The later structural stress tests remain separate comparisons.

Let $\eta=(w_y,\rho,C_D)$ and $q=(v_{0x},v_{0y},v_{0z},b_y,b_z)$. The ideal observation equations are
\begin{equation}
 \mathbf r(t_i;q,\eta)=\mathbf r_i,\quad i=1,2.
\end{equation}
We relax them to componentwise constraints,
\begin{equation}
 \mathcal F_\epsilon=\{(q,\eta)\in\mathcal B:\max_{i,j}|r_j(t_i;q,\eta)-r_{ij}|\le\epsilon\},
 \qquad \epsilon=0.005\ \mathrm m.
\end{equation}
The initial centre remains exactly the problem's $(0,0,0.11)$ m. The 5 mm tolerance is a chosen bounded-error scenario, not an inferred camera accuracy or a Gaussian standard deviation. It would match half a centimetre rounding interval only if rounding were the sole error. Results are also screened at 2, 10 and 20 mm. Exact numerical interpolation can occur; it does not establish exact physical knowledge.

Changing $\eta$ requires solving for $q$ again. A stronger aerodynamic deceleration can be compensated by a different launch velocity and transverse coefficient while the observations remain acceptable. The ball's launch kinetic energy, $E_0=m|\mathbf v_0|^2/2$, is a useful comparison of kicks. It is not the foot's force: that would additionally require contact duration and ball deformation. The full problem has eight scalar unknowns and six measured coordinates. Moreover, only the product $\rho C_D$ appears in the drag term, so these two quantities cannot be separately identified in this reduced model. Effective $b_y,b_z$ already absorb lift efficiency and spin; we do not also multiply them by an independently varied density factor.

\subsection{Reproducible sampling of compatible solutions}
We use seven equally spaced $w_y$ values on $[-8,8]$ m s$^{-1}$, three densities $\{1.1,1.225,1.35\}$ kg m$^{-3}$, and eleven drag values: ten equally spaced points on $[0.15,0.60]$ plus the calm-air fitted value. These are illustrative search settings, not measured wind or ball distributions. For each of the @@COUNT@@ profiles, bounded nonlinear least squares estimates $q$ from two starts:
\begin{equation}
 q^*(\eta)=\arg\min_{q\in\mathcal B_q}\sum_{i=1}^2
 \|\mathbf r(t_i;q,\eta)-\mathbf r_i\|^2.
\end{equation}
The effective transverse-coefficient bounds follow the reference model; launch bounds are tightened by the admissibility checks below. We evaluate the residual during fitting with RK4, step $0.0025$ s, and independently re-evaluate acceptance with \texttt{ode45} from $t=0$ (RelTol $10^{-11}$, AbsTol $10^{-12}$). Whole-flight propagation uses RelTol $10^{-9}$, AbsTol $10^{-11}$ and maximum step $0.002$ s. A release-contact guard avoids reporting the initial ground contact as a landing; the descending ground event and forward goal event are both localized.

We retain a profile representative only when every coordinate residual is within tolerance. Solver failures and rejected representatives are recorded. A failed least-squares representative is not proof that its entire profile has no feasible point: an $L_2$ minimizer need not minimize the largest residual. Likewise, two starts do not establish global optimality. Our displayed family is a sampled \emph{inner subset} of $\mathcal F_\epsilon$, with one representative per sampled profile; it is neither a full feasible envelope nor a confidence region. We do not attach probability weights, and a fraction of candidates crossing the opening is not a scoring probability.

\subsection{Computed family and coupled-kick example}
All @@COUNT@@ optimizations converged without parameter-bound hits. At 5 mm, @@ACCEPT@@ representatives pass the six-coordinate test. They all reach the goal plane before the ground and clear the finite-radius opening in this chosen set. Their launch speeds span @@SPEEDMIN@@--@@SPEEDMAX@@ m s$^{-1}$. Plane-coordinate ranges are $y\in[@@YMIN@@,@@YMAX@@]$ m and $z\in[@@ZMIN@@,@@ZMAX@@]$ m. The separate extrema need not belong to the same shot, and these sampled ranges do not bound untested solutions.
\begin{table}[H]\centering
\caption{Tolerance screening of the same profile representatives. Ranges are sampled spreads, not uncertainty bounds.}
\begin{tabular}{rrrr}\toprule
$\epsilon$ (mm) & Accepted & $y_g$ (m) & $z_g$ (m)\\\midrule
@@ROWS@@
\bottomrule\end{tabular}\end{table}
\begin{table}[H]\centering\small
\caption{Two compatible full-flight solutions, selected for largest pairwise separation within the 5 mm sample. Wind is present from release in both.}
\begin{tabular}{lrrrrrrr}\toprule
 & $w_y$ & $\rho$ & $C_D$ & $|v_0|$ & Error (mm) & $y_g$ & $z_g$\\\midrule
@@PAIRROWS@@
\bottomrule\end{tabular}
\par\smallskip\footnotesize Units: wind and speed in m s$^{-1}$; density in kg m$^{-3}$; plane coordinates in m.
\end{table}
The two plane crossings differ by @@SEPARATION@@ m despite both meeting the same observation tolerance. Their residuals, velocity vectors, effective force coefficients and complete trajectories are retained with the numerical output. This demonstrates non-uniqueness within the stated physical family without inventing an unobserved force switch. The observed spread is modest because these are steady, reduced-force profiles refitted to the observations; complexity alone would not justify a larger or more accurate result.
\begin{figure}[H]\centering
\includegraphics[width=.87\linewidth]{全程约束轨迹簇.png}
\caption{Accepted full-flight trajectories at 5 mm, with two contrasting representatives. Gray curves show selected physical possibilities, not probability density.}
\end{figure}
\begin{figure}[H]\centering
\includegraphics[width=.76\linewidth]{全程约束球门散点.png}
\caption{Zoomed goal-plane scatter for compatible profile representatives; color indicates prescribed $C_D$.}
\end{figure}

\subsection{Numerical evidence and interpretation of ``most likely''}
Python independently integrated all @@COUNT@@ candidates using DOP853 and reproduced the acceptance decisions without disagreements; maximum plane-point difference from MATLAB was $@@PYDIFF@@$ m. Independent refits of the two displayed representatives also agreed. Tests additionally verify wind at release, compensating initial-velocity changes, invariance under fixed $\rho C_D$, three maximum-step settings, drag dissipation relative to air, transverse-force orthogonality, ground-event localization and the ground-frame work balance
\begin{equation}
 \frac{d}{dt}\left(\frac12m|\mathbf v|^2+mgz\right)
 =(\mathbf F_D+\mathbf F_M)\cdot\mathbf v.
\end{equation}
These are numerical checks, not experimental confirmation. The earlier conditional baseline and synthetic tests below serve as additional diagnostics; their Gaussian confidence estimates must not be substituted for the bounded-error family's sampled spread.

Without a measurement likelihood and a probability model for the unknown environment, these data do not specify an unconditional most-likely crossing. To give the requested point estimate, we retain the calm-air maximum-likelihood reference $(11,2.60,1.97)$ m at $0.401$ s, explicitly conditional on its assumptions, alongside the compatible family. A predictive mode marginalized over nuisance parameters is not identified by these two observations alone. The interactive MATLAB view implements the same inverse model: releasing a parameter slider recomputes the entire shot and restarts its animation, while reporting fitted launch conditions and residual acceptance. It does not splice distinct shots into one trajectory.
'''
zh=r'''\section{全程约束轨迹簇：主要方法}
\subsection{耦合参数与观测约束}
同一套受力方程从出脚 $t=0$ 一直作用到第一个终止事件。不存在理想的初始段，也不在 $0.10$ s 时才开启风或阻力。本次降阶扩展保留恒定横风 $\mathbf w=(0,w_y,0)$、恒定空气密度、恒定阻力系数和恒定有效侧向系数，固定 $b_x=0$。阵风、自旋衰减、阻力危机和相机偏差不加入本轨迹簇；这些机制需要更多信息。后面的单因素压力测试仅作独立对照。

设 $\eta=(w_y,\rho,C_D)$，$q=(v_{0x},v_{0y},v_{0z},b_y,b_z)$。理想观测等式为
\begin{equation}
 \mathbf r(t_i;q,\eta)=\mathbf r_i,\quad i=1,2.
\end{equation}
将其放宽为逐坐标容差约束：
\begin{equation}
 \mathcal F_\epsilon=\{(q,\eta)\in\mathcal B:\max_{i,j}|r_j(t_i;q,\eta)-r_{ij}|\le\epsilon\},
 \qquad \epsilon=0.005\ \mathrm m.
\end{equation}
初始球心仍严格采用题给 $(0,0,0.11)$ m。5 mm 是人为声明的有界误差情景，不是由数据推断的相机精度，也不是高斯标准差。只有在唯一误差来自四舍五入时，它才对应保留到厘米的半个舍入区间。另比较 2、10、20 mm 容差。计算中可以出现近乎精确的插值，但不能把这种数值吻合等同于精确物理知识。

每次改变 $\eta$ 都必须重新求 $q$。更强的空气动力减速可由不同的初速度和侧向系数补偿，同时保持观测可接受。出脚动能 $E_0=m|\mathbf v_0|^2/2$ 可用于比较不同射门，但不是脚的作用力；后者还需要接触时间和球体变形信息。完整问题有八个标量未知量和六个位置约束。此外，阻力项中只出现乘积 $\rho C_D$，因此本模型不能把两者分别识别。有效 $b_y,b_z$ 已吸收升力效率和自旋影响，不再额外乘以独立变化的密度因子。

\subsection{可复现的相容解采样}
横风在 $[-8,8]$ m/s 上取七个等间隔值，密度取 $\{1.1,1.225,1.35\}$ kg/m$^3$，阻力系数在 $[0.15,0.60]$ 上取十个等间隔点并加入无风拟合值，共十一值。这些是示例搜索设置，不是实测风速或球型分布。对 @@COUNT@@ 组指定参数分别从两个初值进行有界非线性最小二乘：
\begin{equation}
 q^*(\eta)=\arg\min_{q\in\mathcal B_q}\sum_{i=1}^2
 \|\mathbf r(t_i;q,\eta)-\mathbf r_i\|^2.
\end{equation}
有效侧向系数边界沿用参考模型；初速度边界按下节物理筛选收紧。优化内使用步长 $0.0025$ s 的 RK4；接受前用从 $t=0$ 积分的 \texttt{ode45} 独立回代观测，RelTol=$10^{-11}$、AbsTol=$10^{-12}$。全程预测采用 RelTol=$10^{-9}$、AbsTol=$10^{-11}$、最大步长 $0.002$ s。出脚接触保护避免误把初始地面接触作为落地，下降触地和向前越过球门均进行事件定位。

仅当六个坐标残差全部满足容差时才保留代表解，并记录求解失败和未通过者。最小二乘代表未通过，不能证明该参数剖面不存在其它可行解，因为最小化平方和不等同于最小化最大残差；两个初值也不构成全局最优证明。所显示的族只是 $\mathcal F_\epsilon$ 的一个采样\emph{内部子集}，每个参数剖面仅保留一个代表，既不是全部可行解的包络，也不是置信区间。不为成员赋予概率，穿过门框的成员比例也不是进球概率。

\subsection{轨迹簇结果与耦合出脚示例}
@@COUNT@@ 组优化全部收敛，无参数边界命中。5 mm 容差下 @@ACCEPT@@ 个代表通过六坐标检查，在本组样本内均先到球门平面且满足有限半径净空条件。初速度大小范围为 @@SPEEDMIN@@--@@SPEEDMAX@@ m/s，平面坐标范围为 $y\in[@@YMIN@@,@@YMAX@@]$ m、$z\in[@@ZMIN@@,@@ZMAX@@]$ m。不同坐标的极值不一定属于同一次射门；这些采样范围不能限制未搜索的解。
\begin{table}[H]\centering
\caption{同一批代表解在不同容差下的筛选；坐标范围是采样离散程度，不是不确定性的上下界。}
\begin{tabular}{rrrr}\toprule
$\epsilon$ (mm) & 通过数量 & $y_g$ (m) & $z_g$ (m)\\\midrule
@@ROWS@@
\bottomrule\end{tabular}\end{table}
\begin{table}[H]\centering\small
\caption{5 mm 样本中平面距离最大的两个相容全程解；两者从出脚开始均有横风。}
\begin{tabular}{lrrrrrrr}\toprule
 & $w_y$ & $\rho$ & $C_D$ & $|v_0|$ & 误差(mm) & $y_g$ & $z_g$\\\midrule
@@PAIRROWS@@
\bottomrule\end{tabular}
\par\smallskip\footnotesize 风速和初速单位 m/s；密度 kg/m$^3$；平面坐标 m。
\end{table}
两个解均满足同一观测容差，球门平面穿越点却相距 @@SEPARATION@@ m。完整残差、初速度向量、有效侧向系数和轨迹随数值结果保存。这是在明确物理族内展示非唯一性，而非人为设置最后观测后才出现的受力变化。此处离散程度较小，是因为采用恒定降阶力并重新拟合观测；模型复杂并不自动产生更大或更准确的差异。
\begin{figure}[H]\centering
\includegraphics[width=.87\linewidth]{全程约束轨迹簇.png}
\caption{5 mm 容差内的全程轨迹及两个差异较大的代表。灰线表示所选物理可能性，不是概率密度。}
\end{figure}
\begin{figure}[H]\centering
\includegraphics[width=.76\linewidth]{全程约束球门散点.png}
\caption{相容代表解在球门平面的局部放大散点；颜色为指定阻力系数。}
\end{figure}

\subsection{数值证据与“最可能”的含义}
Python 使用 DOP853 独立积分全部 @@COUNT@@ 个候选，接受判定与 MATLAB 全部一致，球门点最大差异 $@@PYDIFF@@$ m；两个代表的独立重拟合亦一致。验证还覆盖出脚即有风、初速度联动、固定 $\rho C_D$ 的不变性、三档最大步长、相对空气的阻力耗散、侧向力正交性、触地定位以及地面系功与能平衡：
\begin{equation}
 \frac{d}{dt}\left(\frac12m|\mathbf v|^2+mgz\right)
 =(\mathbf F_D+\mathbf F_M)\cdot\mathbf v.
\end{equation}
这些是数值检验，尚非实验确认。后面的无风条件参考模型及合成检验作为补充诊断，其高斯置信估计不能替代有界误差族的采样范围。

若没有观测似然和未知环境的概率模型，仅凭题给信息不能定义无条件的最可能穿越点。为回应题目对点预测的要求，仍报告无风最大似然参考点 $(11,2.60,1.97)$ m、时刻 $0.401$ s，并明确其条件性，同时给出相容轨迹簇。对未知参数积分后的预测密度众数无法由这两次观测单独确定。MATLAB 交互界面采用同一反演模型：松开参数滑块后重新计算整次射门并从头播放，同时显示反求初速度及残差是否合格，不把不同射门拼成一条轨迹。
'''
abstract_en=expand(r'''Two position measurements do not uniquely identify a penalty kick when aerodynamic conditions are unknown. We apply gravity, quadratic drag relative to steady lateral wind, and an effective Magnus term from release throughout the flight. For prescribed wind, air density and drag coefficient, bounded nonlinear least squares jointly estimates initial velocity and two transverse-force coefficients. The six observation equations are relaxed to explicit componentwise tolerances. A grid of @@COUNT@@ parameter profiles yields @@ACCEPT@@ compatible representatives at a chosen 5 mm tolerance; their goal-plane coordinates span $y\in[@@YMIN@@,@@YMAX@@]$ m and $z\in[@@ZMIN@@,@@ZMAX@@]$ m, with two acceptable shots differing by @@SEPARATION@@ m. This family is a sampled inner subset of feasible solutions, not an exhaustive uncertainty envelope or a probability distribution. Density and drag are structurally confounded through their product. A calm-air conditional reference remains $(11,2.60,1.97)$ m at $0.401$ s; the observations alone do not define an unconditional most-likely point. We separately retain a 1 cm Gaussian-noise diagnostic, which gives approximately 0.97 m marginal 95\% half-widths and must not be confused with the bounded-error family's spread. Analytical limiting cases, integration refinement, independent Python calculations and synthetic tests establish numerical consistency. No independent experimental flight is available. The method returns a conditional point, compatible full-flight trajectories, observation residuals and terminal-event classifications, while making its physical and statistical assumptions explicit.''')
abstract_zh=expand(r'''在空气动力条件未知时，两次位置观测不能唯一确定点球轨迹。本模型从出脚开始全程考虑重力、相对恒定横风的二次阻力和有效马格努斯项。指定横风、空气密度和阻力系数后，通过有界非线性最小二乘联合反演初速度与两个有效侧向系数，并把六个观测等式放宽为逐坐标容差约束。@@COUNT@@ 个参数剖面中，@@ACCEPT@@ 个代表满足人为选定的 5 mm 容差，球门坐标范围为 $y\in[@@YMIN@@,@@YMAX@@]$ m、$z\in[@@ZMIN@@,@@ZMAX@@]$ m，其中两个相容射门相距 @@SEPARATION@@ m。该族是可行解的采样内部子集，既不是全部不确定性包络，也不是概率分布。密度和阻力系数通过乘积混淆。无风条件参考仍为 $0.401$ s 时穿过 $(11,2.60,1.97)$ m；题给观测不足以定义无条件的最可能点。另保留 1 cm 高斯噪声诊断，其局部 95\% 边际半宽约 0.97 m，不能与有界误差族的离散程度混同。解析极限、积分细化、Python 独立计算和合成检验支持数值一致性，但尚无独立实验飞行确认。方法返回条件点预测、相容全程轨迹、观测残差和终止事件分类，并明确全部物理及统计假设。''')
conclusion_en=expand(r'''The primary answer is a coupled full-flight inverse problem: environmental parameters and launch conditions must be inferred together under the two observation constraints. At a stated 5 mm coordinate tolerance, @@ACCEPT@@ of @@COUNT@@ sampled profile representatives are compatible, and two differ by @@SEPARATION@@ m at the goal. These ranges demonstrate non-uniqueness but are not global bounds or probabilities. The calm-air point $(11,2.60,1.97)$ m at $0.401$ s remains a conditional reference, not a uniquely established truth or an unconditional predictive mode. Gaussian measurement-error diagnostics are reported separately. Further observations, environmental measurements and independent flight data are needed to narrow the feasible set and calibrate probabilistic forecasts.''')
conclusion_zh=expand(r'''主要答案是一个耦合的全程反问题：环境参数与出脚条件必须在两次观测约束下共同求解。在明确声明的 5 mm 逐坐标容差下，@@COUNT@@ 个参数剖面中 @@ACCEPT@@ 个代表相容，两个代表在球门处相距 @@SEPARATION@@ m。这证明了非唯一性，但采样范围不是全局上下界或概率。$0.401$ s 时的无风点 $(11,2.60,1.97)$ m 仍是条件参考，不是唯一真值或无条件预测众数。高斯测量误差诊断单独报告。要收缩可行集并校准概率预测，还需要更多观测帧、环境测量和独立飞行数据。''')

physical_en=r"""\subsection{Physical admissibility is distinct from observation fit}
A small residual is necessary but insufficient. We additionally require a forward, initially non-descending launch, $10\le|\mathbf v_0|\le40$ m s$^{-1}$, and compliance with the declared environment domain. Optimizer component bounds are tightened to $v_{0x}\in[10,40]$ and $v_{0z}\in[0,25]$ m s$^{-1}$; $v_{0y}$ and $b_y,b_z$ retain their previous bounds. IFAB specifies mass 0.410--0.450 kg and circumference 0.68--0.70 m~\cite{ifab}; the fixed 0.43 kg, 0.11 m radius satisfy these requirements. A study of nine experienced players observed ball speeds of 19.2--32.6 m s$^{-1}$~\cite{nunome}. The adopted 40 m s$^{-1}$ cap includes a margin and is an explicit plausibility choice, not a universal human speed limit. No unsupported numerical upper bound on launch angle is imposed.

To prevent unrealistically large effective side forces, we screen
\begin{equation}
 C_{\perp}(t)=\frac{\|\mathbf b\times\mathbf u(t)\|}{K\|\mathbf u(t)\|^2}\le0.50.
\end{equation}
This converts the fitted coefficient to a dimensionless transverse-force measure rather than pretending to recover true spin. Football wind-tunnel experiments provide the force-coefficient definition and show dependence on rotation~\cite{asai}. The 0.50 cap and broad $C_D\in[0.15,0.60]$ domain are conservative modelling guards, not experimentally established limits for every football. Density and wind ranges remain unmeasured scenario bounds. The force check samples the dense solution at spacings no larger than 0.001 s; it is a numerical admissibility check, not a proof of a continuous-time bound.

Candidates must also remain above the ground before their terminal event, survive beyond the last observation and produce a localized terminal event within 2 s. We do \emph{not} constrain the inferred crossing to lie inside the goal: misses are legitimate outcomes. Observational compatibility, physical-screen acceptance and solver status are stored separately. A guard violation can indicate missing constraints, weak identifiability, a restricted search or a mismatch of assumptions; it does not by itself falsify Newton's equations.

At the default tolerance the maximum transverse coefficient among accepted representatives is @@CLMAX@@. All @@ACCEPT@@ pass the 0.50 guard. Tightening only this guard to 0.40 retains @@CL40@@, and to 0.30 retains @@CL30@@; relaxing to 0.60 retains @@ACCEPT@@. Speed caps 35, 40 and 45 m s$^{-1}$ give the same counts on this fixed sample. This is screening sensitivity, not a reoptimized global feasible envelope. Rejected candidates and reasons remain in the output.

"""
physical_zh=r"""\subsection{物理可接受性与观测拟合是两道检查}
残差小是必要条件，但还不够。新增筛选要求向前、初始不向地下出射，$10\le|\mathbf v_0|\le40$ m/s，并满足明确声明的环境搜索域。优化器分量边界收紧为 $v_{0x}\in[10,40]$、$v_{0z}\in[0,25]$ m/s；$v_{0y}$ 与 $b_y,b_z$ 沿用原边界。IFAB 规定质量 0.410--0.450 kg、周长 0.68--0.70 m~\cite{ifab}，固定质量 0.43 kg 和半径 0.11 m 满足要求。一项九名有经验球员的实验测得球速范围 19.2--32.6 m/s~\cite{nunome}。采用 40 m/s 上限留出余量，这是明确的合理性选择，不是人体球速的普适极限。不加入缺乏依据的数值出射角上限。

为排除过大的有效侧向力，检验无量纲系数：
\begin{equation}
 C_{\perp}(t)=\frac{\|\mathbf b\times\mathbf u(t)\|}{K\|\mathbf u(t)\|^2}\le0.50.
\end{equation}
这将拟合参数转化为侧向力系数，而非假装恢复真实自旋。足球风洞实验给出了力系数定义及其随旋转变化的证据~\cite{asai}。0.50 上限和较宽的 $C_D\in[0.15,0.60]$ 都是保守建模筛选，不是所有足球经实验证明的绝对边界。密度和风速仍为未实测的情景范围。对密集解以不大于 0.001 s 的间隔检查侧力系数；这是数值筛选，不是连续时间严格界的证明。

还要求轨迹在终止前不穿地下、经过最后观测时刻，并在 2 s 内定位终止事件。\emph{不}把落点在门框内作为反演约束，射偏也是合法物理结果。输出分别保存观测相容性、物理筛选和求解器状态。离谱候选可能源于欠约束、弱可辨识、搜索限制或假设不匹配，不能仅由此判定牛顿方程错误。

默认容差下，接受代表的最大侧力系数为 @@CLMAX@@，@@ACCEPT@@ 个均通过 0.50 阈值。仅收紧此阈值至 0.40 时保留 @@CL40@@ 个，至 0.30 时保留 @@CL30@@ 个；放宽至 0.60 仍为 @@ACCEPT@@ 个。速度上限取 35、40、45 m/s 对此固定样本的数量无影响。这是筛选敏感性，未重新优化为全局可行包络；不合格候选及其原因仍保存。

"""
values['CLMAX']=f"{max(a['maxTransverseCoefficient'] for a in f['members'] if a['accepted']):.4f}"
for cap,key in [(.3,'CL30'),(.4,'CL40')]:
    values[key]=str(next(a['accepted'] for a in f['admissibilitySensitivity'] if a['speedCap']==40 and a['coefficientCap']==cap))
en=en.replace(r'\subsection{Computed family and coupled-kick example}',physical_en+r'\subsection{Computed family and coupled-kick example}')
zh=zh.replace(r'\subsection{轨迹簇结果与耦合出脚示例}',physical_zh+r'\subsection{轨迹簇结果与耦合出脚示例}')
for lang,section,abstract,conclusion in [('英文',expand(en),abstract_en,conclusion_en),('中文',expand(zh),abstract_zh,conclusion_zh)]:
    path=ROOT/f'latex/AI/B题{lang}论文.tex';s=path.read_text(encoding='utf-8')
    begin='% FULL_FLIGHT_FAMILY_BEGIN';end='% FULL_FLIGHT_FAMILY_END'
    block=begin+'\n'+section+'\n'+end+'\n'
    if begin in s:
        a=s.index(begin);b=s.index(end)+len(end);s=s[:a]+block+s[b:]
    else:
        anchor=r'\section{Inference, event detection and uncertainty}' if lang=='英文' else r'\section{反演、事件检测与不确定性}'
        if anchor not in s:
            anchor=next(line for line in s.splitlines() if line.startswith(r'\section{') and ('反演' in line or '推断' in line))
        s=s.replace(anchor,block+'\n'+anchor,1)
    a=s.index(r'\begin{abstract}')+len(r'\begin{abstract}');b=s.index(r'\end{abstract}')
    s=s[:a]+'\n'+abstract+'\n'+s[b:]
    marker=r'\section{Conclusions}' if lang=='英文' else r'\section{结论}'
    a=s.index(marker)+len(marker);b=s.index(r'\begin{thebibliography}',a)
    s=s[:a]+'\n'+conclusion+'\n\n'+s[b:]
    s=s.replace('A physics-constrained estimate and its uncertainty','Full-flight inverse constraints and compatible trajectories')
    s=s.replace('一种受物理约束的估计及其不确定性','全程物理约束反演与相容轨迹簇')
    # 明确原六参数及高斯诊断为条件参考，不把它写成新版唯一答案。
    s=s.replace(r'\section{Inference, event detection and uncertainty}',r'\section{Conditional reference: likelihood and measurement diagnostics}')
    s=s.replace(r'\section{Numerical results and verification}',r'\section{Additional reference-model results and verification}')
    if r'\bibitem{nunome}' not in s:
        s=s.replace(r'\end{thebibliography}',r'\bibitem{nunome} H. Nunome, M. Lake, A. Georgakis and L. K. Stergioulas. Impact phase kinematics of instep kicking in soccer. \emph{Journal of Sports Sciences} 24(1), 11--22 (2006). \url{https://doi.org/10.1080/02640410400021450}.'+'\n'+r'\end{thebibliography}')
    path.write_text(s,encoding='utf-8')
print('Full-flight family, abstracts and conclusions synchronized in both languages.')
