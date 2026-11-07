"""由 MATLAB 正式结果 JSON 生成中英文单文件 LaTeX 论文。

只读取新接口字段；候选数、通过数、网格规模、示例编号均从 JSON 动态取得。
"""
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[2]  # 仿真/Matlab版/ -> 项目根
d = json.loads((ROOT / 'data/派生数据/罚球结果.json').read_text(encoding='utf-8'))
assert d.get('mode') == 'full', '需要正式结果 JSON'

f = d['trajectoryFamily']
sel = d['selected']
sp = sel['params']
sf = sel.get('fit') or {}
c = d['constants']
obs = d['observations']
ver = d.get('verification') or {}


# ---------- 数值格式 ----------
def fmt(x, nd=3):
    return f"{x:.{nd}f}" if isinstance(x, (int, float)) else '--'


def sci(x, nd=2):
    if not isinstance(x, (int, float)) or isinstance(x, bool):
        return None
    s = f"{x:.{nd}e}"
    mant, exp = s.split('e')
    return f"{mant}\\times 10^{{{int(exp)}}}"


def mm(x, nd=3):
    return fmt(1000 * x, nd) if isinstance(x, (int, float)) else '--'


summary = min(f['summary'], key=lambda q: abs(q['tolerance'] - sp['tolerance']))
count = f['nCandidates']
n_wind = len(f['windGrid'])
n_dens = len(f['densityGrid'])
n_cd = len(f['CdGrid'])
wind_min, wind_max = f['windGrid'][0], f['windGrid'][-1]
cd_min, cd_max = f['CdGrid'][0], f['CdGrid'][-1]

# 容差筛选表
tol_rows = []
for q in f['summary']:
    tol_rows.append(
        f"{1000 * q['tolerance']:.0f} & {q['nAccepted']}/{count} & "
        f"{q.get('nPlane', '--')} & {q.get('nOpening', '--')} & "
        f"{q.get('nGround', '--')} & "
        f"$[{q['yRange'][0]:.4f},{q['yRange'][1]:.4f}]$ & "
        f"$[{q['zRange'][0]:.4f},{q['zRange'][1]:.4f}]$ " + r'\\'
    )
tol_rows = '\n'.join(tol_rows)

# 两个示例（MATLAB 1-based 编号）
pair = [f['members'][i - 1] for i in f['exampleIndices']]
labels = ['A', 'B']
pair_rows = []
for label, a in zip(labels, pair):
    pair_rows.append(
        f"{label} & {a['windY']:.1f} & {a['density']:.3f} & {a['Cd']:.4f} & "
        f"{a['initialSpeed']:.3f} & {mm(a['maxResidual'])} & "
        f"{a['endpoint'][1]:.4f} & {a['endpoint'][2]:.4f} & "
        f"{a['tEnd']:.4f} " + r'\\'
    )
pair_rows = '\n'.join(pair_rows)

# 当前示例
end = sel.get('endpoint') or [None, None, None]
sel_theta = ", ".join(fmt(x, 4) for i, x in enumerate(sel['theta']) if i != 3)
sel_open_en = 'yes' if sel.get('isOpening') else 'no'
sel_open_zh = '通过' if sel.get('isOpening') else '未通过'

# 检查段落：只有 JSON 给出数值才写数字，否则指向 MATLAB 输出
en_check, zh_check = [], []
mp = ver.get('maxStepDifference')
if isinstance(mp, (int, float)) and not isinstance(mp, bool):
    en_check.append(f"The largest change in the goal-plane endpoint among the three step sizes is ${sci(mp)}$ m.")
    zh_check.append(f"三档步长下球门平面落点的最大变化为 ${sci(mp)}$ m。")
else:
    en_check.append("The three-step comparison is reported in the MATLAB check output.")
    zh_check.append("三档步长的比较见 MATLAB 检查输出。")

ew = ver.get('energyWorkResidual')
if isinstance(ew, (int, float)) and not isinstance(ew, bool):
    en_check.append(f"The maximum magnitude of this work-energy residual is ${sci(ew)}$ J.")
    zh_check.append(f"该功与能累积残差的最大值为 ${sci(ew)}$ J。")
else:
    en_check.append("The work-energy balance is reported in the MATLAB check output.")
    zh_check.append("功与能平衡结果见 MATLAB 检查输出。")

rej = ver.get('implausibleSpeedAndForceRejected')
if isinstance(rej, bool):
    if rej:
        en_check.append("Checks with deliberately excessive launch speed and side force were rejected as intended.")
        zh_check.append("人为设置的过大初速度和侧向力均被检查程序拒绝。")
    else:
        en_check.append("The speed and side-force screening removed no parameter choice.")
        zh_check.append("速度与侧力筛选未剔除参数。")
elif isinstance(rej, (int, float)):
    en_check.append(f"The speed and side-force screening removed {int(rej)} parameter choices.")
    zh_check.append(f"速度与侧力筛选剔除了 {int(rej)} 组参数。")
else:
    en_check.append("The number of parameter choices removed by the speed and side-force screening is reported in the MATLAB check output.")
    zh_check.append("速度与侧力筛选剔除的组数见 MATLAB 检查输出。")
check_en = ' '.join(en_check)
check_zh = ''.join(zh_check)


values = {
    'COUNT': str(count),
    'ACCEPT': str(summary['nAccepted']),
    'TOLMM': f"{1000 * sp['tolerance']:.0f}",
    'NWIND': str(n_wind),
    'WMIN': fmt(wind_min, 3),
    'WMAX': fmt(wind_max, 3),
    'NDENS': str(n_dens),
    'DENSITIES': ", ".join(fmt(x, 3) for x in f['densityGrid']),
    'NCD': str(n_cd),
    'CDMIN': fmt(cd_min, 2),
    'CDMAX': fmt(cd_max, 2),
    'YMIN': fmt(summary['yRange'][0], 4),
    'YMAX': fmt(summary['yRange'][1], 4),
    'ZMIN': fmt(summary['zRange'][0], 4),
    'ZMAX': fmt(summary['zRange'][1], 4),
    'SMIN': fmt(summary['speedRange'][0], 3),
    'SMAX': fmt(summary['speedRange'][1], 3),
    'TOLROWS': tol_rows,
    'PAIRROWS': pair_rows,
    'PAIRSEP': fmt(f['exampleSeparation'], 4),
    'IDXA': str(f['exampleIndices'][0]),
    'IDXB': str(f['exampleIndices'][1]),
    'SPWY': fmt(sp['windY'], 1),
    'SPRHO': fmt(sp['density'], 3),
    'SPCD': fmt(sp['Cd'], 3),
    'SPMSTEP': fmt(sp['maxStep'], 3),
    'SELTHETA': sel_theta,
    'SELCD': fmt(sel['theta'][3], 4),
    'SELMAXRES': mm(sf.get('maxResidual')),
    'SELSPEED': fmt(sel.get('initialSpeed'), 3),
    'SELTEND': fmt(sel.get('tEnd'), 4),
    'SELY': fmt(end[1], 4),
    'SELZ': fmt(end[2], 4),
    'SELSTATUS': str(sel.get('status', '--')).replace('_', r'\_'),
    'SELOPEN_EN': sel_open_en,
    'SELOPEN_ZH': sel_open_zh,
    'CHECKEN': check_en,
    'CHECKZH': check_zh,
    'R': fmt(c['R'], 2),
    'M': fmt(c['m'], 2),
    'G': fmt(c['g'][2], 2),
    'OBS1': ", ".join(fmt(x, 2) for x in obs['positions'][0]),
    'OBS2': ", ".join(fmt(x, 2) for x in obs['positions'][1]),
    'T1': fmt(obs['times'][0], 2),
    'T2': fmt(obs['times'][1], 2),
}


def render(template):
    for key, value in values.items():
        template = template.replace('@@' + key + '@@', value)
    return template


# ---------- 英文正文模板 ----------
template_en = r'''\documentclass[11pt,a4paper]{article}
\usepackage[margin=25mm]{geometry}
\usepackage{amsmath,amssymb,booktabs,graphicx,fancyhdr,hyperref,float}
\usepackage[T1]{fontenc}
\hypersetup{colorlinks=true,linkcolor=blue,urlcolor=blue,citecolor=blue}
\pagestyle{fancy}\fancyhf{}\fancyhead[C]{TEAM NUMBER REQUIRED}\fancyfoot[C]{\thepage}
\setlength{\headheight}{15pt}\setlength{\parskip}{3pt}\setlength{\parindent}{0pt}
\usepackage{caption}
\captionsetup{font=small,labelfont=bf}
\graphicspath{{../../img/}}
\begin{document}
\begin{center}
{\large Team: TEAM NUMBER REQUIRED}\par\vspace{12mm}
{\Large Problem B}\par\vspace{10mm}
{\LARGE Predicting a Penalty Kick from 0.10 Seconds of Motion}\par\vspace{8mm}
\end{center}
\vspace{8mm}
\begin{abstract}
Two early positions of a penalty kick are given, but the wind, air density and drag coefficient are not measured. This paper asks what the rest of the flight can look like under one clear force model. We use gravity, quadratic drag relative to the air, and an effective sideways force written as a cross product, all acting from the kick to the goal plane. The sideways term is not a measurement of spin; it combines spin and lift. For each listed wind, density and drag value we fit the launch velocity and two sideways coefficients to the two positions. A fit is accepted only when all six coordinate residuals lie within a stated allowed error. The lists contain @@COUNT@@ parameter choices; @@ACCEPT@@ are accepted at @@TOLMM@@ mm. Their goal-plane coordinates span $y\in[@@YMIN@@,@@YMAX@@]$ m and $z\in[@@ZMIN@@,@@ZMAX@@]$ m, and two accepted trajectories differ by @@PAIRSEP@@ m. For the current example at wind @@SPWY@@ m s$^{-1}$, density @@SPRHO@@ kg m$^{-3}$ and $C_D=@@SPCD@@$, the centre crosses the plane at $(@@SELY@@,@@SELZ@@)$ m after @@SELTEND@@ s. The spread shows that the two observations do not select one landing point; it is not a probability. MATLAB checks use three step sizes, a work-energy balance and event tests. No independent measured full flight is available, so experimental confirmation remains to be done. Only the listed parameter values are searched, and the ranges in the figures do not guarantee that all possible cases are included.
\end{abstract}
\clearpage

\section{Problem and approach}
This section states what we compute and in what order to read the paper.

The ball centre is given at the kick, $\mathbf r_0=(0,0,@@R@@)$ m, and at $t_1=@@T1@@$ s and $t_2=@@T2@@$ s:
\begin{equation}
\mathbf r_1=(@@OBS1@@)\ \mathrm m,\qquad
\mathbf r_2=(@@OBS2@@)\ \mathrm m.
\end{equation}
We take $x$ toward the goal, $y$ across the goal and $z$ upward. These six numbers are the only position data. Wind, density, drag and launch conditions are unknown, so one trajectory is not determined. We ask instead: for a listed set of environmental values, which launch conditions reproduce both positions within an allowed error, and where do those trajectories cross $x=11$ m? The reading order is as follows: Section 2 writes the force law, Section 3 explains how the observations select parameter choices, Section 4 gives the numbers and figures, Section 5 lists checks and limits, and Section 6 answers the problem.

A trajectory is retained when it meets the allowed coordinate error and the physical checks. Figure~\ref{fig:flow} shows the calculation: choose an environment, solve for launch conditions, and check the resulting flight. Rejected cases are recorded alongside accepted ones.
\begin{figure}[H]\centering
\includegraphics[width=\linewidth]{系统与反演流程.png}
\caption{Coordinate convention and the path from two observations to accepted trajectories. The left panel shows the horizontal geometry; the right panel shows listed environment values, fitted launch conditions, the allowed-error and physical checks, and the recorded terminal events.}
\label{fig:flow}
\end{figure}

\section{Physical model from kick to goal}
@@PHYSICS@@

\section{Finding parameters from the observations}
This section explains how the two positions become parameter choices.

For each listed triple $\eta=(w_y,\rho,C_D)$ we fit $q=(v_{0x},v_{0y},v_{0z},b_y,b_z)$ by bounded nonlinear least squares,
\begin{equation}
q^*(\eta)=\arg\min_{q\in\mathcal B}\sum_{i=1}^{2}
\|\mathbf r(t_i;q,\eta)-\mathbf r_i\|^2,
\label{eq:fit}
\end{equation}
using \texttt{lsqnonlin}~\cite{lsq} from two starting points. The fit uses fourth-order Runge--Kutta integration with step $0.0025$ s. Acceptance is re-evaluated with \texttt{ode45} at RelTol=$10^{-11}$ and AbsTol=$10^{-12}$; whole-flight propagation uses RelTol=$10^{-9}$, AbsTol=$10^{-11}$ and maximum step $0.002$ s. A fit is checked on all six coordinates,
\begin{equation}
|r_j(t_i;q,\eta)-r_{ij}|\le\epsilon,
\qquad i=1,2,\quad j=x,y,z,
\label{eq:tol}
\end{equation}
at $\epsilon=2,@@TOLMM@@,10,20$ mm; @@TOLMM@@ mm is the default. The allowed error is a stated scenario, not a measured camera accuracy.

The search uses @@NWIND@@ wind values on $[@@WMIN@@,@@WMAX@@]$ m s$^{-1}$, @@NDENS@@ density values $\{@@DENSITIES@@\}$ kg m$^{-3}$, and @@NCD@@ drag values on $[@@CDMIN@@,@@CDMAX@@]$, giving @@COUNT@@ parameter choices. These values are example settings, not measured distributions. Only the listed parameter values are searched, and the ranges in the figures do not guarantee that all possible cases are included.

Before a fit is called accepted it also passes a physical screening. The launch must point forward and not downward, the initial speed must lie between 10 and 40 m s$^{-1}$, and along the sampled trajectory the side-force coefficient must satisfy $C_\perp\le0.50$. The ball must stay above the ground before the terminal event, and the event must be located within 2 s. An accepted fit must pass both the six-coordinate test and this screening; observation fit and physical screening are stored separately. The speed range and the 0.50 cap are conservative model screens, not universal limits for players or footballs. We do not force the crossing into the goal: a miss is a recorded outcome, and the goal is used only to describe where the ball goes.

Two further points explain why the result is a set rather than one shot. First, the fitting problem has more unknowns than the six coordinates. Changing the wind or the drag changes the air force, and a different launch velocity can make up for much of that change at the two observed times. The fit therefore finds different launch conditions for different listed environments, even though all of them look the same at $t=0.05$ and $t=0.10$ s. Second, the two observations constrain the early path only, and a small difference in the fitted sideways coefficient keeps acting after the last observation. The accepted set is the result of searching the listed values, not a proof that nothing outside the grid exists.

@@INFERENCE_DETAILS@@

\section{Results}
This section reports the current example, two accepted trajectories with different landings, and the screening figures.

For the current example the listed values are wind @@SPWY@@ m s$^{-1}$, density @@SPRHO@@ kg m$^{-3}$, $C_D=@@SPCD@@$, allowed error @@TOLMM@@ mm, and maximum step @@SPMSTEP@@ s. The fit gives
\begin{equation}
(v_{0x},v_{0y},v_{0z},b_y,b_z)=(@@SELTHETA@@),
\end{equation}
with the first three entries in m s$^{-1}$ and the last two in s$^{-1}$; the prescribed drag coefficient is $C_D=@@SELCD@@$. The largest coordinate residual is @@SELMAXRES@@ mm, and the initial speed is @@SELSPEED@@ m s$^{-1}$. The ball first crosses $x=11$ m at
\begin{equation}
(y_g,z_g)=(@@SELY@@,@@SELZ@@)\ \mathrm m,
\qquad t_g=@@SELTEND@@\ \mathrm s.
\end{equation}
The recorded terminal status is \texttt{@@SELSTATUS@@}, and the finite-radius opening check is @@SELOPEN_EN@@. These numbers belong to this one parameter choice; they are not a unique answer, and no probability is implied.

Table~\ref{tab:tol} repeats the fit for every listed parameter choice at four allowed errors. Reading the @@TOLMM@@ mm row: @@ACCEPT@@ of @@COUNT@@ choices are accepted, their crossings span the stated $y$ and $z$ ranges, and the Plane and Opening columns count arrivals and opening clearances. Ground-first outcomes are listed separately; they are not silently extrapolated below the ground.
\begin{table}[htbp]\centering
\caption{Allowed-error screening of the same parameter choices. Ranges are spreads over the sampled values, not uncertainty bounds.}
\label{tab:tol}
\begin{tabular}{rrrrrrr}\toprule
$\epsilon$ (mm) & Accepted & Plane & Opening & Ground & $y_g$ (m) & $z_g$ (m)\\\midrule
@@TOLROWS@@
\bottomrule\end{tabular}\end{table}

The two displayed examples are taken from the accepted set at @@TOLMM@@ mm; their positions in the MATLAB arrays are @@IDXA@@ and @@IDXB@@. Table~\ref{tab:pair} lists their values. Both reproduce the two positions within @@TOLMM@@ mm, yet their crossings differ by @@PAIRSEP@@ m.

The two rows show the compensation directly. They use different wind and drag values, so the aerodynamic force differs, and the fitted launch speeds are not identical. Each shot also takes its own sideways coefficient to keep the same early positions. These differences also affect the rest of the flight and produce a gap at the goal plane. The gap is modest, and that is expected: these are steady force values fitted to the same two positions, and a more complex force law would not automatically make the answer more certain.
\begin{table}[htbp]\centering\small
\caption{Two accepted trajectories with different goal-plane crossings. Wind acts from the kick in both.}
\label{tab:pair}
\begin{tabular}{lrrrrrrrr}\toprule
 & $w_y$ & $\rho$ & $C_D$ & $|v_0|$ & Error (mm) & $y_g$ & $z_g$ & $t_g$\\\midrule
@@PAIRROWS@@
\bottomrule\end{tabular}
\par\smallskip\footnotesize Units: wind and speed in m s$^{-1}$; density in kg m$^{-3}$; plane coordinates in m; time in s.
\end{table}

Figure~\ref{fig:separation} explains why different goal points can fit the same observations. In the left panel, all six coordinate errors are within the allowed band for both examples. In the right panel, their late paths separate. Their launch conditions compensate for different aerodynamic forces during the measured interval, but this compensation does not keep the later positions identical.
\begin{figure}[htbp]\centering
\includegraphics[width=.96\linewidth]{观测约束与轨迹分离.png}
\caption{The two accepted trajectories satisfy all six coordinate checks (left), but their late paths separate in the horizontal plane (right). The shaded band is the stated allowed error; displacements are shown at true scale.}
\label{fig:separation}
\end{figure}

Figure~\ref{fig:coupling} shows how the inputs and fitted quantities change together. A larger density--drag product generally requires a larger launch speed for these observations (left). The sideways coefficient changes with wind (middle). Increasing the allowed measurement error accepts more parameter choices (right).
\begin{figure}[htbp]\centering
\includegraphics[width=\linewidth]{参数耦合与筛选敏感性.png}
\caption{Parameter coupling and screening: density--drag product against fitted launch speed for three wind slices (left); wind against the fitted $b_z$ for accepted trajectories (middle); accepted count against the allowed error (right). All counts refer to the listed parameter choices.}
\label{fig:coupling}
\end{figure}

The goal-plane view in Figure~\ref{fig:scatter} shows the range of crossing points from the accepted trajectories. Colour identifies the specified drag coefficient; the number of plotted points in a region carries no probability weight.
\begin{figure}[htbp]\centering
\includegraphics[width=.72\linewidth]{全程约束球门散点.png}
\caption{Goal-plane crossings of the accepted trajectories. Colour marks the listed drag coefficient. The spread is a property of the sampled parameter choices, not a probability distribution.}
\label{fig:scatter}
\end{figure}

A separate interactive view allows wind in $[-15,15]$ m s$^{-1}$, density in $[0.90,1.40]$ kg m$^{-3}$ and $C_D\in[0.08,0.80]$; this wider exploration is distinct from the paper grid and uses the same screening.

\section{Numerical checks and scope}
This section lists the checks made inside MATLAB and the limits of the result.

All calculations use the same equations and event definitions. The three step sizes differ only in the maximum step of the integration, so they measure numerical sensitivity rather than a new model. Using the balance derived in Eq.~\eqref{eq:energy}, we check the accumulated residual
\begin{equation}
R_E(t)=E(t)-E(0)-\int_0^t(\mathbf F_D+\mathbf F_M)\cdot\mathbf v\,d\tau.
\end{equation}
Event checks confirm that ground contact, the forward crossing of $x=11$ m and any time-out are located and recorded separately. @@CHECKEN@@

These checks test the code and the internal balance. No independent measured full flight was available, so they do not confirm that this particular ball, field and camera obey the model. The reported crossings and ranges apply only to the listed parameter choices and the stated allowed error. They carry no probability and no guarantee outside the searched values. A wider allowed error accepts more trajectories, but it does not turn the sampled spread into a set that must contain the true landing point.

\section{Conclusion}
This section answers the problem.

From the two early positions, a single trajectory cannot be recovered when the wind and aerodynamic coefficients are unknown. For each listed environment we can still fit launch conditions that reproduce both positions, and those accepted trajectories cross the goal plane at different points. For the current example at wind @@SPWY@@ m s$^{-1}$ and density @@SPRHO@@ kg m$^{-3}$, the centre crosses at $(@@SELY@@,@@SELZ@@)$ m after @@SELTEND@@ s. Across the listed choices, @@ACCEPT@@ of @@COUNT@@ fits pass the @@TOLMM@@ mm test, with crossings spanning $y\in[@@YMIN@@,@@YMAX@@]$ m and $z\in[@@ZMIN@@,@@ZMAX@@]$ m; two accepted trajectories differ by @@PAIRSEP@@ m. These ranges describe the sampled grid only and are not probabilities. More observation times, measured wind and a measured drag coefficient would narrow the answer.

\begin{thebibliography}{9}
\bibitem{upc} University Physics Competition. 2026 contest, Problem B. \url{https://www.uphysicsc.com/2026contest.html}. Accessed 7 November 2026.
\bibitem{ifab} The IFAB. Laws of the Game, Law 2: The Ball. \url{https://www.theifab.com/laws/latest/the-ball/}. Accessed 7 November 2026.
\bibitem{asai} T. Asai, K. Seo, O. Kobayashi and R. Sakashita. Fundamental aerodynamics of the soccer ball. \emph{Sports Engineering} 10, 101--109 (2007). \url{https://doi.org/10.1007/BF02844207}.
\bibitem{nasa} NASA Glenn Research Center. Drag Equation. \url{https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/drag-equation/}.
\bibitem{ode} MathWorks. \texttt{ode45}: Solve nonstiff differential equations. \url{https://www.mathworks.com/help/matlab/ref/ode45.html}. Implementation run with MATLAB R2024a.
\bibitem{lsq} MathWorks. \texttt{lsqnonlin}: Solve nonlinear least-squares problems. \url{https://www.mathworks.com/help/optim/ug/lsqnonlin.html}.
\end{thebibliography}

\smallskip\small
\section*{Report on the use of AI}
AI tools assisted with code generation, numerical checks and manuscript drafting. Human team members must verify all adopted content and take responsibility for it.

\textbf{Team number and the complete AI records are still to be added; this is not a final submission.} The contest rules require the exact input and complete output of the AI tools used; that record must be appended and reviewed before submission.
\end{document}
'''

# ---------- 中文正文模板 ----------
template_zh = r'''\documentclass[11pt,a4paper,UTF8,fontset=windows]{ctexart}
\usepackage[margin=25mm]{geometry}
\usepackage{amsmath,amssymb,booktabs,graphicx,fancyhdr,hyperref,float}
\hypersetup{colorlinks=true,linkcolor=blue,urlcolor=blue,citecolor=blue}
\pagestyle{fancy}\fancyhf{}\fancyhead[C]{队伍：TEAM NUMBER REQUIRED}\fancyfoot[C]{\thepage}
\setlength{\headheight}{15pt}\setlength{\parskip}{3pt}\setlength{\parindent}{0pt}
\usepackage{caption}
\captionsetup{font=small,labelfont=bf}
\graphicspath{{../../img/中文/}}
\begin{document}
\begin{center}
{\large 队伍：TEAM NUMBER REQUIRED}\par\vspace{12mm}
{\Large B 题}\par\vspace{10mm}
{\LARGE 由 0.10 秒运动预测点球}\par\vspace{8mm}
\end{center}
\vspace{6mm}
\begin{center}
{\small 中文译本仅供队友阅读，队号待定，正式提交以英文稿为准。}
\end{center}
\vspace{4mm}
\begin{abstract}
题给点球早期的两个位置，但风、空气密度和阻力系数都没有测量。本文只用一个明确的受力模型，回答球的其余飞行可能是什么样。模型从出脚到球门平面全程考虑重力、相对空气的二次阻力和写成叉乘的有效侧向力；侧向项不是自旋的测量，而是自旋与升力的合并效果。对每一组列出的风速、密度和阻力系数，我们用两个位置反求初速度与两个侧向系数。只有六个坐标残差都在允许误差内，这组参数才被接受。列出的一共有 @@COUNT@@ 组参数，其中 @@ACCEPT@@ 组在 @@TOLMM@@ mm 允许误差下被接受；它们的球门平面坐标覆盖 $y\in[@@YMIN@@,@@YMAX@@]$ m、$z\in[@@ZMIN@@,@@ZMAX@@]$ m，两条被接受的轨迹落点相距 @@PAIRSEP@@ m。当前示例取风速 @@SPWY@@ m/s、密度 @@SPRHO@@ kg/m$^3$、$C_D=@@SPCD@@$，球心在 @@SELTEND@@ s 时穿过球门平面于 $(@@SELY@@,@@SELZ@@)$ m。落点的分散说明两次观测不能选出唯一落点，但它不是概率。MATLAB 检查包括三档步长、功与能平衡和事件检查。没有独立的实测全程飞行，因此这些检查验证计算过程，模型还需实测飞行确认。只搜索列出的参数，图中的范围不能保证包含所有可能情况。
\end{abstract}
\clearpage

\section{问题与求解思路}
本节说明我们计算什么，以及全文按什么顺序阅读。

球心在出脚时的位置为 $\mathbf r_0=(0,0,@@R@@)$ m，另外两个时刻为 $t_1=@@T1@@$ s 与 $t_2=@@T2@@$ s：
\begin{equation}
\mathbf r_1=(@@OBS1@@)\ \mathrm m,\qquad
\mathbf r_2=(@@OBS2@@)\ \mathrm m.
\end{equation}
取 $x$ 指向球门、$y$ 平行球门线、$z$ 向上。这六个数就是全部位置数据。风、密度、阻力和出脚条件都未知，因而一条轨迹并不能被唯一确定。我们改问：对一组列出的环境取值，哪些出脚条件能在允许误差内重现两个位置，这些轨迹又在哪里穿过 $x=11$ m？阅读顺序如下：第 2 节写出受力方程，第 3 节说明观测如何筛出参数取值，第 4 节给出数字和插图，第 5 节列出计算检查与适用范围，第 6 节回答问题。

本文保留同时满足允许误差和物理检查的轨迹。图~\ref{fig:flow}展示计算顺序：先设定环境，再求出脚条件，最后检查整段飞行。未通过检查的结果也会记录。
\begin{figure}[H]\centering
\includegraphics[width=\linewidth]{系统与反演流程.png}
\caption{坐标约定与从两次观测到符合观测轨迹的流程。左图为水平几何示意；右图给出列出的环境取值、反求出的出脚条件、允许误差与物理检查，以及记录下来的终止事件。}
\label{fig:flow}
\end{figure}

\section{从出脚到球门的物理模型}
@@PHYSICS@@

\section{用观测求参数}
本节说明两次位置如何变成一组组参数取值。

对每一组列出的 $\eta=(w_y,\rho,C_D)$，用有界非线性最小二乘反求 $q=(v_{0x},v_{0y},v_{0z},b_y,b_z)$：
\begin{equation}
q^*(\eta)=\arg\min_{q\in\mathcal B}\sum_{i=1}^{2}
\|\mathbf r(t_i;q,\eta)-\mathbf r_i\|^2,
\label{eq:fit}
\end{equation}
使用 \texttt{lsqnonlin}~\cite{lsq}，并从两个初值出发。拟合时用四阶Runge--Kutta方法，步长为$0.0025$ s；接受前另用\texttt{ode45}回代观测，RelTol=$10^{-11}$、AbsTol=$10^{-12}$。全程预测用RelTol=$10^{-9}$、AbsTol=$10^{-11}$、最大步长$0.002$ s。随后逐个检查六个坐标：
\begin{equation}
|r_j(t_i;q,\eta)-r_{ij}|\le\epsilon,
\qquad i=1,2,\quad j=x,y,z,
\label{eq:tol}
\end{equation}
分别在 $\epsilon=2,@@TOLMM@@,10,20$ mm 下筛选，默认取 @@TOLMM@@ mm。允许误差是明确声明的情景，不是测量得到的相机精度。

搜索使用 $[@@WMIN@@,@@WMAX@@]$ m/s 上的 @@NWIND@@ 个风速、$\{@@DENSITIES@@\}$ kg/m$^3$ 的 @@NDENS@@ 个密度值、以及 $[@@CDMIN@@,@@CDMAX@@]$ 上的 @@NCD@@ 个阻力系数，共 @@COUNT@@ 组参数。这些是示例设置，不是实测分布。

一组参数还要通过物理筛选才算被接受：出脚方向向前且不向下，初速度在 10 到 40 m/s 之间，沿采样轨迹的侧力系数满足 $C_\perp\le0.50$；球在终止事件前必须位于地面以上，终止事件要在 2 s 内定位到。观测拟合与物理筛选分别记录，只有两项都通过才被接受。速度范围和 0.50 上限是保守的模型筛选，不是运动员或足球的普适极限。我们不强迫球进门：射偏也是记录下来的结果，球门只用来描述球飞到哪里。

结果之所以是一组轨迹而不是一次射门，有两个原因。第一，未知量比六个坐标多。改变风速或阻力就改变了空气力，不同的初速度又能在两个观测时刻上补偿其中的大部分差别，于是不同的环境取值会反求出不同的出脚条件，而在 $t=0.05$ 与 $t=0.10$ s 处看起来一样。第二，两次观测只约束前段，拟合出的侧向系数只要有一点不同，在最后观测之后仍继续作用。被接受的集合是搜索列出的参数得到的结果，并不证明网格之外没有别的解。

@@INFERENCE_DETAILS@@

\section{计算结果}
本节先给出当前示例的落点与时间，再给出两条符合观测但落点不同的轨迹，最后给出筛选图表。

当前示例取风速 @@SPWY@@ m/s、密度 @@SPRHO@@ kg/m$^3$、$C_D=@@SPCD@@$、允许误差 @@TOLMM@@ mm、最大步长 @@SPMSTEP@@ s。反求结果为
\begin{equation}
(v_{0x},v_{0y},v_{0z},b_y,b_z)=(@@SELTHETA@@),
\end{equation}
前三个量的单位是 m/s，后两个是 s$^{-1}$；阻力系数取指定值 $C_D=@@SELCD@@$。最大坐标残差为 @@SELMAXRES@@ mm，初速度为 @@SELSPEED@@ m/s。球心第一次穿过 $x=11$ m 的位置和时刻为
\begin{equation}
(y_g,z_g)=(@@SELY@@,@@SELZ@@)\ \mathrm m,
\qquad t_g=@@SELTEND@@\ \mathrm s.
\end{equation}
记录的终止状态为 \texttt{@@SELSTATUS@@}，有限半径的球门口检查为 @@SELOPEN_ZH@@。这些数字只属于这一组参数，不是唯一答案，也不代表概率。

表~\ref{tab:tol}把同一套拟合用于所有列出的参数，并分别在四种允许误差下筛选。读 @@TOLMM@@ mm 这一行：@@COUNT@@ 组中有 @@ACCEPT@@ 组被接受，落点覆盖表中 $y$、$z$ 范围；“到平面”和“过门口”两列分别记录到达球门平面及通过球门口检查的组数。先触地的结果单独列出，不会把轨迹悄悄外推到地面以下。
\begin{table}[htbp]\centering
\caption{同一批参数在不同允许误差下的筛选。范围是采样取值的分散程度，不是不确定性的上下界。}
\label{tab:tol}
\begin{tabular}{rrrrrrr}\toprule
$\epsilon$ (mm) & 被接受 & 到平面 & 过门口 & 先触地 & $y_g$ (m) & $z_g$ (m)\\\midrule
@@TOLROWS@@
\bottomrule\end{tabular}\end{table}

图~\ref{fig:separation}所用的两条轨迹取自 @@TOLMM@@ mm 下的被接受结果，在 MATLAB 数组中的编号为 @@IDXA@@ 和 @@IDXB@@。表~\ref{tab:pair}列出它们的参数。两者都在 @@TOLMM@@ mm 内重现两个位置，落点却相差 @@PAIRSEP@@ m。
\begin{table}[htbp]\centering\small
\caption{两条落点不同的符合观测轨迹。两者从出脚起都有横风。}
\label{tab:pair}
\begin{tabular}{lrrrrrrrr}\toprule
 & $w_y$ & $\rho$ & $C_D$ & $|v_0|$ & 最大误差(mm) & $y_g$ & $z_g$ & $t_g$\\\midrule
@@PAIRROWS@@
\bottomrule\end{tabular}
\par\smallskip\footnotesize 风速与初速单位 m/s；密度 kg/m$^3$；平面坐标 m；时间 s。
\end{table}

图~\ref{fig:separation}解释了为什么两次观测不能确定唯一落点：左图中，两条轨迹的六个坐标误差都落在允许范围内；右图中，它们在接近球门时逐渐分开。风和阻力条件变化后，反求出的初速度和侧向力也会变化。这些变化在短时间内能互相补偿，却不能保证后续位置仍然一致。
\begin{figure}[htbp]\centering
\includegraphics[width=.96\linewidth]{观测约束与轨迹分离.png}
\caption{两条符合观测的轨迹都通过六个坐标检查（左），但后程水平轨迹发生分离（右）。浅色带为声明的允许误差；位移按真实比例绘制。}
\label{fig:separation}
\end{figure}

图~\ref{fig:coupling}展示参数之间的联动。对这组观测，密度与阻力系数的乘积越大，通常需要更高的初速度（左）；横风变化时，计算出的侧向系数也会变化（中）；允许误差放宽后，能通过检查的参数组合随之增加（右）。
\begin{figure}[htbp]\centering
\includegraphics[width=\linewidth]{参数耦合与筛选敏感性.png}
\caption{参数耦合与筛选：三个风速切片中密度--阻力乘积与反求初速（左）；横风与被接受轨迹的反求 $b_z$（中）；允许误差与被接受数量（右）。所有数量都只针对列出的参数。}
\label{fig:coupling}
\end{figure}

图~\ref{fig:scatter}直接展示符合观测的轨迹在球门平面上的落点范围。颜色表示指定的阻力系数。点的密集程度由采样取值决定，不能据此判断哪一个落点更可能发生。
\begin{figure}[htbp]\centering
\includegraphics[width=.72\linewidth]{全程约束球门散点.png}
\caption{被接受轨迹在球门平面的落点。颜色表示列出的阻力系数。分散来自采样的参数取值，不是概率分布。}
\label{fig:scatter}
\end{figure}

另有一个交互界面允许横风 $[-15,15]$ m/s、密度 $[0.90,1.40]$ kg/m$^3$、$C_D\in[0.08,0.80]$；这只是范围更大的探索，与论文网格不同，筛选规则相同。

\section{计算检查与适用范围}
本节列出 MATLAB 内部做过的检查，以及这个结果的适用边界。

所有计算使用同一套方程和事件定义。三档步长只改变积分的最大步长，因此检验的是数值敏感性，而不是另一个模型。根据式~\eqref{eq:energy}推导的功与能关系，检查累积残差：
\begin{equation}
R_E(t)=E(t)-E(0)-\int_0^t(\mathbf F_D+\mathbf F_M)\cdot\mathbf v\,d\tau.
\end{equation}
事件检查确认触地、向前穿过 $x=11$ m 以及超时都被分别定位并记录。@@CHECKZH@@

这些检查检验的是程序与内部的能量平衡。没有独立的实测全程飞行，因此它们不能证明这只球、这块场地和这台相机服从该模型。报出的落点和范围只对列出的参数取值和声明的允许误差成立，不带概率，也不能保证搜索范围之外的情况。放宽允许误差会接受更多轨迹，但不会把采样得到的分散范围变成置信区间。

\section{结论}
本节回答问题。

仅有早期两个位置，在风和空气动力系数未知时，无法恢复唯一轨迹。对每一组列出的环境取值，仍可反求能重现两个位置的出脚条件；这些被接受的轨迹在球门平面上的落点并不相同。当前示例取风速 @@SPWY@@ m/s、密度 @@SPRHO@@ kg/m$^3$ 时，球心在 @@SELTEND@@ s 后穿过球门平面于 $(@@SELY@@,@@SELZ@@)$ m。在列出的参数中，@@COUNT@@ 组里有 @@ACCEPT@@ 组通过 @@TOLMM@@ mm 检查，落点覆盖 $y\in[@@YMIN@@,@@YMAX@@]$ m、$z\in[@@ZMIN@@,@@ZMAX@@]$ m，其中两条符合观测的轨迹相距 @@PAIRSEP@@ m。这些范围只描述采样网格，不是概率。更多观测时刻、实测风速和实测阻力系数可以缩小答案。

\begin{thebibliography}{9}
\bibitem{upc} University Physics Competition. 2026 contest, Problem B. \url{https://www.uphysicsc.com/2026contest.html}. 访问于 2026 年 11 月 7 日.
\bibitem{ifab} The IFAB. Laws of the Game, Law 2: The Ball. \url{https://www.theifab.com/laws/latest/the-ball/}. 访问于 2026 年 11 月 7 日.
\bibitem{asai} T. Asai, K. Seo, O. Kobayashi and R. Sakashita. Fundamental aerodynamics of the soccer ball. \emph{Sports Engineering} 10, 101--109 (2007). \url{https://doi.org/10.1007/BF02844207}.
\bibitem{nasa} NASA Glenn Research Center. Drag Equation. \url{https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/drag-equation/}.
\bibitem{ode} MathWorks. \texttt{ode45}: Solve nonstiff differential equations. \url{https://www.mathworks.com/help/matlab/ref/ode45.html}. 使用 MATLAB R2024a 实现.
\bibitem{lsq} MathWorks. \texttt{lsqnonlin}: Solve nonlinear least-squares problems. \url{https://www.mathworks.com/help/optim/ug/lsqnonlin.html}.
\end{thebibliography}

\smallskip\small
\section*{人工智能使用报告}
人工智能工具协助了代码生成、数值检查和文稿起草。所有采用的内容须由队员核验并负责。

\textbf{队号和完整 AI 记录待补，尚非提交终稿。}竞赛规则要求附上所用 AI 工具的完整输入与输出，该记录须在提交前补全并复核。
\end{document}
'''


# 第2/3节的推导在原有六节结构内展开，计算结果仍来自同一JSON。
sc=ver['linearization']['scales']
values.update({
 'UREL':fmt(sc['relativeSpeed'],2),'FD0':fmt(sc['dragForce'],2),'FM0':fmt(sc['sideForce'],2),
 'FG0':fmt(sc['weight'],2),'RE0':sci(sc['Reynolds'],2),'MA0':fmt(sc['Mach'],3),
 'PID':fmt(sc['dragChangeScale'],3),'PIM':fmt(sc['sideChangeScale'],3),'PIG':fmt(sc['gravityChangeScale'],3),
 'FBRATIO':fmt(100*sc['buoyancyToWeight'],2),'FBDISP':fmt(100*sc['buoyancyUnrefittedDisplacement'],2),
 'CORIOLIS':fmt(1000*sc['coriolisDisplacementScale'],2),'SENSCHECK':sci(ver['goalSensitivityMaxDifference'])})
PHYSICS_ZH=r"""
本节依次回答三个问题：球受到什么力，这些力怎样写成方程，哪些简化会影响最后的预测。先建立受力依据，再使用方程计算。
\subsection{研究对象、坐标和初始条件}
研究对象是已经离脚、尚未碰到地面或门柱的足球。取球心位置$\mathbf r=(x,y,z)$，速度$\mathbf v=\dot{\mathbf r}$，$x$轴指向球门，$z$轴向上。球半径$R=@@R@@$ m、质量$m=@@M@@$ kg；初始球心按题目固定，初速度待求：
\begin{equation}
\mathbf r(0)=(0,0,R),\qquad \mathbf v(0)=\mathbf v_0,\qquad
\mathbf g=(0,0,@@G@@)\ \mathrm{m\,s^{-2}}.
\end{equation}
质量和周长$2\pi R$符合IFAB球规格~\cite{ifab}，但不表示测得了本题足球的质量与材料。脚的作用已通过$\mathbf v_0$体现；球离脚后不再承受脚的推力。忽略形变后的球心平动满足牛顿第二定律：
\begin{equation}
m\dot{\mathbf v}=\mathbf F_g+\mathbf F_D+\mathbf F_M,\qquad \mathbf F_g=m\mathbf g.
\label{eq:newton}
\end{equation}
接下来分别建立阻力$\mathbf F_D$和侧向力$\mathbf F_M$，而不是直接从位置曲线猜加速度。
\subsection{从相对气流得到阻力}
空气动力取决于球相对空气的速度。设本次射门期间横风近似不变，$\mathbf w=(0,w_y,0)$，定义
\begin{equation}
\mathbf u=\mathbf v-\mathbf w,\quad U=\|\mathbf u\|,\quad
A=\pi R^2,\quad p_{\rm dyn}=\tfrac12\rho U^2.
\end{equation}
$p_{\rm dyn}$是动压，乘迎风面积得到力的尺度。把球面接缝、边界层和分离等影响收进无量纲阻力系数$C_D$，便有~\cite{nasa}
\begin{equation}
D=p_{\rm dyn}AC_D,\qquad
\mathbf F_D=-D\frac{\mathbf u}{U}
=-\tfrac12\rho AC_DU\mathbf u.
\label{eq:drag}
\end{equation}
负号规定阻力与相对气流反向；$\rho AU^2$的单位是N。式中并非既有“运动阻力”又额外加一份“风阻”：风通过$\mathbf u$改变同一项空气动力。动压给出量级，$C_D$仍需实验标定，不能只用伯努利方程或能量守恒把它推定为某个数。
\subsection{为什么使用有效侧向力}
旋转会改变球两侧的边界层与尾流，使气动力出现垂直于气流的分量。足球风洞实验支持这一机制，也显示系数与旋转、流态有关~\cite{asai}。用通常的力系数定义，侧向力大小为
\begin{equation}
F_M=\tfrac12\rho AU^2C_L.
\end{equation}
为了连接到可计算的向量形式，记$\boldsymbol\omega_\perp$为角速度垂直于$\mathbf u$的部分。若在有限速度范围内近似$C_L\simeq c_sR\|\boldsymbol\omega_\perp\|/U$，其中$c_s$为无量纲系数，并取方向沿$\boldsymbol\omega\times\mathbf u$，则
\begin{equation}
\mathbf F_M\simeq\tfrac12\rho Ac_sR(\boldsymbol\omega\times\mathbf u)
=m(\mathbf b\times\mathbf u),\qquad
\mathbf b=\frac{\rho Ac_sR}{2m}\boldsymbol\omega.
\label{eq:magnus}
\end{equation}
这一步是局部经验近似，不是风洞文献证明的普适线性定律。$\mathbf b$的单位是s$^{-1}$，把真实自旋与未知的气动效率$c_s$合在一起。两帧位置无法分别确定它们，因此程序只求有效$\mathbf b$，不把它解释成实测角速度。

本文取$\mathbf b=(0,b_y,b_z)$，使主要沿$x$方向运动的球，其横向和竖向偏转分别由$b_z$和$b_y$描述。令$b_x=0$是减少未知量的假设，并不说明真实球的这个分量必为零。若要预测自旋本身，还需转动惯量和空气力矩去求$\mathbf I\dot{\boldsymbol\omega}=\boldsymbol\tau$；题目没有这些数据。保持$\mathbf b$恒定只代表在约0.4s内使用恒定有效参数，其误差须留在模型局限中。
\subsection{合成完整运动方程}
将式~\eqref{eq:drag}、\eqref{eq:magnus}代入牛顿方程，并令$K=\rho A/(2m)$，得到程序实际积分的六个一阶方程：
\begin{align}
\dot{\mathbf r}&=\mathbf v,\qquad
\dot{\mathbf v}=\mathbf g-KC_DU\mathbf u+\mathbf b\times\mathbf u,\label{eq:model}\\
\dot v_x&=-KC_DUv_x+b_yv_z-b_z(v_y-w_y),\nonumber\\
\dot v_y&=-KC_DU(v_y-w_y)+b_zv_x,\nonumber\\
\dot v_z&=-g-KC_DUv_z-b_yv_x,\qquad g=\|\mathbf g\|.\nonumber
\end{align}
分量式也说明三种运动相互影响：横风不仅改变$y$方向阻力，$v_y-w_y$还进入$x$方向侧力。改变某一参数后，必须重新计算整个向量状态。
\begin{table}[H]\centering\small
\caption{本节符号及其角色。}\begin{tabular}{lll}\toprule
符号 & 含义 & 单位或来源\\\midrule
$\mathbf r,\mathbf v,\mathbf u$ & 位置、地面速度、相对空气速度 & m；m/s；m/s\\
$m,R,A$ & 质量、半径、迎风面积 & kg；m；m$^2$\\
$\rho,w_y,C_D$ & 密度、横风、阻力系数 & 指定情景，不是实测\\
$\mathbf v_0,b_y,b_z$ & 初速度和有效侧力参数 & m/s；s$^{-1}$，待求\\
$K,C_L,C_\perp$ & 气动尺度、升力与侧力系数 & m$^{-1}$；无量纲\\
\bottomrule\end{tabular}\end{table}
\subsection{用量级判断哪些作用重要}
为判断是否能忽略某项，取当前示例出脚相对气流速度$U_0=@@UREL@@$ m/s、飞行时间$T=@@SELTEND@@$ s。初始重力、阻力、侧力大小分别约为$@@FG0@@$、$@@FD0@@$、$@@FM0@@$ N，气动力与重力是同一量级。进一步比较它们在飞行时间内改变速度的相对大小：
\begin{equation}
\Pi_D=KC_DU_0T=@@PID@@,\quad
\Pi_M=\|\mathbf b\|T=@@PIM@@,\quad
\Pi_g=\frac{gT}{U_0}=@@PIG@@.
\end{equation}
这些量均不远小于1，因此三项都保留。采用示例空气黏度$\mu=1.81\times10^{-5}$ Pa s和声速$c=343$ m/s，可估计
\begin{equation}
\mathrm{Re}=\frac{2\rho RU_0}{\mu}\simeq@@RE0@@,\qquad
\mathrm{Ma}=\frac{U_0}{c}\simeq@@MA0@@.
\end{equation}
声速和黏度为量级估算所用的假设，不是场地测量。较小的马赫数支持忽略压缩性；雷诺数说明阻力系数与流态有关，所以恒定$C_D$是有限工况下的简化，不能因此认定它已被标定。

其余简化也要估计影响。空气浮力与重力之比约为$\rho(4\pi R^3/3)/m=@@FBRATIO@@\%$；在不重新拟合出脚条件的估算中，忽略它可产生约$@@FBDISP@@$ cm的竖向位移差。因此本文不能把毫米级的计算残差说成毫米级的真实准确度。取地球自转角速度$\Omega=7.292115\times10^{-5}$ s$^{-1}$，其导致的位移量级$\Omega\|\mathbf v_0\|T^2$约为$@@CORIOLIS@@$ mm，低于本例允许误差。浮力、非恒定风和变化的气动系数仍会贡献模型误差，更多早期观测才能帮助检验它们。
\subsection{有风时怎样检查功与能}
令$E=\tfrac12m\|\mathbf v\|^2+mgz$。将牛顿方程与地面系速度点乘，重力势能项抵消后得到
\begin{align}
\frac{dE}{dt}&=(\mathbf F_D+\mathbf F_M)\cdot\mathbf v\label{eq:energy}\\
&=-\tfrac12\rho AC_DU^3+(\mathbf F_D+\mathbf F_M)\cdot\mathbf w.\nonumber
\end{align}
推导中用到了$\mathbf F_M\cdot\mathbf u=0$，但$\mathbf v=\mathbf u+\mathbf w$。所以侧向力对相对气流的瞬时功率为零，有风时对地面系的球仍可做功；不能强求球的地面系机械能一直下降。后面的功与能检查正是积分式~\eqref{eq:energy}，用于检验程序是否把力的方向和符号算对。
"""
PHYSICS_EN=r"""
This section builds the force law in order: define the moving system, relate airflow to force, combine the forces, and check the size and limits of the approximations.
\subsection{System, coordinates and initial conditions}
The system is the ball after release from the foot and before contact with the ground or a post. Its centre is $\mathbf r=(x,y,z)$ and $\mathbf v=\dot{\mathbf r}$; $x$ points toward the goal and $z$ upward. With radius $R=@@R@@$ m and mass $m=@@M@@$ kg,
\begin{equation}
\mathbf r(0)=(0,0,R),\quad \mathbf v(0)=\mathbf v_0,\quad
\mathbf g=(0,0,@@G@@)\ \mathrm{m\,s^{-2}}.
\end{equation}
The chosen mass and circumference $2\pi R$ satisfy IFAB specifications~\cite{ifab}; they are not measurements of this particular ball. The kick is represented by the unknown launch velocity. After release there is no continuing force from the foot. Newton's second law for the centre is
\begin{equation}
m\dot{\mathbf v}=\mathbf F_g+\mathbf F_D+\mathbf F_M,\qquad\mathbf F_g=m\mathbf g.
\label{eq:newton}
\end{equation}
The next steps establish drag $\mathbf F_D$ and the transverse force $\mathbf F_M$ instead of prescribing a fitted acceleration curve.
\subsection{Drag from relative airflow}
For steady lateral wind $\mathbf w=(0,w_y,0)$ define
\begin{equation}
\mathbf u=\mathbf v-\mathbf w,\quad U=\|\mathbf u\|,\quad A=\pi R^2,\quad p_{\rm dyn}=\tfrac12\rho U^2.
\end{equation}
Dynamic pressure multiplied by frontal area gives a force scale. The dimensionless drag coefficient includes effects of the surface, boundary layer and flow separation~\cite{nasa}. Hence
\begin{equation}
D=p_{\rm dyn}AC_D,\qquad\mathbf F_D=-D\frac{\mathbf u}{U}
=-\tfrac12\rho AC_DU\mathbf u.
\label{eq:drag}
\end{equation}
The sign makes drag oppose relative motion, and $\rho AU^2$ has units of newtons. Wind changes this same force through $\mathbf u$; it is not an additional copy of drag. Dynamic pressure supplies a scale, while $C_D$ requires calibration. Bernoulli's equation or conservation of energy alone cannot determine its numerical value.
\subsection{Meaning of the effective transverse force}
Rotation changes the boundary layer and wake and can produce a force perpendicular to airflow. Football wind-tunnel measurements support this mechanism while showing dependence on rotation and flow regime~\cite{asai}. Write its magnitude as
\begin{equation}
F_M=\tfrac12\rho AU^2C_L.
\end{equation}
Let $\boldsymbol\omega_\perp$ denote the angular-velocity component perpendicular to $\mathbf u$. Over a limited range, assume $C_L\simeq c_sR\|\boldsymbol\omega_\perp\|/U$ with dimensionless slope $c_s$ and direction along $\boldsymbol\omega\times\mathbf u$. Substitution gives
\begin{equation}
\mathbf F_M\simeq\tfrac12\rho Ac_sR(\boldsymbol\omega\times\mathbf u)
=m(\mathbf b\times\mathbf u),\qquad
\mathbf b=\frac{\rho Ac_sR}{2m}\boldsymbol\omega.
\label{eq:magnus}
\end{equation}
The linear coefficient relation is a modelling approximation, not a universal law established by the cited experiment. The fitted $\mathbf b$, in s$^{-1}$, combines spin with unknown aerodynamic efficiency. Two positions cannot identify these separately, so $\mathbf b$ is not reported as a measured angular velocity.

We set $\mathbf b=(0,b_y,b_z)$: for predominantly forward motion, $b_z$ controls lateral bending and $b_y$ vertical bending. Fixing $b_x=0$ reduces the unknowns but does not prove that the actual component vanishes. Resolving spin would also require inertia and aerodynamic torque in $\mathbf I\dot{\boldsymbol\omega}=\boldsymbol\tau$, which are not given. Constant $\mathbf b$ is a short-flight approximation whose error remains part of the model limitations.
\subsection{The complete equations used for integration}
Substituting the two force expressions into Newton's law and defining $K=\rho A/(2m)$ yields six first-order equations:
\begin{align}
\dot{\mathbf r}&=\mathbf v,\qquad\dot{\mathbf v}=\mathbf g-KC_DU\mathbf u+\mathbf b\times\mathbf u,\label{eq:model}\\
\dot v_x&=-KC_DUv_x+b_yv_z-b_z(v_y-w_y),\nonumber\\
\dot v_y&=-KC_DU(v_y-w_y)+b_zv_x,\nonumber\\
\dot v_z&=-g-KC_DUv_z-b_yv_x,\qquad g=\|\mathbf g\|.\nonumber
\end{align}
These components explain coupling: wind changes both lateral drag and the $x$ component of the transverse force. Changing one input therefore requires recomputing the complete state.
\begin{table}[H]\centering\small
\caption{Symbols used in the force derivation.}\begin{tabular}{lll}\toprule
Symbol & Meaning & Unit or role\\\midrule
$\mathbf r,\mathbf v,\mathbf u$ & Centre, ground velocity, relative velocity & m; m s$^{-1}$; m s$^{-1}$\\
$m,R,A$ & Mass, radius, frontal area & kg; m; m$^2$\\
$\rho,w_y,C_D$ & Density, wind, drag coefficient & prescribed scenarios\\
$\mathbf v_0,b_y,b_z$ & Launch and effective transverse coefficients & fitted; m s$^{-1}$; s$^{-1}$\\
$K,C_L,C_\perp$ & Aerodynamic scale and force coefficients & m$^{-1}$; dimensionless\\
\bottomrule\end{tabular}\end{table}
\subsection{Scale checks and approximation limits}
For the current example, launch relative speed is $U_0=@@UREL@@$ m s$^{-1}$ and flight time $T=@@SELTEND@@$ s. The initial weight, drag and transverse-force magnitudes are approximately $@@FG0@@$, $@@FD0@@$ and $@@FM0@@$ N. Aerodynamic forces are comparable to weight. Their relative velocity changes over one flight are estimated by
\begin{equation}
\Pi_D=KC_DU_0T=@@PID@@,\quad\Pi_M=\|\mathbf b\|T=@@PIM@@,\quad
\Pi_g=\frac{gT}{U_0}=@@PIG@@.
\end{equation}
None is very small, so all three accelerations are retained. With an assumed viscosity $\mu=1.81\times10^{-5}$ Pa s and sound speed $c=343$ m s$^{-1}$,
\begin{equation}
\mathrm{Re}=\frac{2\rho RU_0}{\mu}\simeq@@RE0@@,\qquad\mathrm{Ma}=\frac{U_0}{c}\simeq@@MA0@@.
\end{equation}
Viscosity and sound speed are scale assumptions, not site measurements. The small Mach number supports neglecting compressibility. Reynolds number makes the dependence of $C_D$ on flow regime relevant; it does not calibrate a constant $C_D$ for this ball.

Neglected terms also deserve scale estimates. The buoyancy-to-weight ratio is $\rho(4\pi R^3/3)/m\simeq@@FBRATIO@@\%$. Without refitting launch conditions, omitting buoyancy gives a vertical-displacement scale of about $@@FBDISP@@$ cm over this flight. Thus millimetre fitting residuals cannot be called millimetre physical accuracy. With Earth's rotation rate $\Omega=7.292115\times10^{-5}$ s$^{-1}$, the displacement scale $\Omega\|\mathbf v_0\|T^2$ of roughly $@@CORIOLIS@@$ mm, below the chosen coordinate tolerance. Buoyancy, changing wind and changing coefficients remain sources of model error; additional observations are needed to test them.
\subsection{Work and energy in moving air}
For $E=\tfrac12m\|\mathbf v\|^2+mgz$, taking the scalar product of Newton's law with ground-frame velocity cancels the gravitational potential term:
\begin{align}
\frac{dE}{dt}&=(\mathbf F_D+\mathbf F_M)\cdot\mathbf v\label{eq:energy}\\
&=-\tfrac12\rho AC_DU^3+(\mathbf F_D+\mathbf F_M)\cdot\mathbf w.\nonumber
\end{align}
Here $\mathbf F_M\cdot\mathbf u=0$, but $\mathbf v=\mathbf u+\mathbf w$. The transverse force is perpendicular to relative airflow, yet it can do ground-frame work when wind is present. Ground-frame mechanical energy therefore need not decrease monotonically. Integrating Eq.~\eqref{eq:energy} supplies the work-energy check used later.
"""
DETAILS_ZH=r"""
\subsection{为什么不同参数能通过同样的观测}
前面的数值筛选可以进一步用方程说明。把全部未知量写为
\begin{equation}
p=(v_{0x},v_{0y},v_{0z},b_y,b_z,w_y,\rho,C_D),\qquad
\mathbf f(p)=\begin{bmatrix}\mathbf r(t_1;p)\\\mathbf r(t_2;p)\end{bmatrix}.
\end{equation}
这里$p$有8个分量，而观测只提供6个位置数。更直接的例子是：阻力中只出现$\rho C_D$；保持这个乘积、风和有效$\mathbf b$不变，力场和轨迹就完全相同。这说明仅靠观测不能把密度和阻力系数分别测出来。若横风也改变，就可能用不同的出脚条件使两次位置仍落在允许误差内，后续轨迹却不同。

在一个已经求出的参数附近，记$J=\partial\mathbf f/\partial p$，则
\begin{equation}
\delta\mathbf f\simeq J\delta p,\qquad
\|\mathbf f(p)-\mathbf r_{\rm obs}\|_\infty\le\epsilon.
\end{equation}
使$J\delta p$很小的参数变化，在早期数据上难以区分；这给“参数互相补偿”一个明确含义。它只是局部解释，最终是否符合观测仍须回代完整非线性方程。最小二乘减少的是误差平方和，不保证每个坐标都过关，所以必须保留前面的逐坐标检查。某组最小二乘结果没过关，也不证明这个环境下所有参数都不可能过关。
\subsection{从连续飞行确定球门位置}
球门点不是在一个随意选定的时间读取，而是满足
\begin{equation}
x(t_g;p)=11,\quad v_x(t_g)>0,\qquad
(y_g,z_g)=(y(t_g;p),z(t_g;p)).
\label{eq:goal}
\end{equation}
同时检查下降触地事件$z(t)-R=0,\ v_z(t)<0$，先发生哪个就记录哪个。在球门中心平面，一个简单的有限半径净空检查为
\begin{equation}
|y_g|<3.66-R,\qquad R\le z_g<2.44-R.
\end{equation}
这里严格采用题目球心穿过$x=11$的约定；这个局部几何检查不代替球与门柱的完整碰撞计算，也不把未到达球门的触地点投影成进球点。
\subsection{参数差异怎样传到球门点}
以下只讨论相对速度非零、仍先到达球门且向前穿越的小幅参数变化。不仅$y,z$会变，穿越时刻$t_g$也会变。将式~\eqref{eq:goal}对$p$求导，得到
\begin{equation}
\frac{\partial t_g}{\partial p}=-\frac{\partial x/\partial p}{v_x(t_g)},\qquad
G=\frac{\partial(y_g,z_g)}{\partial p}
=S_{yz}-\frac{1}{v_x(t_g)}\begin{bmatrix}v_y\\v_z\end{bmatrix}S_x.
\label{eq:goalderivative}
\end{equation}
$S_x,S_{yz}$是在同一固定时间$t_g$对位置求导；最后一项补偿“到达时间也变了”的影响。漏掉这一项，就会把固定时刻的位置误差误当成球门平面误差。

这些导数也能沿同一方程计算。令$s=(\mathbf r,\mathbf v)$、$\dot s=F(s,p)$、$S=\partial s/\partial p$，直接求导便有
\begin{equation}
\dot S=F_sS+F_p,\quad
F_s=\begin{bmatrix}0&I\\0&D\end{bmatrix},\qquad
D=-KC_D\left(UI+\frac{\mathbf u\mathbf u^T}{U}\right)+[\mathbf b]_\times.
\end{equation}
$[\mathbf b]_\times\mathbf u=\mathbf b\times\mathbf u$；$S(0)$的位置行全为零，速度对初速度的三个导数构成单位阵，其余为零。不同参数的直接影响例如
\begin{equation}
\frac{\partial\dot{\mathbf v}}{\partial b_y}=(u_z,0,-u_x)^T,\quad
\frac{\partial\dot{\mathbf v}}{\partial b_z}=(-u_y,u_x,0)^T,\quad
\frac{\partial\dot{\mathbf v}}{\partial w_y}=-D\mathbf e_y.
\end{equation}
这些式子分别接回侧力方向和相对气流定义，因此导数与受力模型是一致的。MATLAB另行积分这组变化方程，并用原方程的小幅参数扰动检查式~\eqref{eq:goalderivative}；球门导数的最大差异为$@@SENSCHECK@@$（按所用SI参数单位）。它检查的是公式与程序是否一致，不给参数赋予概率；正式轨迹仍由原来的非线性运动方程求出。
"""
DETAILS_EN=r"""
\subsection{How different parameters can meet the same observations}
The screening has a mathematical explanation. Collect the unknowns and observation map as
\begin{equation}
p=(v_{0x},v_{0y},v_{0z},b_y,b_z,w_y,\rho,C_D),\qquad
\mathbf f(p)=\begin{bmatrix}\mathbf r(t_1;p)\\\mathbf r(t_2;p)\end{bmatrix}.
\end{equation}
There are eight entries in $p$ and six observed coordinates. More directly, drag contains only $\rho C_D$: changing density and drag inversely while keeping wind and effective $\mathbf b$ fixed preserves the complete force field. The two values cannot be measured separately from these positions. Changing wind can instead be compensated by different launch conditions within the allowed coordinate error, while later paths differ.

Near a fitted point, let $J=\partial\mathbf f/\partial p$. Then
\begin{equation}
\delta\mathbf f\simeq J\delta p,\qquad
\|\mathbf f(p)-\mathbf r_{\rm obs}\|_\infty\le\epsilon.
\end{equation}
Changes giving small $J\delta p$ are difficult to distinguish in the early data. This gives a precise local meaning to parameter compensation. Acceptance still requires substitution into the nonlinear equations. Least squares minimizes a sum, not the largest coordinate error, so componentwise screening remains essential. Failure of a fitted choice does not establish that all possible launch conditions in that environment fail.
\subsection{Locating the goal point during continuous flight}
We do not read the endpoint at an arbitrary common time. The crossing satisfies
\begin{equation}
x(t_g;p)=11,\quad v_x(t_g)>0,\qquad
(y_g,z_g)=(y(t_g;p),z(t_g;p)).
\label{eq:goal}
\end{equation}
Descending ground contact $z(t)-R=0,\ v_z(t)<0$ competes with this event, and the first event is recorded. At the centre plane a simple finite-radius opening check is
\begin{equation}
|y_g|<3.66-R,\qquad R\le z_g<2.44-R.
\end{equation}
The predicted point follows the question's centre-crossing convention $x=11$. This local clearance check is not a complete post-collision calculation. A ground-first endpoint is never projected onto the goal as a crossing.
\subsection{How parameter changes reach the goal point}
For small changes that keep relative speed nonzero and preserve a forward goal-first crossing, the crossing time changes along with $y,z$. Differentiating Eq.~\eqref{eq:goal} gives
\begin{equation}
\frac{\partial t_g}{\partial p}=-\frac{\partial x/\partial p}{v_x(t_g)},\qquad
G=\frac{\partial(y_g,z_g)}{\partial p}
=S_{yz}-\frac{1}{v_x(t_g)}\begin{bmatrix}v_y\\v_z\end{bmatrix}S_x.
\label{eq:goalderivative}
\end{equation}
$S_x,S_{yz}$ are position derivatives at the same fixed time $t_g$. The second term corrects for the changed arrival time. Omitting it confuses a fixed-time position change with a change in the goal-plane crossing.

To compute these derivatives consistently, write $s=(\mathbf r,\mathbf v)$, $\dot s=F(s,p)$ and $S=\partial s/\partial p$. Differentiation yields
\begin{equation}
\dot S=F_sS+F_p,\quad
F_s=\begin{bmatrix}0&I\\0&D\end{bmatrix},\qquad
D=-KC_D\left(UI+\frac{\mathbf u\mathbf u^T}{U}\right)+[\mathbf b]_\times.
\end{equation}
Here $[\mathbf b]_\times\mathbf u=\mathbf b\times\mathbf u$. Initially the position derivatives are zero, the velocity derivatives with respect to launch velocity form the identity, and the other entries are zero. Direct parameter effects include
\begin{equation}
\frac{\partial\dot{\mathbf v}}{\partial b_y}=(u_z,0,-u_x)^T,\quad
\frac{\partial\dot{\mathbf v}}{\partial b_z}=(-u_y,u_x,0)^T,\quad
\frac{\partial\dot{\mathbf v}}{\partial w_y}=-D\mathbf e_y.
\end{equation}
These expressions follow directly from the force directions and the definition of relative velocity. MATLAB integrates the variation equations as a diagnostic and checks Eq.~\eqref{eq:goalderivative} against small parameter perturbations of the original flight equations. The largest derivative discrepancy is $@@SENSCHECK@@$ in the adopted SI parameter units. This tests the implementation without assigning probabilities. The reported trajectories continue to be calculated from the nonlinear flight equations.
"""

# DERIVATION_BLOCKS
# PHYSICS_EN/PHYSICS_ZH、DETAILS_EN/DETAILS_ZH 由主审定义；必须先替换占位符再 render。
for lang, template in [('英文', template_en), ('中文', template_zh)]:
    path = ROOT / f'latex/AI/B题{lang}论文.tex'
    filled = template.replace('@@PHYSICS@@', PHYSICS_EN if lang == '英文' else PHYSICS_ZH)
    filled = filled.replace('@@INFERENCE_DETAILS@@', DETAILS_EN if lang == '英文' else DETAILS_ZH)
    text=render(filled).replace('[htbp]','[H]')
    text=text.replace(r'\begin{thebibliography}{9}',r'\begingroup\small\raggedright'+'\n'+r'\begin{thebibliography}{9}').replace(r'\end{document}',r'\endgroup'+'\n'+r'\end{document}')
    text=text.replace(r'\section{Conclusion}',r'\clearpage'+'\n'+r'\section{Conclusion}').replace(r'\section{结论}',r'\clearpage'+'\n'+r'\section{结论}')
    path.write_text(text, encoding='utf-8')

print('中英文单文件论文已由正式 JSON 重新生成。')
