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
    en_check.append(f"The work-energy balance closes to ${sci(ew)}$ J.")
    zh_check.append(f"功与能平衡的残差为 ${sci(ew)}$ J。")
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
This section writes the force law used from the kick to the first terminal event.

The ball is a sphere of radius $R=@@R@@$ m and mass $m=@@M@@$ kg. The air is steady with lateral wind $\mathbf w=(0,w_y,0)$, constant density $\rho$, constant drag coefficient $C_D$, and a constant effective sideways coefficient $\mathbf b=(0,b_y,b_z)$. We use
\begin{equation}
\dot{\mathbf r}=\mathbf v,\qquad
\dot{\mathbf v}=\mathbf g-KC_D|\mathbf u|\mathbf u+\mathbf b\times\mathbf u,
\qquad \mathbf u=\mathbf v-\mathbf w,
\qquad K=\frac{\rho\pi R^2}{2m},
\label{eq:model}
\end{equation}
with $\mathbf g=(0,0,-@@G@@)$ m s$^{-2}$. The drag term is quadratic and opposes the motion relative to the air. The cross term is perpendicular to $\mathbf u$. The coefficient $\mathbf b$ has units s$^{-1}$ and combines rotation and lift efficiency; it is not the angular velocity, and we do not claim to recover the true spin. Setting $b_x=0$ removes a direction that the two positions constrain weakly. Constant coefficients, steady wind, no spin decay and no drag crisis are assumptions, not measurements; they are checked only as a model, not calibrated to this ball. The force form follows standard football aerodynamics~\cite{asai}.

\begin{table}[htbp]\centering
\caption{Symbols, units and roles. Listed environmental values are scenario inputs; launch conditions are inferred.}
\begin{tabular}{lll}\toprule
Symbol & Meaning & Unit or role\\\midrule
$\mathbf r,\mathbf v$ & Ball centre and velocity & m; m s$^{-1}$\\
$\mathbf u=\mathbf v-\mathbf w$ & Velocity relative to the air & m s$^{-1}$\\
$\mathbf w=(0,w_y,0)$ & Steady lateral wind & listed scenario input\\
$\rho$ & Air density & listed scenario input\\
$C_D$ & Drag coefficient & listed scenario input\\
$\mathbf b=(0,b_y,b_z)$ & Effective sideways coefficient & s$^{-1}$; inferred\\
$K$ & $\rho\pi R^2/(2m)$ & m$^{-1}$\\
$\epsilon$ & Allowed coordinate error & m\\
$C_\perp$ & Side-force coefficient $|\mathbf b\times\mathbf u|/(K|\mathbf u|^2)$ & dimensionless\\
$(t_g,y_g,z_g)$ & First crossing of $x=11$ m & s; m; m\\
\bottomrule\end{tabular}
\end{table}

The fixed mass and radius satisfy the IFAB ball specifications~\cite{ifab}. Density, wind and drag values are listed scenarios, not weather measurements. We do not model the goalkeeper, post collisions, deformation during contact, or flight after the first terminal event. The ball starts at the given centre with no vertical offset beyond the radius.

Over the roughly 0.4 s between kick and goal, the three terms can all change the path. Gravity bends the flight downward throughout. Drag grows with the square of the speed relative to the air, so a tailwind or headwind changes the deceleration as well as the direction of the relative flow. The sideways term grows with speed and with the fitted coefficients, and it bends the path without doing work, since $\mathbf b\times\mathbf u$ is perpendicular to $\mathbf u$. A model with gravity alone cannot connect the early positions to a late one, but the data do not measure the size of either aerodynamic term. That gap is why the paper treats the environmental values as listed inputs rather than known quantities.

\section{Finding parameters from the observations}
This section explains how the two positions become parameter choices.

For each listed triple $\eta=(w_y,\rho,C_D)$ we fit $q=(v_{0x},v_{0y},v_{0z},b_y,b_z)$ by bounded nonlinear least squares,
\begin{equation}
q^*(\eta)=\arg\min_{q\in\mathcal B}\sum_{i=1}^{2}
\|\mathbf r(t_i;q,\eta)-\mathbf r_i\|^2,
\label{eq:fit}
\end{equation}
using \texttt{lsqnonlin}~\cite{lsq} from two starting points. The trajectory in the fit is integrated with RK4 at the recorded step. A fit is checked on all six coordinates,
\begin{equation}
|r_j(t_i;q,\eta)-r_{ij}|\le\epsilon,
\qquad i=1,2,\quad j=x,y,z,
\label{eq:tol}
\end{equation}
at $\epsilon=2,@@TOLMM@@,10,20$ mm; @@TOLMM@@ mm is the default. The allowed error is a stated scenario, not a measured camera accuracy.

The search uses @@NWIND@@ wind values on $[@@WMIN@@,@@WMAX@@]$ m s$^{-1}$, @@NDENS@@ density values $\{@@DENSITIES@@\}$ kg m$^{-3}$, and @@NCD@@ drag values on $[@@CDMIN@@,@@CDMAX@@]$, giving @@COUNT@@ parameter choices. These values are example settings, not measured distributions. Only the listed parameter values are searched, and the ranges in the figures do not guarantee that all possible cases are included.

Before a fit is called accepted it also passes a physical screening. The launch must point forward and not downward, the initial speed must lie between 10 and 40 m s$^{-1}$, and along the sampled trajectory the side-force coefficient must satisfy $C_\perp\le0.50$. The ball must stay above the ground before the terminal event, and the event must be located within 2 s. An accepted fit must pass both the six-coordinate test and this screening; observation fit and physical screening are stored separately. The speed range and the 0.50 cap are conservative model screens, not universal limits for players or footballs. We do not force the crossing into the goal: a miss is a recorded outcome, and the goal is used only to describe where the ball goes.

Two further points explain why the result is a set rather than one shot. First, the fitting problem has more unknowns than the six coordinates. Changing the wind or the drag changes the air force, and a different launch velocity can make up for much of that change at the two observed times. The fit therefore finds different launch conditions for different listed environments, even though all of them look the same at $t=0.05$ and $t=0.10$ s. Second, the two observations constrain the early path only. A small difference in the fitted sideways coefficient can build into a larger sideways displacement over the remaining flight, because the force keeps acting after the last observation. The accepted set is the result of searching the listed values, not a proof that nothing outside the grid exists. Only the listed parameter values are searched, and the ranges in the figures do not guarantee that all possible cases are included.

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

All calculations use the same equations and event definitions. The three step sizes differ only in the maximum step of the integration, so they measure numerical sensitivity rather than a new model. The work-energy check compares the change of mechanical energy with the work of drag and the sideways force,
\begin{equation}
\frac{d}{dt}\left(\frac12m|\mathbf v|^2+mgz\right)
=(\mathbf F_D+\mathbf F_M)\cdot\mathbf v.
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
本节写出从出脚到第一个终止事件所使用的受力方程。

把球看成半径 $R=@@R@@$ m、质量 $m=@@M@@$ kg 的球体。空气为恒定横风 $\mathbf w=(0,w_y,0)$，密度 $\rho$、阻力系数 $C_D$ 和有效侧向系数 $\mathbf b=(0,b_y,b_z)$ 都取常数。方程写成
\begin{equation}
\dot{\mathbf r}=\mathbf v,\qquad
\dot{\mathbf v}=\mathbf g-KC_D|\mathbf u|\mathbf u+\mathbf b\times\mathbf u,
\qquad \mathbf u=\mathbf v-\mathbf w,
\qquad K=\frac{\rho\pi R^2}{2m},
\label{eq:model}
\end{equation}
其中 $\mathbf g=(0,0,-@@G@@)$ m/s$^2$。阻力项是二次项，方向与相对空气速度相反；叉乘项与 $\mathbf u$ 垂直。$\mathbf b$ 的单位是 s$^{-1}$，它把转动和升力效率合并在一起，并不是角速度，我们也不声称恢复真实自旋。令 $b_x=0$，因为两次位置对这个方向的约束很弱。系数恒定、风恒定、自旋不衰减、阻力危机不出现，都是模型假设而非测量结果；它们只作为模型被检查，没有对这只球做标定。受力形式沿用足球空气动力学的常用写法~\cite{asai}。

\begin{table}[htbp]\centering
\caption{符号、单位与角色。列出的环境值是情景输入；出脚条件是反求量。}
\begin{tabular}{lll}\toprule
符号 & 含义 & 单位或角色\\\midrule
$\mathbf r,\mathbf v$ & 球心位置与速度 & m；m/s\\
$\mathbf u=\mathbf v-\mathbf w$ & 相对空气速度 & m/s\\
$\mathbf w=(0,w_y,0)$ & 恒定横风 & 列出的情景输入\\
$\rho$ & 空气密度 & 列出的情景输入\\
$C_D$ & 阻力系数 & 列出的情景输入\\
$\mathbf b=(0,b_y,b_z)$ & 有效侧向系数 & s$^{-1}$；反求量\\
$K$ & $\rho\pi R^2/(2m)$ & m$^{-1}$\\
$\epsilon$ & 逐坐标允许误差 & m\\
$C_\perp$ & 侧力系数 $|\mathbf b\times\mathbf u|/(K|\mathbf u|^2)$ & 无量纲\\
$(t_g,y_g,z_g)$ & 首次穿过 $x=11$ m & s；m；m\\
\bottomrule\end{tabular}
\end{table}

固定的质量和半径满足 IFAB 的球规格~\cite{ifab}。密度、风速和阻力系数都是列出的情景取值，不是现场测量。模型不包含守门员、门柱碰撞、触球时的形变以及第一个终止事件之后的飞行。球从题给球心出发，除半径外没有额外的竖直偏移。

从出脚到球门大约 0.4 s，三项力都会改变轨迹。重力全程把轨迹向下弯；阻力随相对空气速度的平方增大，顺风或逆风既改变减速大小，也改变相对气流方向；侧向力随速度和拟合系数增大，它只改变方向而不做功，因为 $\mathbf b\times\mathbf u$ 与 $\mathbf u$ 垂直。只有重力的模型不足以把早期位置接到后期位置，但观测也没有测出两项气动力的大小。正因为如此，后文把环境取值当作列出的输入，而不是已知量。

\section{用观测求参数}
本节说明两次位置如何变成一组组参数取值。

对每一组列出的 $\eta=(w_y,\rho,C_D)$，用有界非线性最小二乘反求 $q=(v_{0x},v_{0y},v_{0z},b_y,b_z)$：
\begin{equation}
q^*(\eta)=\arg\min_{q\in\mathcal B}\sum_{i=1}^{2}
\|\mathbf r(t_i;q,\eta)-\mathbf r_i\|^2,
\label{eq:fit}
\end{equation}
使用 \texttt{lsqnonlin}~\cite{lsq}，并从两个初值出发。拟合中的轨迹用记录步长的 RK4 积分。随后逐个检查六个坐标：
\begin{equation}
|r_j(t_i;q,\eta)-r_{ij}|\le\epsilon,
\qquad i=1,2,\quad j=x,y,z,
\label{eq:tol}
\end{equation}
分别在 $\epsilon=2,@@TOLMM@@,10,20$ mm 下筛选，默认取 @@TOLMM@@ mm。允许误差是明确声明的情景，不是测量得到的相机精度。

搜索使用 $[@@WMIN@@,@@WMAX@@]$ m/s 上的 @@NWIND@@ 个风速、$\{@@DENSITIES@@\}$ kg/m$^3$ 的 @@NDENS@@ 个密度值、以及 $[@@CDMIN@@,@@CDMAX@@]$ 上的 @@NCD@@ 个阻力系数，共 @@COUNT@@ 组参数。这些是示例设置，不是实测分布。

一组参数还要通过物理筛选才算被接受：出脚方向向前且不向下，初速度在 10 到 40 m/s 之间，沿采样轨迹的侧力系数满足 $C_\perp\le0.50$；球在终止事件前必须位于地面以上，终止事件要在 2 s 内定位到。观测拟合与物理筛选分别记录，只有两项都通过才被接受。速度范围和 0.50 上限是保守的模型筛选，不是运动员或足球的普适极限。我们不强迫球进门：射偏也是记录下来的结果，球门只用来描述球飞到哪里。

结果之所以是一组轨迹而不是一次射门，有两个原因。第一，未知量比六个坐标多。改变风速或阻力就改变了空气力，不同的初速度又能在两个观测时刻上补偿其中的大部分差别，于是不同的环境取值会反求出不同的出脚条件，而在 $t=0.05$ 与 $t=0.10$ s 处看起来一样。第二，两次观测只约束前段。拟合出的侧向系数只要有一点不同，在最后观测之后仍继续作用，落到球门平面时就会放大成明显的侧向差距。被接受的集合是搜索列出的参数得到的结果，并不证明网格之外没有别的解。只搜索列出的参数，图中的范围不能保证包含所有可能情况。

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

所有计算使用同一套方程和事件定义。三档步长只改变积分的最大步长，因此检验的是数值敏感性，而不是另一个模型。功与能检查比较机械能的变化与阻力、侧向力所做的功：
\begin{equation}
\frac{d}{dt}\left(\frac12m|\mathbf v|^2+mgz\right)
=(\mathbf F_D+\mathbf F_M)\cdot\mathbf v.
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
\bibitem{ode} MathWorks. \texttt{ode45}: Solve nonstiff differential equations. \url{https://www.mathworks.com/help/matlab/ref/ode45.html}. 使用 MATLAB R2024a 实现.
\bibitem{lsq} MathWorks. \texttt{lsqnonlin}: Solve nonlinear least-squares problems. \url{https://www.mathworks.com/help/optim/ug/lsqnonlin.html}.
\end{thebibliography}

\smallskip\small
\section*{人工智能使用报告}
人工智能工具协助了代码生成、数值检查和文稿起草。所有采用的内容须由队员核验并负责。

\textbf{队号和完整 AI 记录待补，尚非提交终稿。}竞赛规则要求附上所用 AI 工具的完整输入与输出，该记录须在提交前补全并复核。
\end{document}
'''

for lang, template in [('英文', template_en), ('中文', template_zh)]:
    path = ROOT / f'latex/AI/B题{lang}论文.tex'
    text=render(template).replace('[htbp]','[H]')
    text=text.replace(r'\begin{thebibliography}{9}',r'\begingroup\small\raggedright'+'\n'+r'\begin{thebibliography}{9}').replace(r'\end{document}',r'\endgroup'+'\n'+r'\end{document}')
    path.write_text(text, encoding='utf-8')

print('中英文单文件论文已由正式 JSON 重新生成。')
