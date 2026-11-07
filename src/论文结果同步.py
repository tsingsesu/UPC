"""将MATLAB正式结果同步到英文论文；不包含绘图。"""
from pathlib import Path
import json
import numpy as np
from scipy.stats import chi2
ROOT=Path(__file__).resolve().parents[1]
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
