"""Python辅助复核：独立使用DOP853求解，验证MATLAB主模型；不绘图。"""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
import scipy

ROOT = Path(__file__).resolve().parents[2]  # 仿真/Matlab版/ → 项目根
OBS = np.array([[1.5, .45, .44], [2.96, .86, .74]])
R0 = np.array([0., 0., .11])
K = 1.225 * np.pi * .11**2 / (2*.43)
G = 9.80665

def rhs(t, s, p):
    v = s[3:]
    return np.r_[v, np.array([0., 0., -G]) - K*p[3]*np.linalg.norm(v)*v
                 + np.cross([0., p[4], p[5]], v)]

def observations(p):
    sol = solve_ivp(lambda t,s:rhs(t,s,p), (0., .1), np.r_[R0,p[:3]],
                    method='DOP853', t_eval=[.05,.1], rtol=2e-12, atol=2e-13)
    if not sol.success: raise RuntimeError(sol.message)
    return sol.y[:3].T

def crossing(p):
    def goal(t,s): return s[0]-11.
    def ground(t,s): return s[2]-.11
    goal.terminal=True; goal.direction=1
    ground.terminal=True; ground.direction=-1
    init=solve_ivp(lambda t,s:rhs(t,s,p),(0,.1),np.r_[R0,p[:3]],
                   method='DOP853',rtol=2e-12,atol=2e-13).y[:,-1]
    sol=solve_ivp(lambda t,s:rhs(t,s,p),(.1,2.),init,method='DOP853',
                  events=[goal,ground],rtol=2e-12,atol=2e-13,max_step=.005)
    if not sol.success: raise RuntimeError(sol.message)
    if len(sol.t_events[0])!=1: raise RuntimeError('未在触地前穿过球门平面')
    return np.r_[sol.t_events[0][0],sol.y_events[0][0,:3]]

def main():
    fit=least_squares(lambda p:(observations(p)-OBS).ravel(),
                      [30.4,9.4,6.9,.3,0.,-.4],
                      bounds=([10,-20,-5,.05,-1.5,-1.5],[50,20,25,.6,1.5,1.5]),
                      xtol=1e-13,ftol=1e-13,gtol=1e-13, x_scale='jac')
    if not fit.success: raise RuntimeError(fit.message)
    point=crossing(fit.x)
    deriv=[]
    for j in range(6):
        d=1e-5*max(1.,abs(fit.x[j])); p=fit.x.copy(); q=p.copy()
        p[j]+=d; q[j]-=d
        deriv.append((crossing(p)[2:]-crossing(q)[2:])/(2*d))
    jac=np.array(deriv).T
    cov=jac@np.linalg.inv(fit.jac.T@fit.jac)@jac.T*.01**2
    v0=(4*OBS[0]-OBS[1]-3*R0)/.1
    acc=(OBS[1]-2*OBS[0]+R0)/.05**2
    roots=np.roots([acc[0]/2,v0[0],-11.])
    t=min(x.real for x in roots if abs(x.imag)<1e-10 and x.real>0)
    constant=R0+v0*t+.5*acc*t*t
    assert np.max(np.abs(observations(fit.x)-OBS)) < 1e-8
    report={'scipy_version':scipy.__version__,'method':'DOP853 + bounded least_squares',
            'theta':fit.x.tolist(),'observation_max_residual_m':float(np.max(np.abs(fit.fun))),
            'goal_time_and_xyz':point.tolist(),'sigma_m':.01,
            'goal_covariance_m2':cov.tolist(),
            'marginal_95_halfwidth_m':(1.9599639845*np.sqrt(np.diag(cov))).tolist(),
            'constant_acceleration':{'v0':v0.tolist(),'a':acc.tolist(),'t':t,'xyz':constant.tolist()},
            'scope':'独立程序数值复核，不是实验确认；噪声标准差为假设。'}
    def extended(t,s):
        return np.r_[rhs(t,s[:6],fit.x),.43*K*fit.x[3]*np.linalg.norm(s[3:6])**3]
    sol=solve_ivp(extended,(0,point[0]),np.r_[R0,fit.x[:3],0.],
                  method='DOP853',rtol=2e-12,atol=2e-13,max_step=.005)
    energies=.43*np.sum(sol.y[3:6]**2,axis=0)/2+.43*G*sol.y[2]
    report['energy_balance_max_residual_J']=float(np.max(np.abs(energies+sol.y[6]-energies[0])))
    matlab_path=ROOT/'data'/'派生数据'/'罚球结果.json'
    if matlab_path.exists():
        main_data=json.loads(matlab_path.read_text(encoding='utf-8'))
        matlab_point=main_data['prediction']
        report['matlab_difference_yz_m']=(point[2:]-np.array([matlab_point['yGoal'],matlab_point['zGoal']])).tolist()
        assert np.max(np.abs(report['matlab_difference_yz_m']))<1e-7
        assert np.max(np.abs(np.array(main_data['fit']['theta'])-fit.x))<1e-6
    report['trajectory_family']=verify_family(main_data['trajectoryFamily'])
    target=ROOT/'data'/'派生数据'/'独立复核结果.json'
    target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(report,ensure_ascii=False,indent=2))

def verify_family(family):
    """逐候选独立积分；另对两端代表独立重拟合，无绘图。"""
    max_obs=0.; max_goal=0.; max_refit=0.; disagreements=0
    details=[]
    for index,a in enumerate(family['members'],1):
        if a['exitflag']<=0:
            continue
        th=np.asarray(a['theta']); wind=np.array([0.,a['windY'],0.])
        k=a['density']*np.pi*.11**2/(2*.43)
        def force(t,s,theta=th):
            u=s[3:]-wind
            return np.r_[s[3:], [0.,0.,-G]] + np.r_[np.zeros(3),
                -k*theta[3]*np.linalg.norm(u)*u+np.cross([0.,theta[4],theta[5]],u)]
        def plane(t,s): return s[0]-11.
        def ground(t,s): return s[2]-.11 if t>1e-9 else 1e-8
        plane.terminal=True;plane.direction=1
        ground.terminal=True;ground.direction=-1
        sol=solve_ivp(force,(0,2),np.r_[R0,th[:3]],method='DOP853',
                      events=[plane,ground],dense_output=True,rtol=2e-12,atol=2e-13,max_step=.004)
        assert sol.success
        observed=sol.sol([.05,.1])[:3].T
        residual=observed-OBS
        error=float(np.max(np.abs(residual)))
        max_obs=max(max_obs,abs(error-a['maxResidual']))
        times=np.linspace(0.,sol.t[-1],max(2,int(np.ceil(sol.t[-1]/.001))+1))
        vel=sol.sol(times)[3:].T-wind
        cl=np.linalg.norm(np.cross([0.,th[4],th[5]],vel),axis=1)/(k*np.maximum(np.sum(vel**2,axis=1),1e-16))
        limits=family['physicalLimits']; speed=np.linalg.norm(th[:3])
        physical=(limits['speedMin']<=speed<=limits['speedMax'] and th[0]>0 and th[2]>=0
                  and np.max(cl)<=limits['maxTransverseCoefficient'] and sol.t[-1]>=.1)
        assert physical==a['physicalAdmissible']
        assert abs(np.max(cl)-a['maxTransverseCoefficient'])<1e-8
        accepted=error<=family['defaultTolerance'] and physical
        disagreements+=int(accepted != a['accepted'])
        expected='plane' if len(sol.t_events[0]) else ('ground' if len(sol.t_events[1]) else 'timeout')
        assert expected==a['status']
        max_goal=max(max_goal,float(np.linalg.norm(sol.y[:3,-1]-np.array(a['endpoint']))))
        if index in family['exampleIndices']:
            def residual_fun(q):
                par=np.r_[q[:3],a['Cd'],q[3:]]
                trial=solve_ivp(lambda t,y:force(t,y,par),(0,.1),np.r_[R0,q[:3]],
                                method='DOP853',t_eval=[.05,.1],rtol=2e-12,atol=2e-13)
                assert trial.success
                return 1000*(trial.y[:3].T-OBS).ravel()
            fit=least_squares(residual_fun,[30.4,9.4,6.9,-.06,-.36],
                              bounds=([10,-20,0,-1.5,-1.5],[40,20,25,1.5,1.5]),
                              xtol=1e-12,ftol=1e-12,gtol=1e-10,x_scale='jac')
            assert fit.success
            delta=float(np.max(np.abs(fit.x-th[[0,1,2,4,5]])))
            max_refit=max(max_refit,delta)
            details.append({'index':index,'independent_q':fit.x.tolist(),'max_parameter_difference':delta})
    assert max_obs<1e-8 and max_goal<1e-7 and disagreements==0 and max_refit<1e-5
    return {'n_checked':len(family['members']),'max_observation_residual_difference_m':max_obs,
            'max_goal_difference_m':max_goal,'acceptance_disagreements':disagreements,
            'independent_refits':details,'max_refit_parameter_difference':max_refit,
            'scope':'Numerical consistency only; no experimental confirmation.'}

if __name__=='__main__': main()
