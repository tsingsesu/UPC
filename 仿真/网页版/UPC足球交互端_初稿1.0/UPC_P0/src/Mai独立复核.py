"""DOP853 independent audit of an explicit MATLAB run; never writes official data."""
from pathlib import Path
import argparse
import json
import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares


def audit(source):
    data = json.loads(source.read_text(encoding='utf-8-sig'))
    c = data['constants']
    r0 = np.asarray(c['r0'], dtype=float).reshape(3)
    gravity = np.asarray(c['g'], dtype=float).reshape(3)
    radius, mass = c['R'], c['m']
    k = c['rho'] * np.pi * radius**2 / (2 * mass)
    obs = np.asarray(data['observations']['positions'], dtype=float)
    times = np.asarray(data['observations']['times'], dtype=float)
    # This audit deliberately accepts only the existing still-air, constant-Cd baseline.
    cfg = data.get('config', {})
    if cfg:
        assert np.max(np.abs(cfg.get('w', [0, 0, 0]))) == 0
        assert cfg.get('bx', 0) == 0 and cfg.get('cdExponent', 0) == 0
        assert cfg.get('bDecayTau') in (None, float('inf'), 'infinite')
        assert np.allclose(cfg['lb'],[10,-20,-5,.05,-1.5,-1.5],rtol=0,atol=0)
        assert np.allclose(cfg['ub'],[50,20,25,.6,1.5,1.5],rtol=0,atol=0)
        assert np.all(cfg.get('freeMask',[True]*6))
    assert np.allclose(gravity,[0,0,-9.80665],rtol=0,atol=0)
    assert data['numericSettings']['maxFlightTime_s']==2.0

    def rhs(t, state, p):
        v = state[3:6]
        return np.r_[v, gravity-k*p[3]*np.linalg.norm(v)*v+np.cross([0,p[4],p[5]],v)]

    def positions(p):
        sol = solve_ivp(lambda t,s:rhs(t,s,p), (0, times[-1]), np.r_[r0,p[:3]],
                        method='DOP853',t_eval=times,rtol=2e-12,atol=2e-13)
        if not sol.success: raise RuntimeError(sol.message)
        return sol.y[:3].T

    def event(p, start=None, horizon=2):
        initial = r0.copy() if start is None else np.asarray(start)
        a0 = rhs(0,np.r_[initial,p[:3]],p)[5]
        if initial[2] < radius-1e-12: return {'status':'invalid_input'}
        if initial[2] <= radius+1e-12 and (p[2]<0 or (p[2]==0 and a0<=0)):
            return {'status':'ground','eventTime':0.0}
        def goal(t,s): return s[0]-11
        def ground(t,s): return s[2]-radius
        goal.terminal=True; goal.direction=1
        ground.terminal=True; ground.direction=-1
        sol=solve_ivp(lambda t,s:rhs(t,s,p),(0,horizon),np.r_[initial,p[:3]],
                      method='DOP853',rtol=2e-12,atol=2e-13,max_step=.0025,
                      events=[goal,ground])
        if not sol.success: raise RuntimeError(sol.message)
        if len(sol.t_events[0]):
            return {'status':'crossed','eventTime':float(sol.t_events[0][0]),
                    'xyz':sol.y_events[0][0,:3]}
        if len(sol.t_events[1]): return {'status':'ground','eventTime':float(sol.t_events[1][0])}
        return {'status':'timeout'}

    fit=least_squares(lambda p:(positions(p)-obs).ravel(),
                      [30.4,9.4,6.9,.3,0,-.4],
                      bounds=([10,-20,-5,.05,-1.5,-1.5],[50,20,25,.6,1.5,1.5]),
                      xtol=1e-13,ftol=1e-13,gtol=1e-13,x_scale='jac')
    assert fit.success
    point=event(fit.x)
    assert point['status']=='crossed'
    mp=data['prediction']
    diff=point['xyz'][1:]-[mp['yGoal'],mp['zGoal']]
    dt=point['eventTime']-mp['tCross']
    assert np.max(np.abs(diff))<1e-7 and abs(dt)<1e-7
    assert np.max(np.abs(np.asarray(data['fit']['theta'])-fit.x))<1e-6
    assert np.max(np.abs(fit.fun))<1e-8
    # Independent event regressions include the interval formerly skipped by both programs.
    down=event(np.array([30,0,-1,0,0,0]))
    early=event(np.array([30,0,-1,0,0,0]),[0,0,radius+.01])
    zero=event(np.zeros(6))
    assert down['status']=='ground' and down['eventTime']==0
    assert early['status']=='ground' and 0<early['eventTime']<.1
    early_exact=(-1+np.sqrt(1+2*9.80665*.01))/9.80665
    assert abs(early['eventTime']-early_exact)<1e-10
    assert zero['status']=='ground' and zero['eventTime']==0
    assert event(np.array([30,0,6,0,0,0]))['status']=='crossed'
    assert event(np.array([0,0,1,0,0,0]))['status']=='ground'
    def energy_rhs(t,s):
        return np.r_[rhs(t,s[:6],fit.x),mass*k*fit.x[3]*np.linalg.norm(s[3:6])**3]
    sol=solve_ivp(energy_rhs,(0,point['eventTime']),np.r_[r0,fit.x[:3],0],
                  method='DOP853',rtol=2e-12,atol=2e-13,max_step=.005)
    energy=mass*np.sum(sol.y[3:6]**2,axis=0)/2-mass*gravity[2]*sol.y[2]
    energy_error=float(np.max(np.abs(energy+sol.y[6]-energy[0])))
    assert energy_error<1e-8
    return {'passed':True,'scipyVersion':scipy.__version__,'numpyVersion':np.__version__,
            'method':'independent DOP853 from t=0 + bounded least_squares',
            'theta':fit.x.tolist(),'goalTime_s':point['eventTime'],
            'goalXYZ_m':point['xyz'].tolist(),'matlabDifferenceYZ_m':diff.tolist(),
            'matlabDifferenceTime_s':dt,'maxObservationResidual_m':float(np.max(np.abs(fit.fun))),
            'energyBalanceMaxResidual_J':energy_error,
            'eventTests':{'downwardAtGround':down,'earlyContact':early,'stationaryAtGround':zero,
                          'upwardLaunch':'passed','noForwardCrossing':'passed'},
            'scope':'仅数值验证；静风、常Cd、bx=0、无b衰减。不是实验确认。'}


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--matlab',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    report=audit(args.matlab)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    # Exclusive create protects any earlier audit evidence.
    with args.output.open('x',encoding='utf-8') as f:
        json.dump(report,f,ensure_ascii=False,indent=2,allow_nan=False)
    print(json.dumps(report,ensure_ascii=True,indent=2))
