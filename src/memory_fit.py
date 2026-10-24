"""One-parameter binomial parity-decay fit, with explicit lack-of-fit diagnostics."""
import math
import numpy as np


def decay_probability(epsilon,rounds):
    r=np.asarray(rounds,dtype=float)
    return -np.expm1(r*np.log1p(-2*epsilon))/2


def _chi_square_survival_3(value):
    x=max(float(value),0.0)
    return math.erfc(math.sqrt(x/2))+math.sqrt(2*x/math.pi)*math.exp(-x/2)


def fit_epsilon(points,*,metric='fails',minimum_rounds=2):
    rows=sorted([p for p in points if p['rounds']>=minimum_rounds],key=lambda p:p['rounds'])
    if len(rows)<2:raise ValueError('At least two non-R1 points required')
    rounds=np.array([p['rounds'] for p in rows],float)
    n=np.array([p['shots'] for p in rows],float);k=np.array([p[metric] for p in rows],float)
    if np.any(k<0) or np.any(k>n) or np.any(n<=0):raise ValueError('Invalid binomial counts')
    def score(e):
        q=decay_probability(e,rounds);dq=rounds*(1-2*e)**(rounds-1)
        return float(np.sum((k-n*q)*dq/(q*(1-q))))
    lo,hi=1e-14,.499999
    if score(lo)<=0:raise ValueError('No failures: increase shots before fitting')
    for _ in range(90):
        mid=(lo+hi)/2
        if score(mid)>0:lo=mid
        else:hi=mid
    epsilon=(lo+hi)/2;q=decay_probability(epsilon,rounds)
    dq=rounds*(1-2*epsilon)**(rounds-1)
    information=float(np.sum(n*dq*dq/(q*(1-q))))
    d2q=-2*rounds*(rounds-1)*(1-2*epsilon)**(rounds-2)
    residual=k-n*q;variance=q*(1-q)
    observed_information=float(np.sum(n*dq*dq/variance-residual*d2q/variance
                               +residual*dq*dq*(1-2*q)/(variance*variance)))
    if observed_information<=0:raise ValueError('Fit does not have positive curvature')
    influence=dq/(variance*observed_information)
    # Delta-method sampling variance of this finite-window MLE. The endpoints
    # at different R are independent, even when the fixed curve lacks fit.
    stderr=math.sqrt(float(np.sum(influence*influence*n*(k/n)*(1-k/n))))
    observed=k/n;pearson=float(np.sum((k-n*q)**2/(n*q*(1-q))))
    df=len(rows)-1
    # Four requested fit points give df=3; two-point tail fits give df=1.
    pvalue=(_chi_square_survival_3(pearson) if df==3 else math.exp(-pearson/2) if df==2
            else math.erfc(math.sqrt(pearson/2)) if df==1 else None)
    with np.errstate(divide='ignore',invalid='ignore'):
        term=np.where(k>0,k*np.log(np.maximum(k,1)/(n*q)),0)
        term+=np.where(k<n,(n-k)*np.log(np.maximum(n-k,1)/(n*(1-q))),0)
    deviance=float(2*term.sum())
    return {'epsilon':epsilon,'stderr':stderr,'relative_stderr':stderr/epsilon,
            'ci95_normal':[max(0,epsilon-1.96*stderr),epsilon+1.96*stderr],
            'metric':metric,'fit_rounds':rounds.astype(int).tolist(),
            'method':'Binomial maximum likelihood; sampling stderr from observed likelihood curvature and independent endpoint-binomial variances.',
            'fisher_stderr_secondary':1/math.sqrt(information),
            'observed_information':observed_information,
            'model':'p_L(R)=(1-(1-2*epsilon)^R)/2',
            'pearson_chi2':pearson,'degrees_of_freedom':df,'reduced_chi2':pearson/df,
            'goodness_of_fit_p_value':pvalue,'binomial_deviance':deviance,
            'sampling_stderr_warning':'Sampling error of the finite-window fitted parameter; lack of fit limits stationary-rate interpretation and is not cured by more shots.',
            'dispersion_scaled_stderr_diagnostic':stderr*math.sqrt(max(1,pearson/df)),
            'points':[{'rounds':int(r),'shots':int(nn),'fails':int(kk),'observed_p_L':float(y),
                       'fitted_p_L':float(qq),'residual_in_binomial_stderr':float((y-qq)/math.sqrt(qq*(1-qq)/nn)),
                       'epsilon_influence_per_failure_count':float(weight)}
                      for r,nn,kk,y,qq,weight in zip(rounds,n,k,observed,q,influence)]}


def combine_basis_fits(points,fit_X,fit_Z):
    """User-defined epsilon_tot, retaining shared-shot X/Z covariance.

    Cov(kX,kZ) = N*(pXZ-pX*pZ) at each R. R streams are independent.
    No union decay curve is fitted and no independent-marginal assumption is
    made for uncertainty propagation.
    """
    if fit_X['fit_rounds']!=fit_Z['fit_rounds']:raise ValueError('Basis windows differ')
    by_R={p['rounds']:p for p in points};covariance=0.
    for px,pz in zip(fit_X['points'],fit_Z['points']):
        p=by_R[px['rounds']];n=p['shots']
        joint=p['logical_failure_counts']['both_XZ']+p.get('aborted',0)
        covariance+=(joint-p['fails_X']*p['fails_Z']/n)*px['epsilon_influence_per_failure_count']*pz['epsilon_influence_per_failure_count']
    x,z=fit_X['epsilon'],fit_Z['epsilon'];sx,sz=fit_X['stderr'],fit_Z['stderr']
    correlation=covariance/(sx*sz)
    if abs(correlation)>1+1e-9:raise ValueError('Invalid joint endpoint covariance')
    variance=(1-z)**2*sx*sx+(1-x)**2*sz*sz+2*(1-x)*(1-z)*covariance
    epsilon=1-(1-x)*(1-z);stderr=math.sqrt(max(0,variance))
    dx,dz=fit_X['dispersion_scaled_stderr_diagnostic'],fit_Z['dispersion_scaled_stderr_diagnostic']
    scaled_var=(1-z)**2*dx*dx+(1-x)**2*dz*dz+2*(1-x)*(1-z)*correlation*dx*dz
    return {'epsilon':epsilon,'epsilon_tot':epsilon,'epsilon_X':x,'epsilon_Z':z,
            'stderr':stderr,'relative_stderr':stderr/epsilon,
            'ci95_normal':[max(0,epsilon-1.96*stderr),epsilon+1.96*stderr],
            'fit_rounds':fit_X['fit_rounds'],'method':'epsilon_tot=1-(1-epsilon_X)(1-epsilon_Z), joint delta-method propagation including shared-trajectory X/Z covariance.',
            'covariance_epsilon_X_Z':covariance,'correlation_epsilon_X_Z':correlation,
            'stderr_if_covariance_ignored_secondary':math.hypot((1-z)*sx,(1-x)*sz),
            'dispersion_scaled_stderr_diagnostic':math.sqrt(max(0,scaled_var)),
            'union_is_fitted':False,
            'basis_goodness_of_fit':{label:{k:fit[k] for k in ('pearson_chi2','degrees_of_freedom','goodness_of_fit_p_value')}
                                     for label,fit in [('epsilon_X',fit_X),('epsilon_Z',fit_Z)]}}


def R2_trend_diagnostic(points,tail_fit,metric):
    """Hold R=2 out; compare it with the fit to R>=4."""
    p=next(p for p in points if p['rounds']==2);observed=p[metric]/p['shots']
    stderr=math.sqrt(max(observed,1/p['shots'])*(1-observed)/p['shots'])
    prediction=float(decay_probability(tail_fit['epsilon'],2))
    trend_stderr=2*(1-2*tail_fit['epsilon'])*tail_fit['stderr']
    combined=math.hypot(stderr,trend_stderr);difference=observed-prediction
    return {'R2_observed':observed,'Rge4_prediction_at_R2':prediction,'observed_minus_trend':difference,
            'R2_MC_stderr':stderr,'trend_prediction_stderr':trend_stderr,'combined_stderr':combined,
            'sigma_combined':difference/combined,'deviation_in_R2_errorbars':difference/stderr,
            'outside_R2_one_sigma_MC_errorbar':abs(difference)>stderr,
            'criterion':'One error bar means the reported one-sigma R2 MC stderr; the independent trend-fit uncertainty is also propagated separately.'}


def compare_fits(a,b):
    """a minus b, with independent design streams; never compare two bases here."""
    delta=a['epsilon']-b['epsilon'];error=math.hypot(a['stderr'],b['stderr'])
    scaled=math.hypot(a['dispersion_scaled_stderr_diagnostic'],b['dispersion_scaled_stderr_diagnostic'])
    return {'delta_epsilon':delta,'combined_stderr':error,'sigma':delta/error,
            'combined_dispersion_scaled_stderr_diagnostic':scaled,
            'sigma_dispersion_diagnostic':delta/scaled,
            'dispersion_diagnostic_note':'Heuristic Pearson variance inflation, not a substitute for a correctly specified memory model.',
            'relative_difference':delta/b['epsilon'],'significant_at_3sigma':abs(delta)>=3*error,
            'interpretation':'Sampling comparison of finite-window fitted metrics; inspect basis fit quality before a stationary-rate claim.'}
