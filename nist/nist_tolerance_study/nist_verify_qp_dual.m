function [accepted,details] = nist_verify_qp_dual(A,B,raw,originalCandidate,residualScale)
%NIST_VERIFY_QP_DUAL Verify the REPORTED Gurobi bound using its dual iterate.
% Does not replace it by another lower bound. The QP is min ||r||^2 with
% A*w<=r<=B*w and w in the simplex. For nonnegative multipliers alpha,beta,
% its dual function is min_j(A'*alpha-B'*beta)_j-||beta-alpha||^2/4.
accepted=false; details=struct('reason',"no_dual_iterate",'verified',false);
[m,~]=size(A);
if ~isfield(raw,'pi') || numel(raw.pi)<2*m || ...
        any(~isfinite(raw.pi(1:2*m))) || ~isfinite(originalCandidate)
    return
end
alpha=max(0,-double(raw.pi(1:m))); alpha=alpha(:);
beta=max(0,-double(raw.pi(m+1:2*m))); beta=beta(:);
left=A.'*alpha; right=B.'*beta;
gamma=(m+2)*eps/(1-(m+2)*eps);
error=gamma*(abs(A).'*alpha+abs(B).'*beta)+ ...
    8*eps(max(1,abs(left)+abs(right)));
linearMinimum=min(left-right-error);
difference=beta-alpha;
absoluteDifference=abs(difference)+4*eps(max(1,abs(alpha)+abs(beta)));
squared=norm(absoluteDifference,2)^2;
squaredUpper=squared+gamma*squared+8*eps(max(1,abs(squared)));
dualValue=linearMinimum-squaredUpper/4;
dualValue=dualValue-8*eps(max(1,abs(linearMinimum)+squaredUpper/4));
% Avoid overflow in residualScale^2 while comparing on the native QP scale.
nativeCandidate=(originalCandidate/residualScale)/residualScale;
accepted=isfinite(dualValue) && nativeCandidate<=dualValue;
details=struct('reason',"checked_dual_iterate",'verified',accepted, ...
    'reported_candidate_native_scale',nativeCandidate, ...
    'dual_witness_native_scale',dualValue,'residual_scale',residualScale);
end
