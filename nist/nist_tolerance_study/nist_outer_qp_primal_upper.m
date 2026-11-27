function [upper,w] = nist_outer_qp_primal_upper(A,B,result,options)
%NIST_OUTER_QP_PRIMAL_UPPER Feasible primal value for screening, never pruning.
[m,r]=size(A); upper=Inf; w=[];
if ~result.hasSolution || numel(result.x)~=r+m || ...
        any(~isfinite(result.x(1:r))) || any(A>B,'all')
    return
end
w=max(0,double(result.x(1:r))); total=sum(w);
if ~isfinite(total) || total<=0, w=[]; return, end
w=w/total;
quantum=2^48; units=floor(w*quantum);
[~,largest]=max(units);
units(largest)=units(largest)+(quantum-sum(units));
if any(units<0), w=[]; return, end
w=units/quantum;
lo=A*w; hi=B*w;
gamma=(r+2)*eps/(1-(r+2)*eps);
loUpper=lo+gamma*(abs(A)*w)+32*eps(max(1,abs(lo)));
hiLower=hi-gamma*(abs(B)*w)-32*eps(max(1,abs(hi)));
distanceUpper=max(0,max(loUpper,-hiLower));
value=norm(distanceUpper,2)^2;
if isfinite(value)
    upper=value+double(options.BoundSafety)*(1+abs(value));
end
end
