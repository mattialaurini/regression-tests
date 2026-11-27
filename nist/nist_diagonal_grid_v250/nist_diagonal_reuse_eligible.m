function yes=nist_diagonal_reuse_eligible(raw,originalOptions,expectedOptions)
%NIST_DIAGONAL_REUSE_ELIGIBLE Budget changes have an explicit, narrow policy.
yes=false;
if strlength(raw.error_message)>0 || raw.epsilon_absolute~=1e-10, return, end
oldLimit=raw.time_limit; target=expectedOptions.TimeLimit;
if originalOptions.TimeLimit~=oldLimit, return, end
old=originalOptions; expected=expectedOptions;
old.TimeLimit=target;
if ~isequaln(old,expected), return, end
if oldLimit==target
    yes=true;
elseif raw.epsilon_relative==0.01 && oldLimit==600 && target==1000 && logical(raw.certified)
    % Explicitly requested reuse. The 600-s budget remains in the row;
    % this is not claimed to be an identical 1000-s execution, since
    % the total local-search allowance depends on the global budget.
    yes=true;
end
end
