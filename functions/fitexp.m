function [estimates, model, SSE] = fitexp(xdata, ydata)
% fits: A * exp(-xdata/tau) - ydata
model = @expfun;
estimates = fminsearch(model, [.05, 400, 1]);
    function [SSE, FittedCurve] = expfun(params)
        A = params(1);
        tau = params(2);
        C = params(3);
        FittedCurve = A .* exp(-xdata / tau) + C;
        ErrorVector = FittedCurve - ydata;
        SSE = sum(ErrorVector .^ 2);
    end
[SSE, FittedCurve] = model(estimates);
%%confidence interval:
%MSE = SSE ./ (length(xdata)-3); % df = N - number of params
%S = (X*X')^-1*MSE; %X is 
%delta = tinv(0.682,length(xdata)-3).*sqrt(S(2,2)); %(1 std dev.)
%[tau - delta, tau + delta]

end
