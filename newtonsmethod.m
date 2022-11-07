% Newton's method

% This function implements Newton's method for rootfinding, after Atkinson 
% and Han (2004)
%
% Michael McCarthy, November 2022 (michael.mccarthy@wsl.ch)

function root = newtonsmethod(fun,x0,tolX,maxIter,range,maxStepX)

% Define starting error and iteration number
errX = maxStepX;
nIter = 0;

% While the error on x is more than the specified tolerance and the maximum 
% number of iterations has not been reached, continue to make new estimates
% of x
while abs(errX) > tolX && nIter < maxIter
    
    % Calculate f(x) and its derivative f(x)'
    fx = fun(x0);
    dfx = (fun(x0+range)-fun(x0-range))/(2*range);
    
    % Stop if f(x)' is zero
    if dfx == 0
        return
    end
    
    % Make new estimate of x
    x1 = x0-fx/dfx;
    errX = x1-x0;
    
    % Do not allow steps of more than the maximum specified step size
    if errX > maxStepX
        x1 = x0+maxStepX;
    end
    if errX < -maxStepX
        x1 = x0-maxStepX;    
    end
    
    % Update x and iteration number
    x0 = x1;
    nIter = nIter+1;
end

% When the error on x is less than the specified tolerance, the root is the 
% current iteration of x. If the maximum number of iterations was reached,
% the root is the mean of the last two iterations of x
if nIter > maxIter
    root = (x0+x1)/2;
else
    root = x1;
end

end

