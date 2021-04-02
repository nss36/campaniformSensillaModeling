function [y,x,dudt,t] = simulateMinusLowpassPL(u,tau,tvec,a,b,c,d,x0)
    


%     fprintf('tau = %3.6f, a = %3.3f, b = %3.3f, c = %3.3f, d = %3.3f...',tau,a,b,c,d)

    %u: input signal
    %y: output signal
    %x: "hidden" low pass signal
    %tau: time constant of low pass filter
    dudt = centeredDiff(tvec,u);
    
    numSteps = length(tvec);
    
    if isequal(size(tvec),[1,numSteps])
        %good, do nothing
    elseif isequal(size(tvec),[numSteps,1])
        tvec = tvec';
    else
        error('t must be a vector with 1 row and numSteps columns.')
    end
    
    if isequal(size(u),[1,numSteps])
        u = u';
    elseif isequal(size(u),[numSteps,1])
        %good, do nothing
    else
        error('u must be a vector with numSteps rows and 1 column.')
    end
    
    f = @(t,y) 1/tau*sign(interp1(tvec,u,t,'pchip')-y).*(abs(interp1(tvec,u,t,'pchip')-y)).^b;
    options = odeset('RelTol',1e-7,'AbsTol',1e-10);
    [t,x] = ode15s(f,tvec,x0,options); 
    
    if length(tvec) == 2
        %in this case, ode assumes these two values are the start and end
        %times, and will return many values, but we only want two.
        t(2:end-1) = [];
        x(2:end-1) = [];
    end
    
    y = max(0,a*(u - x) + c*u + d);
    
%     fprintf('done.\n')
end