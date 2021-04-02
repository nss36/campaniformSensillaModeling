function ddt = centeredDiff(t,x)
    tLen = length(t);
    xLen = length(x);
    
    if tLen ~= xLen
        error('t and x must have the same length.')
    end
    
    if size(t,1) ~= size(x,1)
        error('t and x must both be columns, or must both be rows.')
    end
    
    ddt = zeros(size(x));
    
    ddt(2:end-1) = (x(1:end-2) - x(3:end))./(t(1:end-2) - t(3:end));
    ddt(1) = (x(2) - x(1))./(t(2) - t(1));
    ddt(end) = (x(end) - x(end-1))./(t(end) - t(end-1));
    
end