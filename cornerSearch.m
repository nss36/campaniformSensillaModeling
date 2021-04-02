function [xf,x,i,converged] = cornerSearch(f,x0,p0,lowerBound,upperBound,maxIter,xTol,fTol,findHighPoint)
    
    x = NaN(maxIter,1);
    x(1) = x0;
    p = abs(p0);
    i = 2;
    xf = NaN;
    converged = false;
    
    if abs(f(x0)) > fTol
        warning('Function must return 0 for initial guess in cornerSearch.')
    elseif abs(f(x0)) <= fTol && abs(f(x0 + xTol)) > fTol && findHighPoint
        converged = true;
        xf = x0;
        i = 0;
        %initial point is already the high corner point.
    elseif abs(f(x0)) <= fTol && abs(f(x0 - xTol)) > fTol && ~findHighPoint
        converged = true;
        xf = x0;
        i = 0;
        %initial point is already the low corner point.
    else
        while ~converged && i <= maxIter
            newPointFound = false;
            j = 1;
            while ~newPointFound && j <= maxIter
                %Compute new guess for the corner.
                if findHighPoint
                    xTry = x(i-1) + p;
                else
                    xTry = x(i-1) - p;
                end
                fprintf('xTry = %6.6f\n',xTry)
                
                %Compute f at the new guess. If the new guess it outside
                %the bounds, save time by returning a nonzero value for f
                %without calling f.
                if xTry > upperBound
                    fTry = 1;
                elseif xTry < lowerBound
                    fTry = 1;
                else
                    fTry = f(xTry);
                end
                
                %If the next guess returns f = 0, then we have a new point
                %and can progress in the algorithm.
                if abs(fTry) <= fTol
                    newPointFound = true;
                    x(i) = xTry;
                else
                    %This guess went too far, so we must decrease p, and we
                    %have not found a new point.
                    p = p/2;
                    if p < xTol
                        %We can stop now, because there is no use finding a
                        %more precise answer.
                        newPointFound = true;
                    end
                end
                j = j + 1;
            end
            %If the step is smaller than the tolerance, then we say we have
            %converged and return the final guess.
            if p < xTol
                converged = true;
                xf = x(i-1);
            end
            i = i + 1;
        end
    end
end