function ceq = loadMatchingFuncMulti(x,i,tau,a,b,c,d,tSamp,uSamp,yDes,x0,toInvert,toPlot)

%     u = [uSamp(i-1);x];
%     
%     t = tSamp(i-1:i);
%     
%     y = simulateMinusLowpassPL(u,tau,t,a,b,c,d,x0);
%     ceq = yDes(i) - y(end);
    
    [~,n] = size(yDes); %num time points, num sensors
    y = zeros(2,n);
    
    %There's only one x, i, tSamp cols.
    %There are n tau, a, b, c, d, uSamp cols.
    for j=1:n
        if toInvert(j)
            u = -[uSamp(i-1,j);x];
            x0act = x0(j);
        else
            u = [uSamp(i-1,j);x];
            x0act = x0(j);
        end
        
        t = tSamp(i-1:i);

        y(:,j) = simulateMinusLowpassPL(u,tau(j),t,a(j),b(j),c(j),d(j),x0act);
    end
    
%     yDes(i,:) - y(end,:)
%     ceq = sum(yDes(i,:) - y(end,:));
%     ceq = sum(abs(yDes(i,:) - y(end,:)));
    ceq = sum((yDes(i,:) - y(end,:)).^2);
%     ceq = max(0,sum(-1 + abs(yDes(i,:) - y(end,:))));

    if toPlot
        figure(100)
        clf
        plot(tSamp,yDes)
        hold on
        plot(tSamp(i),yDes(i,:),'o')
        plot(t,y,'m:','linewidth',2)
        title(sprintf('x = %3.3e, f = %3.3e',x,ceq))
        
%         y
%     
%         
%         keyboard
    end

end