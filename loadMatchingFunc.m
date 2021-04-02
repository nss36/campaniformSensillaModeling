function ceq = loadMatchingFunc(x,i,tau,a,b,c,d,tSamp,uSamp,yDes,x0)

    u = [uSamp(i-1);x];
    
    t = tSamp(i-1:i);
    
    try y = simulateMinusLowpassPL(u,tau,t,a,b,c,d,x0);
    catch
        keyboard
    end
    ceq = yDes(i) - y(end);
    
% %To animate the inversion process, uncomment the code below.
%     figure(100)
%     clf
%     plot(tSamp(1:i),yDes(1:i))
%     hold on
%     plot(tSamp(i),yDes(i,:),'o')
%     plot(t,y,'m:','linewidth',2)
%     title(sprintf('x = %3.3e',x))

end