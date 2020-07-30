

 function [xZ,stbl]= ZeroGPR(x,y,err,count_max,Npoints)
     
x=x(:); 
yfit=y;xfit=x;
bad=isnan(yfit); yfit(bad)=[]; xfit(bad)=[]; %%find ind below see's NaN's as things that change sign therefore flagging them as a zero
N=length(x); 

FindInd = @(x) find(diff(sign(x))); 
indZ=FindInd(yfit);

nZ=length(indZ);  yZ=yfit(indZ); xZ=xfit(indZ); 

stbl=ones(nZ,1); %if stable stbl=1, unstable=0
for k=1:nZ
    count=1; 
     ind=indZ(k);
     if ind >= 3 && ind <= N-2
         xpoint=[x(ind-2), x(ind-1),x(ind),x(ind+1),x(ind+2)];
         ypoint=[y(ind-2), y(ind-1),y(ind),y(ind+1),y(ind+2)];
         p=polyfit(xpoint,ypoint,1);
         stbl(k)=sign(p(1));
     else
         stbl(k)=NaN; 
     %elseif ind >= 2 && ind <= Npoints-1
     %    xpoint=[ x(ind-1),x(ind),x(ind+1)];
     %    ypoint=[ y(ind-1),y(ind),y(ind+1)];
     %    p=polyfit(xpoint,ypoint,1);
     %    stbl(k)=sign(p(1));
    while abs(yZ(k)) > err && count<count_max
        count=count+1; 
        if ind ~= 1 && ind ~= Npoints
            ind_range= ind-1:ind+1; 
        else 
            break
        end
        xdata=xfit(ind_range); ydata=yfit(ind_range);
        model2=fitrgp(xdata(:),ydata(:),'FitMethod','sd','ComputationMethod','v','PredictMethod','sr','InitialStepSize','auto');
        xfit=linspace(min(xdata),max(xdata),Npoints)'; 
        yfit=predict(model2,xfit);
        ind2=FindInd(yfit); 
        if isempty(ind2)~=1 && length(ind2)==1 
        yZ(k)=yfit(ind2); xZ(k)=xfit(ind2); 
        end
    end
    xfit=linspace(min(x) - std(x),max(x)+std(x),Npoints)'; 
    yfit=y;
    
 
end

%indBad= find (xZ>xM || xZ<xm); 


 end 
 

