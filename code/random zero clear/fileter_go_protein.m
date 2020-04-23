function [ X1,golist,prolist,nogo ] = fileter_go_protein(X)
    X1 = X;
    sumrow=sum(X1,2);
    [nogo,~]=find(sumrow<10 | sumrow>200);
    [golist,~]=find(sumrow>=10 & sumrow<=200);
    X1(nogo,:)=0;
    sumcol=sum(X1);
    [~,prolist]=find(sumcol~=0);  
end

