function [precision,recall,fmeasure,auc,aup] = newmetric(ground_truth,predict )
%初始点为（1.0, 1.0）
%计算出ground_truth中正样本的数目pos_num和负样本的数目neg_num

pos_num = sum(ground_truth==1);
neg_num = sum(ground_truth==0);

m=size(ground_truth,1);
[pre,Index]=sort(predict);
ground_truth=ground_truth(Index);
x=zeros(m+1,1);
y=zeros(m+1,1);
p=zeros(m+1,1);
f=zeros(m+1,1);
auc=0;
aup=0;
x(1)=1;y(1)=1;p(1)=1;
TP=sum(ground_truth(2:m)==1);
FP=sum(ground_truth(2:m)==0);


for i=2:m
    if i~=2
        if ground_truth(i-1)==1
            TP=TP-1;
        else
            FP=FP-1;
        end
    end       
if neg_num==0
    x(i)=0;
else
    x(i)=FP/neg_num;
end
if pos_num==0
y(i)=0;
else
y(i)=TP/pos_num;
end
if TP+FP==0
p(i)=0;
else
p(i)=TP/(TP+FP);
end
auc=auc+(y(i)+y(i-1))*(x(i-1)-x(i))/2;
aup=aup+(p(i)+p(i-1))*(y(i-1)-y(i))/2;
end

%x(m+1)=0;y(m+1)=0;
auc=auc+y(m)*x(m)/2;
aup=aup+p(m)*y(m)/2;
for i=1:m+1
    if p(i)+y(i)==0
	f(i)=0;
    else
        f(i)=2*p(i)*y(i)/(p(i)+y(i));
    end
end
f(1)=0;
p(1)=0;
y(1)=0;
[fmeasure,index]=max(f);
precision=p(index);
recall=y(index);
end

