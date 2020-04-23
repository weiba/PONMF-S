function [ p,r,f,auc,aup ] = metric_per( true_matrix,result_matrix)
[~,col]=size(true_matrix);
p_list=zeros(1,col);
r_list=zeros(1,col);
f_list=zeros(1,col);
auc_list=zeros(1,col);
aup_list=zeros(1,col);
for i=1:col
    [precision,recall,fmeasure,auc,aup]=newmetric(true_matrix(:,i),result_matrix(:,i));
    p_list(i)=precision;
    r_list(i)=recall;
    f_list(i)=fmeasure;
    auc_list(i)=auc;
    aup_list(i)=aup;
end
p=mean(p_list);
r=mean(r_list);
f=mean(f_list);
auc=mean(auc_list);
aup=mean(aup_list);
end


