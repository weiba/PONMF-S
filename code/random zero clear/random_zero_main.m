%随机清零
clear();
n=20;
per=5;
n_parameter=0.225;
s_parameter=0.001;
u_parameter=0.4;
v_parameter=0.4;
X=load('./newyeastdata/C/newCgp.txt');
Ln=load('./newyeastdata/P/newPPI.txt');
Ls=load('./newyeastdata/C/newCgogo.txt');
Ln1 = cosin(Ln);                           %计算Ln的cosin相似性，输入的只能是0，1矩阵，返回0-1的矩阵。返回矩阵中可能存在全0行
ln = 0.5*Ln+0.5*Ln1;                   %原来的0，1矩阵与计算后的cosin矩阵相加
Ln = Ln*n_parameter;
Ls = Ls*s_parameter;
[ Xt,golist,prolist,~ ] = fileter_go_protein(X);
[Index_zeroRow,Index_zeroCol] = find(X(:,:)==0);
[Index_PositiveRow,Index_PositiveCol]=find(Xt(:,:)==1);%排除不考虑的GO
totalassociation=length(Index_PositiveRow);
totalper=fix(totalassociation/per);
zero_length = length(Index_zeroRow);

[row,col]=size(X);
PPrecision_list1=zeros(n,1);
RRecall_list1=zeros(n,1);
FFmeasure_list1=zeros(n,1);
AAuc_list1=zeros(n,1);
AAup_list1=zeros(n,1);

PPrecision_list2=zeros(n,1);
RRecall_list2=zeros(n,1);
FFmeasure_list2=zeros(n,1);
AAuc_list2=zeros(n,1);
AAup_list2=zeros(n,1);

PPrecision_list3=zeros(n,1);
RRecall_list3=zeros(n,1);
FFmeasure_list3=zeros(n,1);
AAuc_list3=zeros(n,1);
AAup_list3=zeros(n,1);

PPrecision_list4=zeros(n,1);
RRecall_list4=zeros(n,1);
FFmeasure_list4=zeros(n,1);
AAuc_list4=zeros(n,1);
AAup_list4=zeros(n,1);

PPrecision_list5=zeros(n,1);
RRecall_list5=zeros(n,1);
FFmeasure_list5=zeros(n,1);
AAuc_list5=zeros(n,1);
AAup_list5=zeros(n,1);
for time=1:n
    Precision_list1=zeros(per,1);
    Recall_list1=zeros(per,1);
    Fmeasure_list1=zeros(per,1);
    Auc_list1=zeros(per,1);
    Aup_list1=zeros(per,1);
    
    Precision_list2=zeros(per,1);
    Recall_list2=zeros(per,1);
    Fmeasure_list2=zeros(per,1);
    Auc_list2=zeros(per,1);
    Aup_list2=zeros(per,1);
        
    Precision_list3=zeros(per,1);
    Recall_list3=zeros(per,1);
    Fmeasure_list3=zeros(per,1);
    Auc_list3=zeros(per,1);
    Aup_list3=zeros(per,1);
    
    Precision_list4=zeros(per,1);
    Recall_list4=zeros(per,1);
    Fmeasure_list4=zeros(per,1);
    Auc_list4=zeros(per,1);
    Aup_list4=zeros(per,1);
    
    Precision_list5=zeros(per,1);
    Recall_list5=zeros(per,1);
    Fmeasure_list5=zeros(per,1);
    Auc_list5=zeros(per,1);
    Aup_list5=zeros(per,1);

    p=randperm(totalassociation); %ramdow total association
    for i=1:per
        Xn = X;
        W = X;
        W(W==0)=0.5;
        if i==per
            testset = p(((i-1)*totalper+1):totalassociation);
        else
            testset = p(((i-1)*totalper+1):i*totalper);
        end
        test_length = length(testset);
        true_list = zeros((test_length+zero_length),1);
        for c=1:test_length
            Xn(Index_PositiveRow(testset(c)),Index_PositiveCol(testset(c)))=0;
            W(Index_PositiveRow(testset(c)),Index_PositiveCol(testset(c)))=0;   
            true_list(c,1)=1;
        end
        [result1,result2,result3,result4,result5]=MY_PGdemo1(Xn,W,ln,Ls,n_parameter,s_parameter,u_parameter,v_parameter);
        result_list1 = create_resultlist( result1,testset,Index_PositiveRow,Index_PositiveCol,Index_zeroRow,Index_zeroCol,test_length,zero_length );
        result_list2 = create_resultlist( result2,testset,Index_PositiveRow,Index_PositiveCol,Index_zeroRow,Index_zeroCol,test_length,zero_length );
        result_list3 = create_resultlist( result3,testset,Index_PositiveRow,Index_PositiveCol,Index_zeroRow,Index_zeroCol,test_length,zero_length );
        result_list4 = create_resultlist( result4,testset,Index_PositiveRow,Index_PositiveCol,Index_zeroRow,Index_zeroCol,test_length,zero_length );
        result_list5 = create_resultlist( result5,testset,Index_PositiveRow,Index_PositiveCol,Index_zeroRow,Index_zeroCol,test_length,zero_length );
        [precision1,recall1,fmeasure1,auc1,aup1 ] = newmetric( true_list,result_list1 );
        [precision2,recall2,fmeasure2,auc2,aup2 ] = newmetric( true_list,result_list2 );
        [precision3,recall3,fmeasure3,auc3,aup3 ] = newmetric( true_list,result_list3 );
        [precision4,recall4,fmeasure4,auc4,aup4 ] = newmetric( true_list,result_list4 );
        [precision5,recall5,fmeasure5,auc5,aup5 ] = newmetric( true_list,result_list5 );
      
        Precision_list1(i,1)=precision1;
        Recall_list1(i,1)=recall1;
        Fmeasure_list1(i,1)=fmeasure1;
        Auc_list1(i,1)=auc1;
        Aup_list1(i,1)=aup1;
        
        Precision_list2(i,1)=precision2;
        Recall_list2(i,1)=recall2;
        Fmeasure_list2(i,1)=fmeasure2;
        Auc_list2(i,1)=auc2;
        Aup_list2(i,1)=aup2;
        
        Precision_list3(i,1)=precision3;
        Recall_list3(i,1)=recall3;
        Fmeasure_list3(i,1)=fmeasure3;
        Auc_list3(i,1)=auc3;
        Aup_list3(i,1)=aup3;
        
        Precision_list4(i,1)=precision4;
        Recall_list4(i,1)=recall4;
        Fmeasure_list4(i,1)=fmeasure4;
        Auc_list4(i,1)=auc4;
        Aup_list4(i,1)=aup4;
        
        Precision_list5(i,1)=precision5;
        Recall_list5(i,1)=recall5;
        Fmeasure_list5(i,1)=fmeasure5;
        Auc_list5(i,1)=auc5;
        Aup_list5(i,1)=aup5;     
    end
    PPrecision_list1(time,1)=mean(Precision_list1);
    RRecall_list1(time,1)=mean(Recall_list1);
    FFmeasure_list1(time,1)=mean(Fmeasure_list1);
    AAuc_list1(time,1)=mean(Auc_list1);
    AAup_list1(time,1)=mean(Aup_list1);
    
    PPrecision_list2(time,1)=mean(Precision_list2);
    RRecall_list2(time,1)=mean(Recall_list2);
    FFmeasure_list2(time,1)=mean(Fmeasure_list2);
    AAuc_list2(time,1)=mean(Auc_list2);
    AAup_list2(time,1)=mean(Aup_list2);
    
    PPrecision_list3(time,1)=mean(Precision_list3);
    RRecall_list3(time,1)=mean(Recall_list3);
    FFmeasure_list3(time,1)=mean(Fmeasure_list3);
    AAuc_list3(time,1)=mean(Auc_list3);
    AAup_list3(time,1)=mean(Aup_list3);
   
    PPrecision_list4(time,1)=mean(Precision_list4);
    RRecall_list4(time,1)=mean(Recall_list4);
    FFmeasure_list4(time,1)=mean(Fmeasure_list4);
    AAuc_list4(time,1)=mean(Auc_list4);
    AAup_list4(time,1)=mean(Aup_list4);
    
    PPrecision_list5(time,1)=mean(Precision_list5);
    RRecall_list5(time,1)=mean(Recall_list5);
    FFmeasure_list5(time,1)=mean(Fmeasure_list5);
    AAuc_list5(time,1)=mean(Auc_list5);
    AAup_list5(time,1)=mean(Aup_list5);

end

data1=[n_parameter,s_parameter, mean(PPrecision_list1), mean(RRecall_list1), mean(FFmeasure_list1),mean( AAuc_list1),mean(AAup_list1)]
data2=[n_parameter,s_parameter, mean(PPrecision_list2), mean(RRecall_list2), mean(FFmeasure_list2),mean( AAuc_list2),mean(AAup_list2)]
data3=[n_parameter,s_parameter, mean(PPrecision_list3), mean(RRecall_list3), mean(FFmeasure_list3),mean( AAuc_list3),mean(AAup_list3)]
data4=[n_parameter,s_parameter, mean(PPrecision_list4), mean(RRecall_list4), mean(FFmeasure_list4),mean( AAuc_list4),mean(AAup_list4)]
data5=[n_parameter,s_parameter, mean(PPrecision_list5), mean(RRecall_list5), mean(FFmeasure_list5),mean( AAuc_list5),mean(AAup_list5)]




