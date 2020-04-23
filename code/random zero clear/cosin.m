function [ B ] = cosin( A )
%只是用于特定数据，就是0，1矩阵的余弦相似性
[m,n] = size(A);
B = zeros(m,m);
for i=1:m
    for j=i:m
        up=dot(A(i,:),A(j,:));
        down=sqrt(sum(A(i,:)))*sqrt(sum(A(j,:)));
        B(i,j)=up/down;
    end
    B(i,i)=0;
end
B = B+B';
end

