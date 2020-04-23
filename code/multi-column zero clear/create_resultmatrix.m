function [ result ] = create_resultmatrix( result,testset,prolist,nogo )
leave_col=prolist(testset);
result = result(:,leave_col);
result(nogo,:)=[];
end

