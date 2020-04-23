function [ result_list ] = create_resultlist( result,testset,Index_PositiveRow,Index_PositiveCol,Index_zeroRow,Index_zeroCol,test_length,zero_length )
    result_list = zeros((test_length+zero_length),1);      
    for i =1:test_length
        result_list(i,1)=result(Index_PositiveRow(testset(i)),Index_PositiveCol(testset(i)));
    end
    for i=1:zero_length
        result_list((i+test_length),1)=result(Index_zeroRow(i),Index_zeroCol(i));
    end
end

