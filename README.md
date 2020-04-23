#PONMF-S

Predicting protein functions by using non-negative matrix factorization with multi-networks co-regularization

Developer:Lun Li from Kunming University of Science and Technology.

## Instructions to PONMF-S(version 2.0.0)


Requirement
-----------------------------------------------------------------------------------------------------------------------
* 8GB memory
* MATLAB R2016a or later

------------------------------------------------------------------------------------------------------------------------




Input data descriptions
------------------------------------------------------------------------------------------------------------------------


The folder './data/newyeastdata/P' store the yeast data for BP process.
The folder './data/newyeastdata/C' store the yeast data for CCprocess.
The folder './data/newyeastdata/F' store the yeast data for MFprocess.
The folder './data/humandata/P' store the human data for BP process.
The folder './data/humandata/C' store the human data for CC process.
The folder './data/humandata/F' store the human data for MFprocess.

The file './data/newyeastdata/P/newPgp.txt' is the matrix of protein and Goterm association.The entry Xij equaling 1 represents connection between the i-th go and j-th protein related to the list above, and 0 otherwise.
The file './data/newyeastdata/P/newPgogo.txt' is the adjanceny matrix of the interaction newtwork of go and go.
The file './data/newyeastdata/P/newPPI.txt' is the adjanceny matrix of the interaction newtwork of protein and protein.
------------------------------------------------------------------------------------------------------------------------





parameter interpretation
------------------------------------------------------------------------------------------------------------------------

n_parameter:The tuning parameter of interactome information for integrating ppi network information into the model.

s_parameter:The tuning parameter of interactome information for integrating gogo network information into the model.

u_parameter:The tuning parameter of Frobenius norm based regularization to prevent over-fitting problem.

v_parameter:The tuning parameter of Frobenius norm based regularization to prevent over-fitting problem.

k:protein or goterm representation in the K dimension.



Output data descriptions
------------------------------------------------------------------------------------------------------------------------
result1:PONMF-S
result2:PONMF-S1
result3:PONMF-S2
result4:PONMF-S3
result5:PONMF-S4
------------------------------------------------------------------------------------------------------------------------
