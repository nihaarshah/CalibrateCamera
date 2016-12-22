% function [mat] = test()
% for row = 1:3
% for col = 1:3
% mat(row,col) = rand;
% end
% end
% mat
function [cosindex] = test()
nBestConsensus = 10;
for row = 1:2:((2*nBestConsensus)-1)
    cosindex = (row+1)/2
    
end

end