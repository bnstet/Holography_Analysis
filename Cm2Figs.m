function [mov,nFrames] =Cm2Figs(A, C, b,f, startIdx,endIdx)

Ashape = size(A); Cshape = size(C); bshape = size(b); fshape = size(f);
nCells = Ashape(1);
nFrames = min( endIdx - startIdx + 1, Cshape(1) - startIdx + 1) ;
endIdx = startIdx + nFrames - 1;
width = Ashape(3);
height = Ashape(2);
nbg = bshape(1);

Atemp = permute(A, [3 2 1]);
Atemp = reshape(Atemp,[], nCells);
Ctemp = transpose(C(startIdx:endIdx,:));


btemp = permute(b, [3 2 1]);
btemp = reshape(btemp,[], nbg);
ftemp = transpose(f(startIdx:endIdx,:));

mov = Atemp*Ctemp + btemp*ftemp;
mov = reshape(mov,width, height, nFrames);

end

