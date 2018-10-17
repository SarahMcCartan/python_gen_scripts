
%% Function to catch inf/nan error in a transition matric / pm /  km  & replace with 0

function f = infNan(outF,ns,mat,k,temp,pep)
for aa = 1: ns
    for bb = 1:ns
        if isnan(mat(aa,bb)) == 1 || isinf(mat(aa,bb)) ==1
            mat(aa,bb) = 0;
            fileID = fopen(sprintf('%s\\%d%serrorLog%d.txt',outF, temp,pep,k), 'w');
            formatSpec = '%d_TM(%d,%d) = NaN/Inf\n';
            fprintf(fileID,formatSpec,k,aa,bb)
        end
    end
end

f = mat;
end
