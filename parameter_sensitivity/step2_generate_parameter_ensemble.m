function step2_generate_parameter_ensemble(filename,nensemble)

% input, 1. params_range.txt: file contain the parameter range
% input, 2. 7: how many ensemble parameter (2^7) want to generate
% output, parameters are saved in ALM_param.txt

Nrows = numel(textread(filename,'%1c%*[^\n]')); % determine how many parameter to tune

fileID = fopen(filename);
C = textscan(fileID,'%s %f %f',Nrows);
fclose(fileID);

sample_normalized = HOSobol(nensemble,Nrows,1); % sobol samples [0 1]

Y = [];
for i = 1 : Nrows
    Y = [Y repmat(C{1,3}(i,1) - C{1,2}(i,1),length(sample_normalized),1).* sample_normalized(:,i) + repmat(C{1,2}(i,1),length(sample_normalized),1)]; % actual sobol samples [lower upper] 
end

fid = fopen('ALM_param.txt','w');
for j = 1 : Nrows
    fprintf(fid,'%s',cell2mat(C{1,1}(j,1)));
    if j < Nrows
        fprintf(fid,'\t');
    else
        fprintf(fid,'\n');
    end
end
for i =1 : length(Y)
    for j = 1 : Nrows
        fprintf(fid,'%e',Y(i,j));
        if j < Nrows
            fprintf(fid,'\t');
        else
            fprintf(fid,'\n');
        end
    end
end

end