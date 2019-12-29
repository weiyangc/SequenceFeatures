function baseCM_19features(inputpath, outputpath, OffsetWidth)
%   Read in DNA sequences from text file, compute the co-occurrence matrix and texture features.
%   Each DNA sequence must be stored as a text. The text cannot contain any extra characters, only A, C, G, T sequence data.
%   The used function GLCM_Features4 can be downloaded from: https://ww2.mathworks.cn/matlabcentral/fileexchange/22354-glcm_features4-m-vectorized-version-of-glcm_features1-m-with-code-changes


% Parameters %
% inputpath: the path of image, like '/path/';
% outputpath: results output path, like '/outputpath/';
% OffsetWidth: used in computing GLCM, default is 1;


if nargin == 2
    OffsetWidth = 1;
end

files_list3 = dir([inputpath '*.txt']); 

file_num = length( files_list3 ); 


% result file
%%%%%%%%  
result_txt = [outputpath 'baseCM_19features.txt'];
fid_txt = fopen(result_txt,'w');
fprintf(fid_txt, 'features\tAutocorrelation\tContrast\tCorrelation\tClusterP\tClusterS\tDissimilarity\tEnergy\tEntropy\tHomogeneity\tMaximumP\tSumSquares\tSumAverage\tSumVariance\tSumEntropy\tDifferenceVariance\tDifferenceEntropy\tInformation1\tInverseINN\tInverse3\n');


result_GLCM = [outputpath 'baseCM_GLCM.txt'];
fid_GLCM = fopen(result_GLCM,'w');
fprintf(fid_GLCM, 'GLCM\tG1\tG2\tG3\tG4\tG5\tG6\tG7\tG8\tG9\tG10\tG11\tG12\tG13\tG14\tG15\tG16\n');
%%%%%%%%%



for i_file3 = 1:length( files_list3 )
	sequnce1 = textread([inputpath files_list3(i_file3).name], '%s');
	s1_len = length( sequnce1{1} );
	
	% translate A,C,G,T to 1,2,3,4
	s1_numeric = [];   % store one sequence's numeric result
	for i_base = 1:s1_len
		if sequnce1{1}(i_base) == 'A'
			s1_numeric = [s1_numeric 1];
		elseif sequnce1{1}(i_base) == 'C'
			s1_numeric = [s1_numeric 2];
		elseif sequnce1{1}(i_base) == 'G'
			s1_numeric = [s1_numeric 3];
		elseif sequnce1{1}(i_base) == 'T'
			s1_numeric = [s1_numeric 4];
		else
			disp('1 base is not A,C,G,T;');
		end
	end
	
	
	
	glcm = graycomatrix(s1_numeric, 'GrayLimits',[1 4], 'NumLevels', 4, 'Offset', [0 1], 'Symmetric', true);
	
	% 
	% the function GLCM_Features4 can be downloaded from: https://ww2.mathworks.cn/matlabcentral/fileexchange/22354-glcm_features4-m-vectorized-version-of-glcm_features1-m-with-code-changes
	stats = GLCM_Features4(glcm, 0);   
	
	
	
		% Autocorrelation: stats.autoc;
		% Contrast: stats.contr;
		% Correlation: stats.corrm;
		% ClusterP: stats.cprom;    % 
		% ClusterS: stats.cshad;   % 
		% Dissimilarity: stats.dissi;   % 
		% Energy: stats.energ;
		% Entropy: stats.entro;
		% Homogeneity: stats.homom;
		% MaximumP: stats.maxpr;
		% SumSquares: stats.sosvh;
		% SumAverage: stats.savgh;
		% SumVariance: stats.svarh;
		% SumEntropy: stats.senth;
		% DifferenceVariance: stats.dvarh;
		% DifferenceEntropy: stats.denth;
		% Information1: stats.inf1h;
		% InverseINN: stats.indnc;
		% Inverse3: stats.idmnc;
	
	
	% write one sequence's results to the result file
	onlyname = files_list3(i_file3).name(1 : length(files_list3(i_file3).name)-4 );	
	fprintf(fid_txt, '%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n' ,  ...
	                  onlyname, stats.autoc, stats.contr, stats.corrm,  stats.cprom, stats.cshad, stats.dissi, stats.energ, stats.entro, stats.homom, stats.maxpr,  ...
					  stats.sosvh, stats.savgh, stats.svarh, stats.senth, stats.dvarh, stats.denth, stats.inf1h, stats.indnc, stats.idmnc );
	
	% write GLCM to file
	glcm2 = glcm'; 
	glcmVector = glcm2(:)';  % one row
	fprintf(fid_GLCM, '%s', onlyname );
	for i_glcm = 1:length(glcmVector)
		fprintf(fid_GLCM, '\t%f', glcmVector(i_glcm) );
	end
	fprintf(fid_GLCM, '\n');
	
end

% %

fclose(fid_txt);
fclose(fid_GLCM);

end








