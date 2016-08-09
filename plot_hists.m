myFolder = 'C:\Users\Cris\Desktop\Fermion-Transport\Measurements\conditional mutual info disorder 7-18-2016';

filePattern = fullfile(myFolder, '*.mat');
theFiles = dir(filePattern);
for k = 1 : length(theFiles)
  baseFileName = theFiles(k).name;
  fullFileName = fullfile(myFolder, baseFileName);
  [~, filename, ~] = fileparts(fullFileName);
  data = load(fullFileName);
  data = cell2mat(struct2cell(data));
  
  for i=1:size(data(:), 1)
    if data(i) == 0
        data(i) = 1.4211E-14;
    end
  end
  
  site = 51;
  nbins = 20;
  
  h = figure;
  hist(data(:, site));
  saveas(h, ['Hist_', filename, '.jpeg']);
end