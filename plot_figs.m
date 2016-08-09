myFolder = 'C:\Users\Cris\Desktop\Fermion-Transport\Measurements\New folder';

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
  data = mean(abs(data), 1);
  sizeC_values = 2*[1:99];
  %figure;
  disp(filename);
  scatter(sizeC_values, log(abs(data)), 25, 'filled')
  xlabel('Size of region B');
  ylabel('Log of CMI(A:C|B)');
  legend('beta = 0.01', 'beta = 0.05', 'beta = 0.1', 'beta = 0.5', 'beta = 1', 'beta = 2', 'beta = 5');
  title('CMI for Mu = 0, Disorder = 0.1');
  hold on;
  %savefig([filename, '.fig']);
end