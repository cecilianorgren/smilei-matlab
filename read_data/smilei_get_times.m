function out = smilei_get_times(filePath)
fileInfo = h5info(filePath);
nOutput = numel(fileInfo.Groups.Groups);

for iOutput = 1:nOutput
  time(iOutput) = fileInfo.Groups.Groups(iOutput).Attributes(1).Value;
end

out = time;