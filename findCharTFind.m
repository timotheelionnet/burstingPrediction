%%Find indices of characterized TFs in "all TF" list
function [charTFind, burstingData2] = findCharTFind(charTFs, tfList, burstingData)
charTFind = [];
for v = 1:length(charTFs)
charTFind = [charTFind, strmatch(charTFs(v), tfList, 'exact')];
end

retained = [];
for f = 1:numel(charTFind)
retained = [retained, strmatch(tfList(charTFind(f)), burstingData.TFnames,'exact')];
end
burstingData2 = burstingData(retained,:);

end