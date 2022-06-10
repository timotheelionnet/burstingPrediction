% inputs are a folder full of z-stacks, each with a collection of objects,
% for instance nulei, each nucleus with its own ID (i.e. intensity value).
% For each input image, the script loops through objects, crops them, 
% 2D projects them and saves in own image per each.
% The script removes any object other than the chosen one in the cropped stack

%% change folders locations!!!
% where the input z-stacks are located - should me masks, each object (e.g. nuclei)
% with a unique ID.
dirName = '\\research-cifs.nyumc.org\research\lionnt01lab\lionnt01labspace\Nich_M\Imaging\Nikon\20211019\masks\312b';

% where the output of the nuclei will be saved 
outDir = '\\research-cifs.nyumc.org\research\lionnt01lab\lionnt01labspace\Nich_M\Imaging\Nikon\20211019\masks\crop\312b\C2-HoxA5';

%%
% how many voxels padding around the object you want
% (will not add padding if the object touches the image border)
paddingSize = 15;

% how you recoginze the files to load - They need to have a specific extension
imgExtension = '.tiff';
%%
if ~exist(outDir,'dir')
    mkdir(outDir);

end
%mask here
ImgMetadata = imfinfo('\\research-cifs.nyumc.org\research\lionnt01lab\lionnt01labspace\Nich_M\Imaging\Nikon\20211019\masks\312b\C1-MAX_HCR_test-HOX-312-0505_1_MMStack_Pos0.ome_mask.tiff');

flist = dir(dirName); % get list of files in input dir
for i=1:numel(flist)
    if endsWith(flist(i).name,imgExtension)
        
        % load stack
        fname = fullfile(dirName,flist(i).name);
        %disp(['Opening ',fname,' ...']);
        img = timtiffread(fname);
        
        % get list of objects IDs
        nucleiList = unique(img(:));
        
        % loop through objects, crop, 2D project each and save as its own image
        for j=1:numel(nucleiList)
            
            % collect coordinates of current object
            idxNucleus = find(img == nucleiList(j));
            [x,y] = ind2sub(size(img),idxNucleus);
            
            % crop stack around current object, adding desired patdding
            minx = max(1,min(x)-paddingSize);
            miny = max(1,min(y)-paddingSize);
            
            maxx = min(size(img,1),max(x)+paddingSize);
            maxy = min(size(img,2),max(y)+paddingSize);
            
            arr = zeros(numel(minx:maxx),numel(miny:maxy),numel(ImgMetadata));
            disp(['Starting cell ID ',num2str(j),' ...']);
            for v=1:numel(ImgMetadata)
                %image to be cropped here
                fullimg = imread('\\research-cifs.nyumc.org\research\lionnt01lab\lionnt01labspace\Nich_M\Imaging\Nikon\20211019\splits\312b\c2\C2-MAX_HCR_test-HOX-312-0505_1_MMStack_Pos0.ome.tif',v);
                %disp(['Opening time point ', num2str(v),' ...']);
                
                crop = fullimg(minx:maxx,miny:maxy);
                arr(:,:,v) = crop(:,:);v;
            end
                fname = [strrep(flist(i).name,imgExtension,''), ...
                        ['CellID',num2str(nucleiList(j)),'.tif']];   
                %disp(['Saving cell ID ', num2str(j),' ...']);    
                save_as_tiff(arr,fullfile(outDir,fname),'.tif');
        end
    end
end