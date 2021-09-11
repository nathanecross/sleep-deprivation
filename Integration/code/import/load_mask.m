
function mask = load_mask(masks_dir, regionOfInterest)

 files = struct2table(dir(masks_dir)); files = files.name;
 mask_filename = files(contains(files,char(regionOfInterest)));
 mask_path = [masks_dir char(mask_filename)];
 mask = load_nii(mask_path);
 mask = logical(mask.img);
 
 end
 

