% Euclidean distance between dipole north & south points from Dynamo models.
%
% Michael Wozny 20220610
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USAGE: 1. index_number_all = model numbers with the dipoles to be measured
%        2. cat_tomo_vol_dir = path to dir with /volume_#/models/*.omd
%        3. Enter pixel_size and bin to report in nm/Angstroms, set to 1 for px
%        4. Plot dipole_data.txt with dipole_length.r
%
% For lists binned by dipole length:
%        1. set binDataByLengths = 1
%        2. set lengthBins = to the number of bins

% User input for general usage
pathToOutput = 'sortByLength/'
index_number_all = [6,8,10,12,16,18,24,26,28,30,34,40,42,44,48,52,54,56,58,60,62,64,66,70,72,76,78,80,82,84,88,92,94,96,98,100,102,104,106,110,112,114,118,120,122,124,126,128,130,132,134];
cat_tomo_vol_dir = 'Cat/tomograms';
pixel_size = 0.2684;
bin = 2;

% User input for sorted .csv output
binDataByLengths = 1;
lengthBins = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT: dipole_data.txt
% 1 - dipole length
% 2 - dynamo tomo idx
% 3 - tomo running number (ordered by dynamo tomo idx from lowest to highest)
% 4 - dipole model center value empty? (0 = crop point, 1 = not a crop point)
% 5 - dipole id (running number of dipoles, includes those without crop)
% 6 - crop point x coordinate
% 7 - crop point y coordinate
% 8 - crop point z coordinate
% 9 - particle id (running number of dipoles with crop, 0 = not a crop point)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dipole_length = cell(1,length(index_number_all));

for q = 1:length(index_number_all)
    % input model index number corresponding to the tomograms in Catalogue .vll
    index_number = int2str(index_number_all(q));
    model_path = strcat(cat_tomo_vol_dir,'/volume_',index_number,'/models/');
    model_path_files = dir(model_path);
    model_file = model_path_files(3).name;
    model = dread(strcat(model_path,model_file));
    
    % loop over model.dipoles, store center value of each dipole in
    % consecutive order as model_table
    norm_dist = [];
    
    % keep track of tomo id and crop points to find particle id later
    particle_list = [];
    empty_dipole_log = [];
    emptyDipole = 0;
    for k = 1:length(model.dipoles)
        if isempty(model.dipoles{1,k}.center) || isempty(model.dipoles{1,k}.north) || isempty(model.dipoles{1,k}.south)
            disp('q is: ')
            disp(q)
            disp('k is: ')
            disp(k)
            %norm(model.dipoles{1,k}.north - model.dipoles{1,k}.south)
            disp('tomgram index number is:')
            disp(index_number_all(q))
            disp('isempty return:')
            disp(isempty(model.dipoles{1,k}.center))
            disp('center value:')
            disp(model.dipoles{1,k}.center)
            emptyDipole = emptyDipole + 1;
            disp('empty dipole:')
            disp(emptyDipole);
            empty_dipole_holder(1,1) = k;
            empty_dipole_holder(1,2) = index_number_all(q);
            continue
        end
        kk = k - emptyDipole;
        norm_dist(kk,1) = norm(model.dipoles{1,k}.north - model.dipoles{1,k}.south);
        
        % build for all_particle_list which is used for crop points per
        % tomogram (plot of number of bridges per tomo)
        particle_list(k,1) = norm(model.dipoles{1,k}.north - model.dipoles{1,k}.south);
        particle_list(k,2) = index_number_all(q);
        particle_list(k,3) = q;
        particle_list(k,4) = isempty(model.dipoles{1,k}.center);
        particle_list(k,5) = k;
        
        % crop points xyz
        particle_list(k,6) = model.crop_points(k,1);
        particle_list(k,7) = model.crop_points(k,2);
        particle_list(k,8) = model.crop_points(k,3);
        
        % particle id
        if q == 1
            particle_list(k,9) = kk;
        elseif q > 1
            particle_list(k,9) = kk + (particle_all_list{q-1}(size(particle_all_list{q-1},1),9));
        else
            disp('something funny happened...')
        end
    end
    
    % round to 4 digits and transpose arrary
    norm_dist = round(norm_dist,4);
    
    particle_list(:,1) = round (particle_list(:,1),4);

    dipole_length{q} = (norm_dist * (pixel_size * bin));
    particle_list(:,1) = (particle_list(:,1) * (pixel_size * bin));
    particle_all_list{q} = particle_list;
    empty_dipole_log_all{q} = empty_dipole_log;
end

% concatenate NEAREST_NEIGHBOUR_ALL to have all nearest neighbour
% measurements together
all = cat(1, dipole_length{:});

% build for all_particle_list which is used for crop points per
% tomogram (plot of number of bridges per tomo)
all_particle_list = cat(1, particle_all_list{:});

% for display
disp(strcat(['Median = ',num2str(median(all))]));

% report bin with max counts
binWidth = 1;
counts = histcounts(all,'BinWidth',binWidth,'BinLimits',[0,100]);
[row,col]=find(ismember(counts,max(counts)));
disp(strcat(['Hist peak = ',int2str((col * binWidth)-binWidth),'-',int2str(col * binWidth)]));

% plot nearest neighbour measurements as a histogram
figure('Renderer', 'painters', 'Position', [10 10 1200 600]);
histogram(all,'BinWidth',binWidth,'BinLimits',[0,100]);
xlabel('Dipole Length (nm)');
set(gca,'XTick',0:5:100);
ylim([0 45]);

% OUTPUT, general

csvwrite(strcat(pathToOutput,'dipole_data.txt'),all_particle_list);

% OUTPUT, sorted by length
if binDataByLengths == 1
    
  all_particle_list_sortByLength = sortrows(all_particle_list,1);
  
  for l=1:lengthBins
      disp('lower, upper bin edges');
      if l == 1
          lowerBinEdge = 1
          upperBinEdge = floor(((size(all_particle_list_sortByLength,1)/lengthBins)))-1
      elseif l == lengthBins
          lowerBinEdge = floor(((size(all_particle_list_sortByLength,1)/lengthBins))) * (l-1)
          upperBinEdge = size(all_particle_list_sortByLength,1)
      elseif l > 1 && l < lengthBins
          lowerBinEdge = floor(((size(all_particle_list_sortByLength,1)/lengthBins))) * (l-1)
          upperBinEdge = (floor(((size(all_particle_list_sortByLength,1)/lengthBins))) * l)-1
      else
          disp('something funny with the length bins');
      end
      all_particle_list_sortByLength_binned = all_particle_list_sortByLength(lowerBinEdge:upperBinEdge,:);
      all_particle_list_sortByLength_binned = sortrows(all_particle_list_sortByLength_binned,9);
      disp('min length');
      disp(min(all_particle_list_sortByLength_binned(:,1)));
      disp('max length');
      disp(max(all_particle_list_sortByLength_binned(:,1)));
      disp(strcat(['Median = ',num2str(median(all_particle_list_sortByLength_binned(:,1))')]));
      sortByLengthBinFileName = strcat(pathToOutput,'dipole_data_lengthBin_',int2str(l),'.txt')
      csvwrite(sortByLengthBinFileName,all_particle_list_sortByLength_binned);
  end
end