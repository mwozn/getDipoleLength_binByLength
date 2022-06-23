% Classify subTOM MOTL according to dipole_data_lengthBin_#.txt
%
% Michael Wozny 20220610
%
% According to subTOM convention, column 20 is reserved for classes:
% 1 = ALWAYS
% 2 = NEVER (also negative numbers)
% 3 = first distinct class
% 4 = second distinct class
% ... = ... distinct class

allmotl = 'allmotl_1.em';
allmotl = dread(allmotl);

% Reset class to 2 (only use subtomos for dipoles in input files)

allmotl(20,:) = 2;

% Get the dipole_data_lengthBin files from dipole_length_sortByLength.m
input_lengthBins = dir('*dipole_data_lengthBin_*');

% Match and cassify
for k=1:length(input_lengthBins)
    lengthBin = csvread(strcat('dipole_data_lengthBin_',int2str(k),'.txt'));
    for m = 1:length(allmotl)
        for q=1:length(lengthBin)
            % if MOTL subtomo ID is same as input subtomo ID
            if allmotl(4,m) == lengthBin(q,9)
               allmotl(20,m) = 2 + k;
            end
        end
    end
end

% OUTPUT
dwrite(allmotl, 'allmotl_classByLength_1.em');