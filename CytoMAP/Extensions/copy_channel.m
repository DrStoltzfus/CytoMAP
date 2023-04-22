filename = '13-1 auLNslide Settings2 LN1.ims';
h5disp(filename);

% Try to copy Channel 11 -> Channel 12
% Which would apply to all the simple channel computation like summary of
% severial channels.
pre = '/DataSet/ResolutionLevel';
ending = '/TimePoint 0/Channel 11';
% Give the name to the new channell - Channel 12
new_ending = '/TimePoint 0/Channel 12';
Datasets = ['/Data', '/Histogram'];
for i = 0:5
    for Dataset = ["/Data", "/Histogram"]
        % Write Data
        location = strcat(pre, " ", num2str(i), ending, Dataset)
        data = hdf5read(filename, location);
        % Specific the location for new data
        new_location = strcat(pre, " ", num2str(i), new_ending);
        dset_tests.Location = new_location;
        name = Dataset.split('/');
        dset_tests.Name = name(2);
        dset_tests
        hdf5write(filename, dset_tests, data, 'WriteMode','append');
%         
%         % Write attrs
%         for attr = ['ImageSizeX', 'ImageSizeY', 'ImageSizeZ', 'HistogramMin', 'HistogramMax']
%             location = strcat(pre, " ", num2str(i), ending);
%             new_location = strcat(pre, " ", num2str(i), new_ending);
%             val = h5readatt(filename, location, attr);
%             hdf5writeatt(filename, new_location, attr, val)
%         end
    end
end

h5disp(filename);