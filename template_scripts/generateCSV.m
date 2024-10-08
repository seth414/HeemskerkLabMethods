% read in position data
path2position = 'positions.mat';
load(path2position);

% file name for saving
expName = 'BMP_2024';

% channel used for tracking
ch_i = 1;

% generate csv files for trackmate
for i = 1:length(positions)
    colony = positions(i);
    ID = (1:sum(colony.ncells))';
    pos = vertcat(colony.cellData.XY); % make sure position's unit is correct
    chs = vertcat(colony.cellData.nucLevel);
    POSITION_X = pos(:, 1);
    POSITION_Y = pos(:, 2);
    FRAME = repelem((1:colony.nTime) - 1, colony.ncells)'; % trackmate is 0-based
    MEAN_INTENSITY_CH = chs(:, ch_i);
    dataTable = table(ID, POSITION_X, POSITION_Y, FRAME, MEAN_INTENSITY_CH);
    writetable(dataTable, [expName '_colony_' num2str(i) '_tracking.csv'], 'WriteRowNames', false);
end