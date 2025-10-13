% Original file path
fname = 'C:\Users\oranm\OneDrive - Imperial College London\PHD\HSSCC\QUANTA\HOLE12_EBSD_100625\oran_hole12_hdf5.h5';

% Output CTF file path
ctf_fname = 'C:\Users\oranm\OneDrive - Imperial College London\PHD\HSSCC\QUANTA\HOLE12_EBSD_100625\locA_73.36.ctf';

% Read EBSD data
phi1  = h5read(fname, '/OVERNIGHT_EBSD_hdf5/EBSD/Data/phi1');
PHI   = h5read(fname, '/OVERNIGHT_EBSD_hdf5/EBSD/Data/PHI');
phi2  = h5read(fname, '/OVERNIGHT_EBSD_hdf5/EBSD/Data/phi2');
x     = h5read(fname, '/OVERNIGHT_EBSD_hdf5/EBSD/Data/X BEAM');
y     = h5read(fname, '/OVERNIGHT_EBSD_hdf5/EBSD/Data/Y BEAM');
phase = h5read(fname, '/OVERNIGHT_EBSD_hdf5/EBSD/Data/Phase');

% Read step sizes
xStep = h5read(fname, '/OVERNIGHT_EBSD_hdf5/EBSD/Header/XSTEP');
yStep = h5read(fname, '/OVERNIGHT_EBSD_hdf5/EBSD/Header/YSTEP');

% Flatten arrays (column vectors)
phi1 = phi1(:); PHI = PHI(:); phi2 = phi2(:);

x = x(:); y = y(:); phase = phase(:);

% Prepare dummy columns
num_points = length(phase);
Bands = zeros(num_points, 1);    % dummy 0
Error = zeros(num_points, 1);    % dummy 0
MAD = zeros(num_points, 1);      % dummy 0
BC = 255 * ones(num_points, 1);  % max contrast 255
BS = 255 * ones(num_points, 1);  % max slope 255

% Combine all columns (order must match header)
data_matrix = [ ...
    phase, ...
    x, ...
    y, ...
    Bands, ...
    Error, ...
    phi1, ...
    PHI, ...
    phi2, ...
    MAD, ...
    BC, ...
    BS ...
];

% Updated header including dummy columns
header_lines = {
    'Channel Text File'
    'Prj unnamed'
    'Author	[Unknown]'
    'JobMode	Grid'
    ['XCells	', num2str(numel(unique(x)))]
    ['YCells	', num2str(numel(unique(y)))]
    ['XStep	', num2str(xStep, '%.6f')]
    ['YStep	', num2str(yStep, '%.6f')]
    'AcqE1	0'
    'AcqE2	0'
    'AcqE3	0'
    'Euler angles refer to Sample Coordinate system (CS0)!	Mag	3551.000000	Coverage	100	Device	0	KV	10.000000	TiltAngle	70	TiltAxis	0'
    'Phases	2'
    '3.305000;3.305000;3.305000	90.000000;90.000000;90.000000	titanium-beta	11	229'
    '2.950000;2.950000;4.683000	90.000000;90.000000;120.000000	titanium-alpha	9	194'
    'Phase	X	Y	Bands	Error	Euler1	Euler2	Euler3	MAD	BC	BS'
};

% Write to CTF file
fid = fopen(ctf_fname, 'w');
for i = 1:numel(header_lines)
    fprintf(fid, '%s\n', header_lines{i});
end
for i = 1:size(data_matrix, 1)
    fprintf(fid, '%d\t%.4f\t%.4f\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%d\t%d\n', data_matrix(i, :));
end
fclose(fid);

disp('CTF file written with dummy columns for Bands, Error, MAD, BC, BS, and real EBSD data.');
