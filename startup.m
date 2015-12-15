function startup

p = fileparts(mfilename('fullpath'));
addpath(p)
addpath(fullfile(p, 'util'))
addpath(fullfile(p, 'analysis'))
p = fileparts(p);
addpath(fullfile(p, 'ald'))
addpath(fullfile(p, 'ald/tools'))
addpath(fullfile(p, 'commons/schemas'))
addpath(fullfile(p, 'commons/lib'))
addpath(fullfile(p, 'commons/visual-stimuli'))
addpath(fullfile(p, 'dimitri'))
addpath(fullfile(p, 'jake'))
addpath(fullfile(p, 'oopsi'))
addpath(fullfile(p, 'vangoghstim'))
