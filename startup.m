function startup

p = fileparts(mfilename('fullpath'));
addpath(p)
p = fileparts(p);
addpath(fullfile(p, 'jake'))
addpath(fullfile(p, 'dimitri'))
addpath(fullfile(p, 'commons/schemas'))
addpath(fullfile(p, 'commons/lib'))
addpath(fullfile(p, 'commons/visual-stimuli'))
addpath(fullfile(p, 'vangoghstim'))
addpath(fullfile(p, 'oopsi'))
