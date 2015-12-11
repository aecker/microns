% Receptive field analysis on van Gogh noise
% AE 2015-12-11

maps = fetchn(vangogh.RF, 'map');


%%
m = cellfun(@(m) abs(m(:, :, 2)), maps, 'uni', false);
m = cat(3, m{:});

% compute significance score
%   we take the spatial variance of the second bin, which usually has the
%   best signal divided by that of the last bin, which is typically noise.
ii = 10 : 40;
jj = 60 : 110;
v = cellfun(@(m) var(reshape(m(ii, jj, :), [], size(m, 3))), maps, 'uni', false);
v = cat(1, v{:});
r = v(:, 2) ./ v(:, end);
% hist(r, 0.5 : 0.1 : 5)

[~, order] = sort(r, 'descend');

%%
total = numel(order);
M = 5;
N = 4;
i = 0;
while i < total
    k = mod(i, M * N);
    if ~k
        figure
    end
    subplot(M, N, k + 1)
    map = maps{order(i + 1)}(:, :, 2);
    imagesc(map - median(map(:)))
    caxis([-1 1] * max(abs(caxis)))
    colormap(ne7.vis.doppler)
    axis equal off
    i = i + 1;
end

