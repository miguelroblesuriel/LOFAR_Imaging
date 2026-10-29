% For every saved map figure (*.fig in the current folder):
%   - read the image data back (no need to rerun the imaging),
%   - replot it with the sky convention (x = l, East on the left, y = m),
%   - mark the catalog positions of Cygnus A and Cassiopeia A,
%   - detect the two strongest peaks and compute the position errors,
%   - save the new figure (.fig and .png) and a CSV table with all errors.
clear all
close all
clc

%% Settings
outDir  = 'with_refs';
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% Data and grid (same FracAngle as the imaging script)
load lofar_DSP_data_1.mat
c      = 299792458;
lambda = c / freq;
dist = sqrt(  (meshgrid(poslocal(:,1)) - meshgrid(poslocal(:,1)).').^2 ...
            + (meshgrid(poslocal(:,2)) - meshgrid(poslocal(:,2)).').^2);
D   = max(dist(:));
res = lambda / D;                       % resolution limit in (l,m) units
dl  = 0.25 * res;
l   = 0:dl:1;
l   = [-fliplr(l(2:end)), l]';
m   = l;
len = length(l);

%% Reference positions of the sources
names = {'Cyg A', 'Cas A'};
% Stellarium gives (azimuth from North via East, elevation) in degrees.
dms = @(d, mi, s) d + mi/60 + s/3600;
az = [dms(108,11,51.8), dms(48,9,52.7)];     % Cyg A, Cas A
el = [dms( 67,10,20.3), dms(46,53, 3.7)];

% Direction cosines for this data set: l = East, m = North
ref.l = cosd(el) .* sind(az);
ref.m = cosd(el) .* cosd(az);
fprintf('Reference positions (Stellarium):\n');
for j = 1:2
    fprintf('  %s: az = %.2f, el = %.2f deg -> l = %.4f, m = %.4f\n', ...
            names{j}, az(j), el(j), ref.l(j), ref.m(j));
end

%% Loop over all saved figures
files = dir('*.fig');
rows  = {};                              % results table
for q = 1:numel(files)
    [~, name] = fileparts(files(q).name);

    % --- read data from the figure ------------------------------------------
    fig  = openfig(files(q).name, 'invisible');
    hImg = findobj(fig, 'Type', 'image');
    if isempty(hImg)
        warning('%s: no image found, skipped.', files(q).name);
        close(fig); continue
    end
    hImg  = hImg(1);
    img   = hImg.CData;
    xdata = hImg.XData;
    ax    = hImg.Parent;
    cmap  = colormap(ax);
    clims = ax.CLim;
    close(fig);

    if ~isequal(size(img), [len len])
        warning('%s: image size does not match the grid, skipped.', files(q).name);
        continue
    end
    % Old figures: imagesc(img) (XData = [1 len]), rows = l, cols = m.
    % Figures already saved with l/m axes: XData starts at l(1), rows = m.
    if abs(xdata(1) - l(1)) < 1e-6
        imgLM = img.';                   % convert to rows = l, cols = m
    else
        imgLM = img;
    end

    % --- find the two strongest, separate peaks ---------------------------------
    [L, M] = ndgrid(l, m);               % L(a,b) = l(a), M(a,b) = m(b)
    work = imgLM;
    pk = zeros(2, 3);                    % [l  m  value]
    for t = 1:2
        [val, ix] = max(work(:), [], 'omitnan');
        pk(t,:) = [L(ix), M(ix), val];
        work(hypot(L - L(ix), M - M(ix)) < 5*res) = NaN;   % mask this peak
    end

    % --- match peaks to sources (assignment with the smallest total error) -------
    e = zeros(2,2);                      % e(i,j): peak i to source j
    for i = 1:2
        for j = 1:2
            e(i,j) = hypot(pk(i,1) - ref.l(j), pk(i,2) - ref.m(j));
        end
    end
    if e(1,1) + e(2,2) <= e(1,2) + e(2,1), pm = [1 2]; else, pm = [2 1]; end
    % pm(j) = index of the peak assigned to source j

    err = zeros(1,2); angErr = zeros(1,2);
    for j = 1:2
        pl = pk(pm(j),1); pmm = pk(pm(j),2);
        err(j) = hypot(pl - ref.l(j), pmm - ref.m(j));
        sp = [pl;       pmm;       sqrt(max(0, 1 - pl^2 - pmm^2))];
        sr = [ref.l(j); ref.m(j);  sqrt(max(0, 1 - ref.l(j)^2 - ref.m(j)^2))];
        angErr(j) = acosd(min(1, sp.' * sr));          % angle between directions
    end

    % --- new figure ---------------------------------------------------------------
    nf = figure('Visible', 'off');
    h  = imagesc(l, m, imgLM.');
    set(h, 'AlphaData', ~isnan(imgLM.'));
    axis xy equal tight;
    set(gca, 'XDir', 'reverse');         % sky view: East on the left
    colormap(cmap); caxis(clims); colorbar;
    xlabel('East \leftarrow  l  \rightarrow West');
    ylabel('South \leftarrow  m  \rightarrow North');
    hold on;
    hr = plot(ref.l, ref.m, 'wo', 'MarkerSize', 14, 'LineWidth', 1.5);
    hp = plot(pk(:,1), pk(:,2), 'rx', 'MarkerSize', 10, 'LineWidth', 1.5);
    for j = 1:2
        plot([ref.l(j) pk(pm(j),1)], [ref.m(j) pk(pm(j),2)], 'w:', 'LineWidth', 1);
        text(ref.l(j), ref.m(j) - 0.08, names{j}, 'Color', 'w', ...
             'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    end
    legend([hr hp], {'Catalog position', 'Detected peak'}, 'Location', 'southwest');
    title(sprintf('%s | err Cyg A: %.4f (%.1f \\lambda/D), Cas A: %.4f (%.1f \\lambda/D)', ...
          strrep(name, '_', '\_'), err(1), err(1)/res, err(2), err(2)/res), ...
          'FontSize', 9);

    savefig(nf, fullfile(outDir, [name '_ref.fig']));
    try
        exportgraphics(nf, fullfile(outDir, [name '_ref.png']), 'Resolution', 200);
    catch
        saveas(nf, fullfile(outDir, [name '_ref.png']));
    end
    close(nf);

    rows(end+1,:) = {name, err(1), err(1)/res, angErr(1), ...
                           err(2), err(2)/res, angErr(2)}; %#ok<SAGROW>
    fprintf('%-18s  CygA: %.4f  CasA: %.4f   (units of l,m)\n', name, err(1), err(2));
end

%% Summary table
T = cell2table(rows, 'VariableNames', {'figure', ...
    'CygA_err_lm', 'CygA_err_lambdaD', 'CygA_err_deg', ...
    'CasA_err_lm', 'CasA_err_lambdaD', 'CasA_err_deg'});
disp(T);
writetable(T, fullfile(outDir, 'position_errors.csv'));
