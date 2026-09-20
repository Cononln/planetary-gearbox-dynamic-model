function plot_my_3d_waterfall(ax, xdata, layers, offsets, colors, opts)
%PLOT_MY_3D_WATERFALL layered 3D waterfall, axis layout of the thesis
%   xdata  : waveform abscissa (frequency Hz for spectra, time s for time)
%   layers : cell{L} of amplitude vectors
%   offsets: [1xL] layer positions on the parameter axis (depth/angle)
%   colors : [Lx3]
%   opts: mode 'spectrum'|'time'; lineW; view [-37.5 16];
%         xLim; showBase; yTickLabels; yDir
%
%   Axis assignment copied from the reference thesis (fig 3-11..3-14):
%     spectrum mode: X = frequency, running along the bottom-right edge;
%                    Y = parameter (depth/angle) on the bottom-left edge,
%                    minimum value at the front corner, increasing toward
%                    the back (YDir normal); one line per parameter value.
%     time mode:     X = parameter (depth/angle) on the bottom-right edge,
%                    minimum value at the front, increasing toward the
%                    back (XDir normal); Y = time on the bottom-left edge
%                    with the maximum time at the front corner
%                    (YDir reverse); one wall per parameter value.
%   Time walls are the raw waveform traces (thin lines, no fill): the
%   fault impulses and their growth with depth/angle read directly.
%   Spectrum layers are plain lines. Layers with the largest parameter are
%   furthest back and are drawn first.

    if nargin < 6, opts = struct(); end
    if ~isfield(opts, 'mode'),        opts.mode = 'spectrum'; end
    if ~isfield(opts, 'lineW'),       opts.lineW = 0.55; end
    if ~isfield(opts, 'view'),        opts.view = [-37.5 16]; end
    if ~isfield(opts, 'xLim'),        opts.xLim = [min(xdata), max(xdata)]; end
    if ~isfield(opts, 'showBase'),    opts.showBase = false; end
    if ~isfield(opts, 'yTickLabels'), opts.yTickLabels = {}; end
    if ~isfield(opts, 'timeTicks'),   opts.timeTicks = []; end
    if ~isfield(opts, 'timePointLimit'), opts.timePointLimit = 10000; end
    if ~isfield(opts, 'preserveExtrema'), opts.preserveExtrema = true; end
    if ~isfield(opts, 'fmHz'),        opts.fmHz = 168; end
    if ~isfield(opts, 'yDir')
        if strcmp(opts.mode, 'time')
            opts.yDir = 'normal';    % 0 s at the front, 1 s receding
        else
            opts.yDir = 'normal';    % min parameter at the front corner
        end
    end

    if ~iscell(layers)
        layers = num2cell(layers, 1);
    end
    if iscolumn(offsets)
        offsets = offsets';
    end
    x = xdata(:)';
    L = numel(layers);
    npf = opts.timePointLimit; % display-only decimation; source values unchanged

    % ---- amplitude range from the raw signals (both modes) ----
    zAll = -inf; zMinAll = inf;
    for k = 1:L
        zAll = max(zAll, max(layers{k}));
        zMinAll = min(zMinAll, min(layers{k}));
    end

    hold(ax, 'on');

    % largest parameter is furthest back: draw it first
    [~, ord] = sort(offsets, 'descend');
    for jj = 1:L
        k = ord(jj);
        switch opts.mode
            case 'time'
                % Display-only downsampling.  Keeping the local minimum and
                % maximum of each time bin prevents narrow crack pulses from
                % disappearing in the exported raster figure.
                if numel(x) <= npf
                    idx = 1:numel(x);
                elseif opts.preserveExtrema
                    nBins = max(1, floor(npf/2));
                    edges = round(linspace(1, numel(x)+1, nBins+1));
                    idx = zeros(1, 2*nBins);
                    nKeep = 0;
                    zFull = layers{k}(:)';
                    for ib = 1:nBins
                        lo = edges(ib);
                        hi = min(numel(x), max(lo, edges(ib+1)-1));
                        [~, iMin] = min(zFull(lo:hi));
                        [~, iMax] = max(zFull(lo:hi));
                        pair = unique(sort([lo+iMin-1, lo+iMax-1]), 'stable');
                        idx(nKeep + (1:numel(pair))) = pair;
                        nKeep = nKeep + numel(pair);
                    end
                    idx = idx(1:nKeep);
                else
                    idx = round(linspace(1, numel(x), npf));
                end
                xd = x(idx); xd = xd(:)';
                zd = layers{k}(:)'; zd = zd(idx); zd = zd(:)';
                plot3(ax, offsets(k) + zeros(size(xd)), xd, zd, ...
                      '-', 'Color', colors(k, :), 'LineWidth', 0.4);
            otherwise
                % spectrum line along X = frequency at constant Y = offsets(k)
                plot3(ax, x, offsets(k) + zeros(size(x)), layers{k}(:)', ...
                      '-', 'Color', colors(k, :), 'LineWidth', opts.lineW);
        end
        if opts.showBase
            plot3(ax, [min(x) max(x)], offsets(k) + [0 0], ...
                  zMinAll + [0 0], '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 0.4);
        end
    end

    % ---- axis layout ----
    ySpan = max(offsets) - min(offsets);
    if ySpan == 0, ySpan = 1; end
    pad = 0.06*ySpan;
    if strcmp(opts.mode, 'time')
        zAbs = max(abs([zMinAll, zAll]));
        ax.ZLim = 1.05 * [-zAbs, zAbs];
        ax.XLim = [min(offsets) - pad, max(offsets) + pad];   % X = parameter
        ax.YLim = opts.xLim;                                  % Y = time
    else
        ax.ZLim = [0, 1.08*zAll];
        ax.XLim = opts.xLim;                                  % X = frequency
        ax.YLim = [min(offsets) - pad, max(offsets) + pad];   % Y = parameter
    end
    ax.YDir = opts.yDir;
    xr = diff(ax.XLim); yr = diff(ax.YLim); zr2 = diff(ax.ZLim);
    if strcmp(opts.mode, 'time')
        % Match the clearer pre-transfer-path time figure: the condition
        % axis is long, time recedes with enough depth, and amplitude stays
        % tall enough to show positive and negative transients.
        ax.PlotBoxAspectRatio = [2.15, 0.92, 1.00];
        ax.PlotBoxAspectRatioMode = 'manual';
        ax.DataAspectRatioMode = 'auto';
    else
        yBox = 1.45;
        zBox = 1.05;
        ax.DataAspectRatio = [xr/2.5, yr/yBox, zr2/zBox];
        ax.DataAspectRatioMode = 'manual';
        ax.PlotBoxAspectRatioMode = 'auto';
    end
    view(ax, opts.view);
    ax.Box = 'on';
    ax.XAxis.FontSize = 9; ax.YAxis.FontSize = 9; ax.ZAxis.FontSize = 9;

    % parameter ticks: on X for time mode, on Y for spectrum mode;
    % the time axis only shows its two end values (thesis style)
    if ~isempty(opts.yTickLabels)
        if strcmp(opts.mode, 'time')
            ax.XTick = offsets;
            ax.XTickLabel = opts.yTickLabels;
            if isempty(opts.timeTicks)
                ax.YTick = [min(opts.xLim), max(opts.xLim)];
            else
                ax.YTick = opts.timeTicks;
            end
        else
            ax.YTick = offsets;
            ax.YTickLabel = opts.yTickLabels;
        end
    end
end
