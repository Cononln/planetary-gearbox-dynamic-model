function apply_my_figure_style(fig, ax, kind)
%APPLY_MY_FIGURE_STYLE unified global style for all user-gearbox figures
%   kind: '2d' -> 145.2499 x 63.928 mm ; '3d' -> 145.2499 x 85 mm
%   fonts: Times New Roman (latin/digits), SimSun (CJK, auto)
%   sizes: tick 9pt, label 10pt, legend 9pt ; axis line 0.8
%   white background, no title inside axes (title set by layout)
%   Every axes in the figure is styled (subplots/insets included), and
%   legend text also gets the CJK-aware font (fixes hollow-box glyphs).

    if nargin < 3
        kind = '2d';
    end
    switch lower(kind)
        case '2d'
            wmm = 145.2499; hmm = 63.928;
        case '3d'
            wmm = 145.2499; hmm = 85;
        otherwise
            error('kind must be 2d or 3d');
    end
    fig.Units = 'centimeters';
    fig.Position = [2 2 wmm/10 hmm/10];
    fig.Color = 'w';

    axs = findobj(fig, 'Type', 'axes');
    for k = 1:numel(axs)
        a = axs(k);
        a.FontName = 'Times New Roman';
        a.FontSize = 9;
        a.LineWidth = 0.8;
        a.TickDir = 'in';
        a.TickLength = [0.012 0.012];
        a.Box = 'on';
        a.XGrid = 'off'; a.YGrid = 'off'; a.ZGrid = 'off';
        labs = {'XLabel', 'YLabel', 'ZLabel', 'Title'};
        for j = 1:numel(labs)
            h = a.(labs{j});
            h.FontSize = 10;
            h.FontName = pick_font(h.String);
        end
    end
    if strcmpi(kind, '2d')
        ax.Title.String = '';
    end

    lg = findobj(fig, 'Type', 'legend');
    for k = 1:numel(lg)
        lg(k).FontSize = 9;
        lg(k).LineWidth = 0.8;
        lg(k).Box = 'on';
        lg(k).FontName = pick_font(strjoin(cellstr(lg(k).String), ' '));
    end
end
