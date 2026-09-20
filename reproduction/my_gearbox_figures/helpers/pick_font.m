function s = pick_font(txt)
%PICK_FONT auto font: SimSun for CJK text, Times New Roman otherwise
    if isempty(txt)
        s = 'Times New Roman';
        return
    end
    if any(double(char(txt)) > 127)
        s = 'SimSun';
    else
        s = 'Times New Roman';
    end
end
