function [cf, Pos] = whi_fig_set_position(cf, PosStr)
% [cf, Pos] = whi_fig_set_position(cf, PosStr)
%
% Moves and resizes the figure according to: Top, Bottom, Left, Right, All

    ScreenSize = get( 0, 'ScreenSize' );

    if strcmpi(PosStr, 'Top')
        set(cf, 'Position', [ScreenSize(1) ScreenSize(4)/2 ScreenSize(3) ScreenSize(4)/2]);
    elseif strcmpi(PosStr, 'Bottom')
        set(cf, 'Position', [ScreenSize(1) ScreenSize(2) ScreenSize(3) ScreenSize(4)/2]);
    elseif strcmpi(PosStr, 'Left')
        set(cf, 'Position', [ScreenSize(1) ScreenSize(2) ScreenSize(3)/2 ScreenSize(4)]);
    elseif strcmpi(PosStr, 'Right')
        set(cf, 'Position', [ScreenSize(3)/2 ScreenSize(2) ScreenSize(3)/2 ScreenSize(4)]);
    elseif strcmpi(PosStr, 'All')
        set(cf, 'Position', ScreenSize);
    else
        warning('chk:arg', 'Not recognized position');
    end
    
    Pos = get(cf, 'Position');
    

end
