function [state, dpd_st, adapt_st] = Control_Panel(state, par, dpd_st, adapt_st, sa_st)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


if isempty(fieldnames(state))

    % control panel figure
    state.ctrl1_fn = 101;
    state.ctrl1_fh = findobj(allchild(groot), 'flat', 'type', 'figure', 'number', state.ctrl1_fn);
    if ~isempty(state.ctrl1_fh); clf(state.ctrl1_fh); end
        
    state.ctrl1_fh = figure(state.ctrl1_fn);

    bgcolor = state.ctrl1_fh.Color;

    % add sliders

    state.k_sl = uicontrol('Parent',state.ctrl1_fh, 'Style','slider', ...
        'Units','normal', 'Position',[0.25,0.9,0.7,0.05], ...
        'value', 0.39, 'min', 0, 'max', 2);

    state.l_sl = uicontrol('Parent',state.ctrl1_fh, 'Style','slider', ...
        'Units','normal', 'Position',[0.25,0.8,0.7,0.05], ...
        'value', 1, 'min', 0, 'max', 5);

    state.d_sl = uicontrol('Parent',state.ctrl1_fh, 'Style','slider', ...
        'Units','normal', 'Position',[0.25,0.7,0.7,0.05], ...
        'value', 0.39, 'min', 0, 'max', 2);

    state.f_sl = uicontrol('Parent',state.ctrl1_fh, 'Style','slider', ...
        'Units','normal', 'Position',[0.25,0.6,0.7,0.05], ...
        'value', 18, 'min', 10, 'max', 30); % 2.3

    state.k_tx = uicontrol('Parent',state.ctrl1_fh, 'Style','text', ...
        'Units','normal', 'Position',[0,0.9,0.25,0.05],...
        'String', 'k','BackgroundColor',bgcolor);

    state.l_tx = uicontrol('Parent',state.ctrl1_fh, 'Style','text', ...
        'Units','normal', 'Position',[0,0.8,0.25,0.05],...
        'String', 'k','BackgroundColor',bgcolor);

    state.d_tx = uicontrol('Parent',state.ctrl1_fh, 'Style','text', ...
        'Units','normal', 'Position',[0,0.7,0.25,0.05],...
        'String', 'k','BackgroundColor',bgcolor);

    state.f_tx = uicontrol('Parent',state.ctrl1_fh, 'Style','text', ...
        'Units','normal', 'Position',[0,0.6,0.25,0.05],...
        'String', 'k','BackgroundColor',bgcolor);

    % add monitors
    state.mchpwr_sl = uicontrol('Parent',state.ctrl1_fh, 'Style','slider', ...
        'Units','normal', 'Position',[0.25,0.5,0.7,0.05], ...
        'value', -30, 'min', -30, 'max', 0);

    state.dpd_sl = uicontrol('Parent',state.ctrl1_fh, 'Style','slider', ...
        'Units','normal', 'Position',[0.25,0.4,0.7,0.05], ...
        'value', 0, 'min', 0, 'max', 2);

    state.dpd2_sl = uicontrol('Parent',state.ctrl1_fh, 'Style','slider', ...
    'Units','normal', 'Position',[0.25,0.3,0.7,0.05], ...
    'value', 0, 'min', 0, 'max', 2);

    state.adapt1_sl = uicontrol('Parent',state.ctrl1_fh, 'Style','slider', ...
        'Units','normal', 'Position',[0.25,0.2,0.7,0.05], ...
        'value', 0, 'min', 0, 'max', 2);

    state.adapt2_sl = uicontrol('Parent',state.ctrl1_fh, 'Style','slider', ...
        'Units','normal', 'Position',[0.25,0.1,0.7,0.05], ...
        'value', 0, 'min', 0, 'max', 2);

    state.mchpwr_tx = uicontrol('Parent',state.ctrl1_fh, 'Style','text', ...
        'Units','normal', 'Position',[0,0.5,0.25,0.05],...
        'String', 'k','BackgroundColor',bgcolor);

    state.dpd_tx = uicontrol('Parent',state.ctrl1_fh, 'Style','text', ...
        'Units','normal', 'Position',[0,0.4,0.25,0.05],...
        'String', 'k','BackgroundColor',bgcolor);

    state.dpd2_tx = uicontrol('Parent',state.ctrl1_fh, 'Style','text', ...
        'Units','normal', 'Position',[0,0.3,0.25,0.05],...
        'String', 'k','BackgroundColor',bgcolor);

    state.adapt1_tx = uicontrol('Parent',state.ctrl1_fh, 'Style','text', ...
        'Units','normal', 'Position',[0,0.2,0.25,0.05],...
        'String', 'k','BackgroundColor',bgcolor);

    state.adapt2_tx = uicontrol('Parent',state.ctrl1_fh, 'Style','text', ...
        'Units','normal', 'Position',[0,0.1,0.25,0.05],...
        'String', 'k','BackgroundColor',bgcolor);

end


% update texts
state.k_tx.String = sprintf('k=%.2f', state.k_sl.Value);
state.l_tx.String = sprintf('l=%.2f', state.l_sl.Value);
state.d_tx.String = sprintf('d=%.2f', state.d_sl.Value);
state.f_tx.String = sprintf('f=%.2f', state.f_sl.Value);

% update values for DPD and adaptation
for i = 1:length(par)
    dpd_st{i}.in_scaler = state.k_sl.Value;
    dpd_st{i}.out_scaler = state.l_sl.Value;
    adapt_st{i}.in_scaler = state.d_sl.Value;
    adapt_st{i}.out_scaler = state.f_sl.Value;
end

% update monitors
mchpwr = sa_st.mchpwr(1,sa_st.iter);
z_max = dpd_st{1}.z_max;
x_max = dpd_st{1}.x_max;
u_max = adapt_st{1}.z_max;
y_max = adapt_st{1}.y_max;

state.mchpwr_tx.String = sprintf('mchpwr=%.2f', mchpwr);
state.dpd_tx.String = sprintf('z_max=%.2f', z_max);
state.dpd2_tx.String = sprintf('x_max=%.2f', x_max);
state.adapt1_tx.String = sprintf('z_max=%.2f', u_max);
state.adapt2_tx.String = sprintf('y_max=%.2f', y_max);

state.mchpwr_sl.Value = mchpwr;
state.dpd_sl.Value = z_max;
state.dpd2_sl.Value = x_max;
state.adapt1_sl.Value = u_max;
state.adapt2_sl.Value = y_max;

