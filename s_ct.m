% S_CT - Starts the s_ct User Interface
%
%   Author:         Gianluca Iori (gianthk.iori@gmail.com)
%   BSRT - Charite Berlin
%   Created on:   30/01/2018
%   Last update:  30/01/2018
%
%   this class is part of the synchro toolbox    
%   ______________________________________________________

% generate valid variable name
i = 1;
while exist(['sctdata_' num2str(i)]) == 1
    i = i+1;
end
eval([genvarname(['sctdata_' num2str(i)]) ' = sctdata;']);

% launch the s_ctGUI
eval(['s_ctGUI_L(sctdata_' num2str(i) ')']);
clear i;

