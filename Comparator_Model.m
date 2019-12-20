function [edge_time, y_sel] = Comparator_Model(y)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% real-valued feedback with comparator using direct learning architecture
% the digital circuit extracts the time when feedback signal crosses the
% set value, this value is afterwards given into DPD calculation 

% generate low speed comparator voltages
comp_val = -0.95:0.15:0.95;
N = ceil(length(y)/length(comp_val));
comp_sig = reshape(repmat(comp_val,N,1),length(comp_val)*N,1);
comp_sig = comp_sig(1:length(y));   

% ---------------------------------------------------------------------
% comparator model
fb_sig = real(y);

% compare the output of the PA with the comparison signal
comp_out = real(fb_sig) > comp_sig;

% detect edges on the comparator output
comp_edges = [comp_out(2:end)-comp_out(1:end-1); 0];
edge_ind = find(abs(comp_edges)==1);
num_samples = length(edge_ind)

% calculate exact edge time using linear interpolation
edge_time = zeros(size(edge_ind));
for i = 1:length(edge_ind)
    ei = edge_ind(i);
    edge_time(i) = ei + (fb_sig(ei)-comp_sig(ei))/(fb_sig(ei)-fb_sig(ei+1));
end

y_sel = comp_sig(edge_ind);
% ---------------------------------------------------------------------

return;









% direct learning DPD with one quadrature and comparator
[U, z_sel] = DDR2_Timed_Matrix(z, state.P, state.M, edge_time, 1024);
M = [real(U) -imag(U)];

% use Gauss-Newton method to update coefficients
b = [real(state.coefs); imag(state.coefs)];
b = b - state.mu*((M'*M)\(M'*(real(z_sel)-y_sel)));
state.coefs = b(1:(size(b,1)/2)) + 1j*b((size(b,1)/2)+1:end);
    
%     % ----- calculate selected points for AM/AM and AM/PM -----
%     % plot the AM/AM and AM/PM characteristics with separated measured points
%     % interpolate the imaginary component to get precise value of
%     % the imaginary value at the sample time
%     fb_sig_i = imag(tx_signal_pa);
%     fb_sig_comp = zeros(size(edge_ind));
%     dpd_out_comp = zeros(size(edge_ind));
%     tx_sig_comp = zeros(size(edge_ind));
%     for i = 1:length(edge_ind)
%         ei = edge_ind(i);
%         t = edge_time(i) - edge_ind(i);
%         fb_sig_comp(i) = comp_sig(ei) + 1j*(fb_sig_i(ei)+(fb_sig_i(ei+1)-fb_sig_i(ei))*t);
%         dpd_out_comp(i) = dpd_out(ei)+(dpd_out(ei+1)-dpd_out(ei))*t;
%         tx_sig_comp(i) = tx_signal(ei)+(tx_signal(ei+1)-tx_signal(ei))*t;
%     end
%     
%     % plot AM/AM characteristisc
%     figure(6);
%     plot(abs(dpd_out), abs(tx_signal_pa), 'o',  'DisplayName','Real PA');
%     hold on;
%     plot(abs(dpd_out_comp), abs(fb_sig_comp), 'x',  'DisplayName','Selected samples');
%     plot(abs(tx_signal), abs(dpd_out), '.',  'DisplayName','DPD output');
%     plot(abs(tx_signal), abs(tx_signal_pa), '.',  'DisplayName','Linearised');
%     hold off;
%     title('AM/AM characteristics');
%     legend('show', 'Location','southeast');
%        
%     % calculate NMSE
%     tx_signal_pan = lscov(tx_signal_pa, tx_signal) * tx_signal_pa;
%     nmserr = nmse(tx_signal, tx_signal_pan);


