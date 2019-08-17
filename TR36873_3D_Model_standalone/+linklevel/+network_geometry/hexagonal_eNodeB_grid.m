function [ eNodeB_positions ] = hexagonal_eNodeB_grid(inter_BS_distance, nr_BS_rings)
% Creates an hexagonal set of eNodeBs according to the loaded parameters
% BTS generation is taken from the xf_generate_hsdpa_network from Martin Wrulich's MIMO-HSDPA simulator
% (c) Josep Colom Ikuno, Martin Wrulich INTHFT, 2008

% total nr. of eNodeBs
% each ring consists of 6 NodeBs at the corners of the hexagon
% each edge then has "i-1" further NodeBs, where "i" is the ring index
% total_nr_eNodeB = sum(6*(1:n_rings))+1;

[tmp_gridx,tmp_gridy] = meshgrid(-nr_BS_rings:nr_BS_rings,...
    (-nr_BS_rings:nr_BS_rings)*sin(pi/3)); % regular grid
if mod(nr_BS_rings,2) == 0
    tmp_shift_idx = 2:2:2*nr_BS_rings+1; % shift all even rows
else
    tmp_shift_idx = 1:2:2*nr_BS_rings+1; % shift all odd rows
end

tmp_gridx(tmp_shift_idx,:) = tmp_gridx(tmp_shift_idx,:) + 0.5; % shift

rot = @(w_) [cos(w_),-sin(w_);sin(w_),cos(w_)]; % rotation operator
tmp_hex = zeros(7, 2);
for i_ = 1:7
    % border of the network
    tmp_hex(i_,:) = ((nr_BS_rings+0.5)*rot(pi/3)^(i_-1)*[1;0]).'; % #ok<AGROW>
end

tmp_valid_positions = inpolygon(tmp_gridx,tmp_gridy,tmp_hex(:,1),tmp_hex(:,2));
tmp_x = tmp_gridx(tmp_valid_positions);
tmp_y = tmp_gridy(tmp_valid_positions);

eNodeB_positions = [tmp_x*inter_BS_distance, tmp_y*inter_BS_distance];

