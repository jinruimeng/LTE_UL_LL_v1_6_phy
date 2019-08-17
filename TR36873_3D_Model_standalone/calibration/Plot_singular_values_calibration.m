% Plot smallest and largest singular values of the channel calculated per
% subcarrier
% % Load all results files from folder
foldername = './results';
filelist  = dir(sprintf('%s/*.mat',foldername));

H_matrix    = [];


% extract results from each results file.
for ii = 1:length(filelist)
    filename = filelist(ii).name;
    load(fullfile(foldername,filename),'UEs');
    fprintf(filelist(ii).name); fprintf('\n');
    % Filter UEs that were simulated.
    active_UE_idcs    = ~[UEs(:).deactivate_UE];
    H_mat_0 = [UEs(active_UE_idcs)];
    H_matrix = [H_matrix, H_mat_0];
end

for i = 1:length(H_matrix)
    H_mat(:,:,:,:,:,i)= H_matrix(:,i).H_0_per_subcarrier;
    
end


for i = 1:length(H_matrix)
    H_mat(:,:,:,:,:,i)= H_matrix(:,i).H_0_per_subcarrier;
end



H_mat = permute(H_mat, [1 2 5 6 3 4]);
for h=1:size(H_mat,4)
    for k = 1:size(H_mat,3)
        H_her(:,:,k,h)= H_mat(:,:,k,h)';
    end
end


for h=1:size(H_mat,4)
    for k = 1:size(H_mat,3)
        H_all(:,:,k,h)= H_mat(:,:,k,h)*H_her(:,:,k,h);
    end
end



H = sum(H_all,3)/6;
H=permute(H, [1 2 4 3]);

for h=1:size(H,3)
    Chan_H_eigen(:,h)=eig(H(:,:,h));
end


EIG = real(Chan_H_eigen);

for i = 1:length(EIG)
    largest_sing_value(i) = max(EIG(:,i));
    smallest_sing_value(i) = min(EIG(:,i));
end


% Plot
figure; grid on; [f1, x1] = ecdf(10*log10(largest_sing_value)); plot(x1,f1,'r');
xlabel('Largest singular value (10log10)');ylabel('ECDF');legend('UMa');
figure; grid on; [f2, x2] = ecdf(10*log10(smallest_sing_value)); plot(x2,f2,'r');
xlabel('Smallest singular value (10log10)');ylabel('ECDF');legend('UMa');


