function [C_clus_mat,G_clus_mat,common_tri] = clustering(vpp_tri_cur,ndv_tri_cur,grains,cell)
% vpp_tri_cur: Nvpp by 3 matrix, with grain ids for each triangle
% ndv_tri_cur: normalized relative velocity for each edges (mean is 1 for
% each layer)
% grains: grains
% C_clus_mat: [no_grains, center_of_gravity, min_y, max_y] - for each
% cluster_id
% G_clus_mat: [cluster_id, grain_x, grain_y, cluster_min_y, cluster_max_y,
% grain_cluster_ell1,grain_cluster_ell2] - for each grain_id
% grain_cluster_ell = 2*min(cluster_max_y-grain_y+0.5,
% grain_y-cluster_min_y+0.5)

% L=5;

% x and y shall be corrected
cgys = grains(:,3);
cgxs = grains(:,2);
ds_per_tri = max(ndv_tri_cur,[],2);
nos_grains = size(grains,1);
H = cell(1,3);
%                     common_tri = vpp_tri_cur(ds_per_tri < quantile(ds_per_tri(:),0.2),:);
% main condition#####
common_tri = vpp_tri_cur(ds_per_tri < 1,:);

static_tri = 1:size(common_tri,1);
g_clus = zeros(1,size(common_tri,1));
cur_c_no = 1;
cluster_finish = 0;
% possible_merge = {};
while (~cluster_finish)
    cur_loc = find(g_clus==0,1); % clus lt 1
    g_clus(cur_loc) = cur_c_no; % here 
    tri_to_check = cur_loc;
%     possible_merge{cur_c_no} = [];

    % reference grain for current cluster is the 1st grain of cur loc tri

    while (~isempty(tri_to_check))
        new_tri_to_check = [];
        for tri_loc = tri_to_check
%                 g_clus(tri_loc) = cur_c_no; % not here

            for tedg = 1
                gA = common_tri(tri_loc,tedg);
                gB = common_tri(tri_loc,mod((tedg+1)-1,3)+1);
                gC = common_tri(tri_loc,mod((tedg+2)-1,3)+1);              
                
                
                com_sum = sum(common_tri==gA,2) + sum(common_tri==gB,2) + sum(common_tri==gC,2);
                % if com_sum == 1, add that tri to possible merge
%                     poss_merg_tri_loc = setdiff(find(com_sum==1),find(g_clus));
%                     possible_merge{cur_c_no} = [possible_merge{cur_c_no},poss_merg_tri_loc(:)'];
                % >=2 for common edg, >=1 for common grain
                new_tri_loc = setdiff(find(com_sum>=1),find(g_clus)); 
                
                % correct these new_tri_locs
                new_tri_loc_grns = common_tri(new_tri_loc,:);
                g_to_correct = union(new_tri_loc_grns(:),[gB,gC]);
                [dxBcorr,dyBcorr] = gBcorr(gA,g_to_correct,grains,cell);
                cgxs(g_to_correct) = cgxs(gA) + dxBcorr;
                cgys(g_to_correct) = cgys(gA) + dyBcorr;
                
                
                
                if ~isempty(new_tri_loc)
                    new_tri_to_check = union(new_tri_to_check, new_tri_loc);
                    g_clus(new_tri_loc) = cur_c_no; % here
                end
            end
        end
        tri_to_check = new_tri_to_check(:)';
    end

    cur_c_no = cur_c_no + 1;

    if (sum(g_clus==0)==0) 
        cluster_finish=1;
    end
end

if isempty(g_clus)
    no_clus = 0;
else
    no_clus = max(g_clus);
end
C_clus_mat = zeros(no_clus,4);
% C_clus_mat: [no_grains, center_of_gravity, min_y, max_y]
G_clus_mat = ones(nos_grains,7);
% G_clus_mat: [cluster_id, grain_x, grain_y, cluster_min_y, cluster_max_y, grain_cluster_ell1, grain_cluster_ell2]
grains_done = [];
for l = 1:no_clus
    cur_clus_grns_tri = common_tri(g_clus==l,:);
    cur_clus_grns = unique(cur_clus_grns_tri(:));
    num_cur_grns = length(cur_clus_grns);
    cur_grains_ys = cgys(cur_clus_grns);
    cur_grains_xs = cgxs(cur_clus_grns);
    
    % figure;plot(grains(:,2),grains(:,3),'ok');hold on;plot(cur_grains_xs,cur_grains_ys,'+b');
    
    if isempty(intersect(grains_done,cur_clus_grns))
        grains_done = union(grains_done,cur_clus_grns);
    else
        disp('setdiff(grains_done,cur_clus_grns): ')
        disp(setdiff(grains_done,cur_clus_grns));
        error('non-empty setdiff(grains_done,cur_clus_grns)');
        
    end
    
    no_grains = length(cur_clus_grns);
    center_of_gravity = mean(cur_grains_ys);
    c_min_y = min(cur_grains_ys)-0.5;
    c_max_y = max(cur_grains_ys)+0.5;
    min_y = zeros(num_cur_grns,1);
    max_y = zeros(num_cur_grns,1);
    for c_grn = 1:num_cur_grns
        c_grn_x = cur_grains_xs(c_grn);
%         if c_grn_x > 1 && c_grn_x < L-1
            c_selec_grns = (cur_grains_xs > c_grn_x-1) & (cur_grains_xs < c_grn_x+1);
%         elseif c_grn_x < 1
%             c_selec_grns = (cur_grains_xs > L+c_grn_x-1) & (cur_grains_xs > 0) & (cur_grains_xs < c_grn_x+1);
%         elseif c_grn_x > L-1
%             c_selec_grns = (cur_grains_xs > c_grn_x-1) & (cur_grains_xs < L) & (cur_grains_xs < c_grn_x+1-L);
%         end
        min_y(c_selec_grns) = min(cur_grains_ys(c_selec_grns))-0.5;
        max_y(c_selec_grns) = max(cur_grains_ys(c_selec_grns))+0.5;
    end
    
    
    ell1 = max_y-cur_grains_ys(:);
    ell2 = cur_grains_ys(:)-min_y;
    grain_cluster_ell1 = 2.*min([ell1,ell2],[],2);
    grain_cluster_ell2 = ell1+ell2;
    grain_cluster_ell2(grain_cluster_ell2>H/2) = H/2;
    
    C_clus_mat(l,:) = [no_grains, center_of_gravity, c_min_y, c_max_y];
    G_clus_mat(cur_clus_grns,:) = [repmat(l,num_cur_grns,1),...
        cur_grains_xs,cur_grains_ys,...
        min_y,max_y,grain_cluster_ell1,grain_cluster_ell2];
    
end

% single grain clusters:
remaining_grains = setdiff(1:nos_grains,grains_done);
next_clus_no = no_clus+1;
for cur_grn = remaining_grains
    cur_clus_grns = cur_grn;
    num_cur_grns = length(cur_clus_grns);
    cur_grains_ys = cgys(cur_clus_grns);
    cur_grains_xs = cgxs(cur_clus_grns);
    
    C_clus_mat(next_clus_no,:) = [1, cur_grains_ys, cur_grains_ys-0.5, cur_grains_ys+0.5];
    G_clus_mat(cur_clus_grns,:) = [repmat(next_clus_no,num_cur_grns,1),...
        cur_grains_xs,cur_grains_ys,...
        cur_grains_ys-0.5,cur_grains_ys+0.5,1,1];
    next_clus_no = next_clus_no + 1;
    
end


end

function [dxBcorr,dyBcorr] = gBcorr(gA,gB,grains,cell)
    gAx = grains(gA,2);
    gAy = grains(gA,3);
    dx = grains(gB,2)-gAx;
    dy = grains(gB,3)-gAy;
    H = cell(1,3);
    L = cell(1,2);
    shift = cell(1,7);
%     shift = 0;
    dyBcorr1 = (dy<-H/2).*H - (dy>H/2).*H;
    dxBcorr1 = (dy<-H/2).*shift - (dy>H/2).*shift;
    dx1 = dx + dxBcorr1;
    dxBcorr2 = (dx1<-L/2).*L - (dx1>L/2).*L;
    dxBcorr = dx+dxBcorr1 + dxBcorr2;
    dyBcorr = dy+dyBcorr1;
   
end
    
    


    
    
    