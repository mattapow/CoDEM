classdef Cluster < PostP
    %Cluster Identify clusters
    
    properties
        
        triplets = [NaN NaN NaN]
        idx = NaN % each row is the cluster index of that particle
        
        rCenters = NaN % Distance between two arbitrary points
        correl = NaN % Probability of grains at distance r being in the same cluster
        
        lTriplet = NaN % Length scale from triplet jamming
        lVf = NaN % length scale from velocity fluctuations
        lContact = NaN % Length scale from contact lifetime / shear rate
        
        k = NaN
    end
    
    methods
        
        function obj = Cluster(dirPath)
            obj = obj@PostP(dirPath);
        end
        
        function l = getL(obj)
            % L for system length/size
            l = obj.L;
        end
        
        %% Velocity Fluctuations
        function obj = getLVf(obj)
            % Get the average length scale at all saved timesteps
            l = zeros(obj.tEnd-2, 1);
            for filenum = 2:obj.tEnd-1
                l(filenum-1) = getLVf1(obj, filenum);
            end
            obj.lVf = mean(l);
        end
        
        function l = getLVf1(obj, filenum)
            % Get the length scale based on velocity fluctuations in one
            % file
            % \deltav = \omega l
            
            % Read grain position and velocity
            fp = strcat(obj.dirPath, '/grain/grain_', string(filenum));
            formatspec = '%d %f %f %f %f %*f %*f %*f %*f %*f %*f';
            grain = Grain();
            grain = grain.readGrain(fp, formatspec);            
            
            % Read shear rate
            fp = strcat(obj.dirPath, '/para/parameter_', string(filenum));
            formatspec = '%*f %*f %*f %*f %*f %f %f';
            para = Para();
            para = para.readPara(fp, formatspec);
            
            % Read grain velocity gradient based on vonoroi neighbours
            fp = strcat(obj.dirPath, '/grain_post_process/grain_pp', string(filenum));
            formatspec = '%*f %*f %*f %*f %f %f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f';
            grainPp = GrainPp();
            grainPp = grainPp.readGrainPp(fp, formatspec);
            
            % <dv^2>, dv = v_x - \dot{\gamma}*y for each grain
            dV2 = mean((grain.V(:,1) - para.shear_rate*grain.X(:,2)).^2);
            
            % <dw^2>, dw = w - \dot{\gamma}/2 for each grain
            % NB: w = dv/dy relies on vonoroi neighbours to get dy
            gradV = reshape(grainPp.gradV,2,2,[]);
            w = zeros(size(gradV,3), 1);
            for i = 1:size(gradV, 3)
                w(i) = norm(.5*(gradV(:,:,i) + gradV(:,:,i)'));
            end
            dw2 = mean((w - para.shear_rate/2).^2);
            
            l = sqrt(dV2/dw2);
            
        end
        
        %% Triplet Jamming
        function Nbrs = getNbrs(obj, filenum)
            % Nbrs are the grain vonoroi neighbours as saved in grainpp
            
            % read grain neighbours
            fp = strcat(obj.dirPath, '/grain_post_process/grain_pp', string(filenum));
            formatspec = '%*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %f %f %f %f %f %f %f %f %f %f %f %f';
            grainPp = GrainPp();
            grainPp = grainPp.readGrainPp(fp, formatspec);
            Nbrs = grainPp.voro_neighbour_ID;
            
        end
        
        function obj = getTriplets(obj, filenum)
            % returns num_tri by 3 matrix with each row corresponding to the
            % three neighbouring grains that form a triangle
            
            nbrs = getNbrs(obj, filenum);
            num_grain = max(max(nbrs));
            
            % cur_nbrs = grains_vpp(:,22:33)+1; % grain_loc = grain_id + 1
            n_tris = 2*size(nbrs,1);
            vpp_tri = zeros(n_tris,3);
            vp_count = 1;
            
            % Read all grain positions and radii
            grain = Grain();
            fp = strcat(obj.dirPath, '/grain/grain_', string(filenum));
            formatspec = '%d %f %f %*f %*f %*f %*f %*f %*f %*f %f';
            grain = grain.readGrain(fp, formatspec);
            
            for g1=1:num_grain
                g1Nbrs = nbrs(g1,nbrs(g1,:)>g1); % no need to consider any previous grain nos
                for i = 1:length(g1Nbrs)
                    g2 = g1Nbrs(i);
                    g2Nbrs = nbrs(g2,nbrs(g2,:)>g1);% no need to consider any previous grain nos
                    % select third common grain, no need to consider previous g1Nbrs.
                    gg3 = setdiff(intersect(g1Nbrs,g2Nbrs),g1Nbrs(1:(i-1)));
                    for g3 = gg3
                        % if the three grains each touch at least one other
                        if (Cluster.tripletTouch(grain, g1, g2, g3))
                            vpp_tri(vp_count,:) = sort([g1,g2,g3]);
                            vp_count = vp_count + 1;
                        end
                    end
                end
            end
            
            vpp_tri(sum(vpp_tri,2)==0,:)=[];
            obj.triplets = vpp_tri;
            
        end
        
        function ds = getDs(obj, triplets, filenum)
            % Compute all triplet inner displacement at time 'time'
            % triplets is a nTripletx3 list of triplets
            % ds_row = sum of the distances between grains iwithin a triplet
            
            % Read grain positions
            fp = strcat(obj.dirPath, '/grain/grain_', string(filenum));
            formatspec = '%d %f %f %*f %*f %*f %*f %*f %*f %*f %*f';
            grain = Grain();
            grain = grain.readGrain(fp, formatspec);
            X = grain.X;
            
            % Reads time, L shift
            fp = strcat(obj.dirPath, '/cell/cell_', string(filenum));
            formatspec = '%f %f %f %*f %*f %*f %f %*f %*f %*f %*f %*f %*f %*f';
            cell = Cell();
            cell = cell.readCell(fp, formatspec);
            
            % Grain positions
            x = X(:,1);
            y = X(:,2);
            
            % Inter-grain distance of within triplets
            dx = x(triplets)-x(triplets(:,[2,3,1]));
            dy = y(triplets)-y(triplets(:,[2,3,1]));
            [dx, dy] = cell.CLP(dx, dy);
            ds = hypot(dx, dy);
            
        end
        
        function t = getT(obj, filenum)
            % Read cell time
            fp = strcat(obj.dirPath, '/cell/cell_', string(filenum));
            formatspec = '%f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f';
            cell = Cell();
            cell = cell.readCell(fp, formatspec);
            t = cell.time;
        end
        
        function ndVTri = getDVTri(obj, filenum)
            % Compute the triplet inner speeds
            % Based on central difference in time of triplet displacements
            
            if (filenum == 1 || filenum == obj.tEnd); error('Cannot compute velocity at edge cases.'); end
            if (any(isnan(obj.triplets))); error('Must compute triplets first.'); end
            
            % Compute the sum of the distance within triplets
            ds_prev = getDs(obj, obj.triplets, filenum-1);
            ds_next = getDs(obj, obj.triplets, filenum+1);
            
            % Get the times at both before + after snapshots
            t_prev = getT(obj, filenum-1);
            t_next = getT(obj, filenum+1);
            dt = t_next - t_prev;
            
            % Compute the speed within the triplets
            dV = abs(ds_next - ds_prev) / dt;
            
            % Normalise the speed
            %             ndVTri = dV / mean(dV);
            ndVTri = dV / quantile(dV, .5);
            %             ndVTri = dV;
        end
        
        %         function [obj, G_clus_mat, ndVTri] = getClusters(obj, filenum)
        function [obj, G_clus_mat] = getClusters(obj, filenum)
            
            % Get the triplets
            obj = getTriplets(obj, filenum);
            
            % Get the inner triplet speeds
            ndVTri = getDVTri(obj, filenum);
            
            % Read current grains and cell
            grains_cur = load(sprintf('%s/grain/grain_%d',obj.dirPath,filenum));
            cell_cur = load(sprintf('%s/cell/cell_%d',obj.dirPath,filenum));
            
            % Compute the clusters
            [C_clus_mat,G_clus_mat, ~] = clustering(obj.triplets, ndVTri, grains_cur, cell_cur);
            obj.idx = G_clus_mat(:, 1);
            
        end
        
        function obj = clearBig(obj)
            % Clear triplets and idx to reduce the amount of data saved
            obj.triplets = [NaN NaN NaN];
            obj.idx = NaN;
        end
        
        function obj = getLTriplet(obj)
            % Get the correlation between distance between two particles and whether
            % they are in the same cluster or not
            
            % Presize pdf sums
            dr = 1;
            obj.rCenters = 0:dr:100;
            edges = [obj.rCenters-dr, max(obj.rCenters)+dr];
            sumHSameClus = zeros(1,length(obj.rCenters));
            sumHAll = zeros(1,length(obj.rCenters));
            
            % Sample some number of timesteps
            nSample = 100; % 200
            fprintf('Number of sample timesteps: %d\n', nSample);
            tSample = randsample(2:obj.tEnd-1,nSample);
            %             progress = 0;
            
            for filenum = tSample
                %                 fprintf('%.1f%%\n', progress/nSample*100);
                %                 progress = progress + 1;
                
                % Find the clusters
                clust=Cluster(obj.dirPath);
                [~, clustG] = getClusters(clust, filenum);
                
                % Readability - cluster ID, X, Y for each grain
                clustID = clustG(:,1);
                grainX = clustG(:,2);
                grainY = clustG(:,3);
                L = clust.L(1);
                H = clust.L(2);
                shift = clust.getShift(filenum);
                
                % Repmat over Lees-Edwards BC
                clustGBig = [[clustID,     grainX,       grainY];...
                    [clustID+1e7, grainX+L,     grainY];...
                    [clustID+2e7, grainX-L,     grainY];...
                    [clustID+3e7, grainX+shift, grainY+H];...
                    [clustID+4e7, grainX-shift, grainY-H]];
                
                % Random selection of pairs of grains
                nSelect = 1e7;
                nGrains = size(clustGBig,1);
                pair1 = randsample(nGrains,nSelect,true);
                pair2 = randsample(nGrains,nSelect,true);
                
                % Distance between selected pairs
                dx = clustGBig(pair2,2)-clustGBig(pair1,2);
                dy = clustGBig(pair2,3)-clustGBig(pair1,3);
                r = hypot(dx, dy);
                
                % Get grain pairs in the same cluster
                sameClusGrain = clustGBig(pair1,1)==clustGBig(pair2,1);
                sameClus = sameClusGrain & (r<=L);
                
                % Get the pdf (counts) in each rBin for: every pair, and pairs in the same cluster
                [hAll,~] = histcounts(r,edges);
                [hSameClus,~] = histcounts(r(sameClus),edges);
                
                % Add these counts to the sum over timesteps
                sumHSameClus = sumHSameClus + hSameClus;
                sumHAll = sumHAll + hAll;
                
                %         title(max(r.*sameClus))
                %         drawnow
            end
            
            % Correlation of between distance between two particles and whether they are in the same cluster or not
            obj.correl = sumHSameClus./sumHAll;
            obj.correl(1) = 1;
            obj.correl(isnan(obj.correl)) = 0; % remove nan
            
            %             semilogy(rBins,correl)
            %             hold on
            %             xlim([0,30])
            %             drawnow
            
            % Take (pow)th moment
            pow = 1;
            obj.lTriplet= (sum((obj.rCenters+1).^pow.*obj.correl)./sum(obj.correl))^(1/pow);
            
            %             csvwrite('rVsCor.csv',[rCenters',correl']);
        end
        
        %% Kmean clustering
        function plotK(obj)
            hold on
            for kk = 5:100
                [~, ~, sumd] = kMeansOne(obj, kk);
                scatter(kk, sumd+kk)
            end
            hold off
            xlabel('k')
            ylabel('Inner cluster distance')
        end
        
        function [obj, sumd] = kMeansOne(obj, k, filenum)
            % obj = kMeans(obj, t)
            % cluster the grains at time t into k clusters according to
            % grain: position, vorticity (gradV.asym),?shear_rate (gradV.sym)
            
            % Read grain gradV, V
            fpGrainPp = strcat(obj.dirPath, '/grain_post_process/grain_pp', string(filenum));
            formatspec =  '%*f %*f %*f %*f %f %f %f %f %f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f';
            grainPp = GrainPp();
            grainPp = grainPp.readGrainPp(fpGrainPp, formatspec);
            dvdy = grainPp.gradV;
            %             V = grainPp.V;
            %             strain = (dvdy(:,2) + dvdy(:,3))/2;
            %             vorticity = (dvdy(:,2) - dvdy(:,3))/2;
            
            % Read grain positions
            fpGrain = strcat(obj.dirPath, '/grain/grain_', string(filenum));
            formatspec = '%d %f %f %*f %*f %*f %*f %*f %*f %*f %*f';
            grain = Grain();
            grain = grain.readGrain(fpGrain, formatspec);
            X = grain.X;
            
            % cluster them
            data = [dvdy X];
            [idx,~,sumd] = kmeans(data,k);            %#ok<PROPLC>
            sumd = norm(sumd);
            obj.idx = idx; %#ok<PROPLC>
            obj.k = k;
        end
        
        %% Visu
        function visuClust(obj, filenum)
            % Visulise grains clusters at a single timesteps
            
            % Filepaths
            fpGrain = strcat(obj.dirPath, '/grain/grain_', string(filenum));
            fpCell = strcat(obj.dirPath, '/cell/cell_', string(filenum));
            fpPara = strcat(obj.dirPath, '/para/parameter_', string(filenum));
            
            % Read in Grain positions and radii
            formatspec = '%d %f %f %*f %*f %*f %*f %*f %*f %*f %f';
            grain = Grain();
            grain = grain.readGrain(fpGrain, formatspec);
            
            % Plot the clustered grains
            if isnan(obj.idx); error('Must run getCluster for this instance.'); end
            grain.plotGrainClust(obj.idx);
            
            % Plot the delaunay triangulation
            %             hold on; grain.plotVonoroi; hold off
            
            
            % Read in Cell time, length and shift
            formatspec = '%f %f %f %*f %*f %*f %f %*f %*f %*f %*f %*f %*f %*f';
            cell = Cell();
            cell = cell.readCell(fpCell, formatspec);
            
            % Read in the confining Pressure, shear Rate and cohesion
            % parameters
            formatspec = '%*f %*f %*f %*f %f %f %f';
            para = Para();
            para = para.readPara(fpPara, formatspec);
            
            ax = gca;
            ax.XLim = [0-1 cell.L(1)+1];
            ax.YLim = [-cell.L(2)/2-1 +cell.L(2)/2+1];
            ttl = sprintf('P = %.0f, Shear Rate = %.2f\nC = %.0f, nShears = %.0f', para.P, para.shear_rate, para.cohesion, cell.time*para.shear_rate);
            title(ttl)
            
            % Plot the cell boundaries
            cell.plotCell();
            
            axis equal
            axis tight
            
        end
        
    end
    
    methods (Static)
        function output = tripletTouch(grain, g1, g2, g3)
            % returns true if the three grains each touch at least one other
            output = true;
            
            % Get triplet grain positions and radii
            r(3) = grain.R(g3);
            r(2) = grain.R(g2);
            r(1) = grain.R(g1);
            X(3,:) = grain.X(g3,:);
            X(2,:) = grain.X(g2,:);
            X(1,:) = grain.X(g1,:);
            
            % for each grain in triplet
            for i = 1:3
                % Get the other two grains
                Nbr1 = mod(i+1, 3)+1;
                Nbr2 = mod(i+3, 3)+1;
                
                % Determine whether this grain overlaps its neighbours
                touch1 = r(i)+r(Nbr1) >= norm(X(Nbr1,:)-X(i,:));
                touch2 = r(i)+r(Nbr2) >= norm(X(Nbr2,:)-X(i,:));
                
                % if this grain doesn't overlap either neighbour, return false
                if ~(touch1 || touch2)
                    output = false;
                    return;
                end
            end
            
        end
    end
    
    
end

