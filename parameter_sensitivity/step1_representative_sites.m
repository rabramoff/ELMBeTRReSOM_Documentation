function nsites = step1_representative_sites(n_ensemble)

% input, 3: how many representative sites ~ n*100
% output, sites infomation is saved in global_sparse_grid.txt

%% read in MAT,MAP,PFT map
% addpath('/Users/QZhu/Desktop/ILAMB_sample/DATA/');

ilamb_pftmask = ncread('ilambgpp_mask.nc','mask');
ilamb_pftmask = [ilamb_pftmask(361:720,:); ilamb_pftmask(1:360,:)];
pftmap = ncread('pft_360x720_2000.nc','PCT_NAT_PFT');
Temp = ncread('tas_0.5x0.5.nc','tas');
Prec = ncread('pr_0.5x0.5.nc','pr');
MAT = NaN(720,360); %R creates a giant matrix of NaN
MAP = NaN(720,360); 
for i = 1 : 720
    for j = 1 : 360
        if ~isnan(sum(Temp(i,j,:))) && ~isnan(sum(Prec(i,j,:)))
            MAT(i,j) = mean(Temp(i,j,:))-273.15;
            MAP(i,j) = mean(mean(reshape(Prec(i,j,:),12,[])))*3600*24*365;
        end
    end
end
MAP(MAP>4500) = 4500; %fixes the MAT and MAP units and fills in the matrix
save  MAT_climatology.mat MAT;
save  MAP_climatology.mat MAP;
clear Temp Prec;
%imagesc(MAT);
%imagesc(MAP);

% Dominate pft definition:
% 1. vegetated surface > 90%
% 2. dominant pft cover > 50%
pft_domi = NaN(720,360);
for i = 1 : 720
    for j  = 1 : 360
        if (max(reshape(pftmap(i,j,2:17),1,[])) > 50 && ilamb_pftmask(i,j) == 1) % pftmap(i,j,1) < 10 &&  
            [M,I] = max(reshape(pftmap(i,j,2:17),1,[]));
            pft_domi(i,j) = I;
        end
    end
end
% figure 1
figure;imagesc(pft_domi>0& ilamb_pftmask== 1);title('vegetated surface > 90%, dominant pft cover > 50%')
pft_domi = [pft_domi(361:720,:); pft_domi(1:360,:)];

burn = nanmean(ncread('burntArea_0.5x0.5.nc','burntArea'),3);

MAT_MAP_PFT = [];
lon_lat = [];
for i = 1 : 720
    for j  = 1 : 360
        if ~isnan(pft_domi(i,j)) && ~isnan(MAT(i,j)) && ~isnan(MAP(i,j))
            MAT_MAP_PFT = [MAT_MAP_PFT; MAT(i,j) MAP(i,j) pft_domi(i,j) burn(i,j)];
            lon_lat = [lon_lat; [i j]];
        end
    end
end
%{
Xedges = [-20:1:30]; % T
Yedges = [0:100:4000]; % P
burn_aggregate = zeros(length(Xedges)-1,length(Yedges)-1);
burn_num = zeros(length(Xedges)-1,length(Yedges)-1);
burn_aggregate_pft = zeros(length(Xedges)-1,length(Yedges)-1,17);
burn_num_pft = zeros(length(Xedges)-1,length(Yedges)-1,17);
for i = 1 : 720
    for j = 1 : 360
        if ~isnan(pft_domi(i,j)) && ~isnan(MAT(i,j)) && ~isnan(MAP(i,j))
            mat_index = min(int32((MAT(i,j)+ 20)) + 1,50);
            map_index = min(int32((MAP(i,j)- mod(MAP(i,j),100))/100)+ 1, 40);
            burn_aggregate(mat_index,map_index)= burn_aggregate(mat_index,map_index) +burn(i,j);
            burn_num(mat_index,map_index)= burn_num(mat_index,map_index) + 1;
            
            burn_aggregate_pft(mat_index,map_index,pft_domi(i,j))=burn_aggregate_pft(mat_index,map_index,pft_domi(i,j))+burn(i,j);
            burn_num_pft(mat_index,map_index,pft_domi(i,j))=burn_num_pft(mat_index,map_index,pft_domi(i,j)) + 1;
        end
    end
end
burn_aggregate= burn_aggregate./burn_num;
burn_aggregate_pft=burn_aggregate_pft./burn_num_pft;
%}
%subplot(1,2,1);
%imagesc(burn_aggregate);hold;
%xlabel('MAP');ylabel('MAP');title('All')
%subplot(1,2,2);imagesc(burn_num);hold;
%xlabel('MAP');ylabel('MAP');title('All')
%for i = [6 13 14 15]
%    figure;
%    subplot(1,2,1);
%    imagesc(burn_aggregate_pft(:,:,i));hold;
%    xlabel('MAP');ylabel('MAT');title(['PFT' num2str(i)])
%    subplot(1,2,2);imagesc(burn_num_pft(:,:,i));hold;
%    xlabel('MAP');ylabel('MAT');title(['PFT' num2str(i)])
%end

%h = histogram2(MAT_MAP_PFT(:,1),MAT_MAP_PFT(:,2),Xedges,Yedges);
%ndhist(MAT_MAP_PFT(:,1:2),'Ctrs', centers,'FaceAlpha',.25);
%xlabel('MAT'); ylabel('MAP');
%set(gcf,'renderer','opengl');
%set(gca,'xdir','reverse');
%set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');

%ndhist(lon_lat,[50, 40],'FaceAlpha',.25);
%xlabel('lon'); ylabel('lat');
%set(gcf,'renderer','opengl');
%set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');

%% for each pft, geographic location map vs MAT/MAP map
pft_count=NaN(16,1);
for i = 1:16
    %figure;
    tmp = lon_lat(MAT_MAP_PFT(:,3) == i,1:2);
    %if ~isempty(tmp)
    %    subplot(1,2,1);
    %    hist3(tmp,[40, 40],'FaceAlpha',.25);
    %    xlabel('lon'); ylabel('lat');
    %    set(gcf,'renderer','opengl');
    %    set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
    %end
    pft_count(i,1) = length(tmp);
    %tmp = MAT_MAP_PFT(MAT_MAP_PFT(:,3) == i,1:2);
    %if ~isempty(tmp)
    %    subplot(1,2,2);
    %    hist3(tmp,[40, 40],'FaceAlpha',.25);
    %    xlabel('MAT'); ylabel('MAP');
    %    set(gcf,'renderer','opengl');
    %    set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
    %end
end

%% pft hist gram
pft_count_pct = pft_count./sum(pft_count)*100;
% figure 2
figure;
bar(pft_count);xlim([0 17]);hold;
for i = 1 : length(pft_count)
    text(i-0.35,pft_count(i)+200,sprintf('%.2f',pft_count_pct(i)),'FontSize',12,'FontWeight','bold');
end
ylabel('Frequency (#)','FontSize',18,'FontWeight','bold');
set(gca,'Xtick',[1:16],'XTickLabel',{'NET Temp','NET Bor','NDT Bor','BET Tro','BET Temp','BDT Tro','BDT Temp','BDT Bor','BES Temp','BDS Temp','BDS Bor','C3 arc','C3','C4','CropR','CropI'})
xtickangle(45)

%% get most representative gridcells
% PFT         PCT      representative
% 1  NET Temp
% 2  NET Bor
% 3  NDT Bor
% 4  BET Tro
% 5  BET Temp
% 6  BDT Tro
% 7  BDT Temp
% 8  BDT Bor
% 9  BES Temp
% 10 BDS Temp
% 11 BDS Bor
% 12 C3 arc
% 13 C3
% 14 C4
% 15 Crop
% total
figure;
repre_site_ori = ceil(pft_count_pct); % BES Temp
for ensemble_case = 1 : n_ensemble
    repre_site = repre_site_ori*ensemble_case;
    repre_MAT = NaN(15,max(repre_site));
    repre_MAP = NaN(15,max(repre_site));
    repre_lon = NaN(15,max(repre_site));
    repre_lat = NaN(15,max(repre_site));
    for i = 1:15
        tmp = MAT_MAP_PFT(MAT_MAP_PFT(:,3) == i,1:2); % get lon lat of each pft class
        centers = {(-20:1:30),(0:100:4000)};
        [n,c] = hist3(tmp,'Ctrs', centers);
        [a b] = find(n==(max(max(n)))); % get the most representative T/P cliamte region
        if length(a) > 1
            a = a(1);
        end
        if length(b) >  1
            b = b(1);
        end
        repre_MAT(i,1) = c{1,1}(a); % get T of the most representative T/P cliamte region
        repre_MAP(i,1) = c{1,2}(b); % get P of the most representative T/P cliamte region
        if repre_site(i) > 1 % next: get secondary most representative T/P climate region, and so on
            for j = 1 : repre_site(i) - 1
                n(a,b) = 0;
                [a b] = find(n==(max(max(n))));
                if length(a) > 1
                    a = a(1);
                end
                if length(b) >  1
                    b = b(1);
                end
                repre_MAT(i,1+j) = c{1,1}(a);
                repre_MAP(i,1+j) = c{1,2}(b);
            end
        end
    end
    
    % calculate the actual land surface covered by representative T/P regions
    repre_site_cumsum = cumsum(repre_site);
    repre_site_cumsum = [0; repre_site_cumsum];
    pft_mask = NaN(720,360);
    pft_climate_mask = NaN(720,360);
    for i = 1 : 15
        for j = repre_site_cumsum(i)+1 : repre_site_cumsum(i+1)
            for l = 1 : 720
                for ll = 1 : 360
                    if abs(MAT(l,ll)-repre_MAT(i,j-repre_site_cumsum(i)))<=0.5 && abs(MAP(l,ll)-repre_MAP(i,j-repre_site_cumsum(i)))<=50
                        if (i == pft_domi(l,ll)) % multiple pft can share the same MAT/MAP climate
                            pft_mask(l,ll) = pft_domi(l,ll);
                            pft_climate_mask(l,ll) = j;
                        end
                    end
                end
            end
        end
    end
    
    % MAT_MAP representative grid -> sparse grid, pick up the nearest grid
    % cells to the center of those representative grid cell cluster
    sparse_grid = [];
    for i = 1 : 15
        for j = repre_site_cumsum(i)+1 : repre_site_cumsum(i+1)
            [a,b] = find(pft_climate_mask==j);
            [av,ai] = min(sum(([a,b] - repmat([median(a),median(b)],length(a),1)).^2,2));
            sparse_grid = [sparse_grid; [a(ai)/2-180.25 b(ai)/2-90.25 a(ai) b(ai) double(i) double(j)]]; %also index -> lon lat
        end
    end
    % sparse grid -> ilamb global mask
    %save('ilamb_mask.mat','pft_climate_mask'); % also see sparse_grid(:,6)
    %movefile('ilamb_mask.mat',['ilamb_mask' num2str(ensemble_case) '.mat'])
    
    %nccreate('ilammask_0.5x0.5_c170903.nc','pftmask','Dimensions', {'lat',360,'lon',720});
    %ncwriteatt('ilammask_0.5x0.5_c170903.nc','pftmask','long_name','representative pft');
    %ncwrite('ilammask_0.5x0.5_c170903.nc','pftmask',pft_climate_mask');
    
    % plot
    addpath('m_map');
    Lonmax=180;
    Lonmin=-180;
    Latmax=90;
    Latmin=-90;
    nlon=720;
    nlat=360;
    tmp=NaN(720,360);
    
    for i = 1 : nlat
        latitude(i,1:nlon) = Latmin + 0.5*double(i-1);
    end
    for i = 1 : nlon
        longitude(1:nlat,i) = Lonmin + 0.5*double(i-1);
    end
    
    subplot(n_ensemble,1,ensemble_case,...
        'FontSize',10,'FontWeight','bold');
    m_proj('Robinson','lon',[Lonmin Lonmax],'lat',[Latmin Latmax]);
    m_coast('line');
    hold;
    m_pcolor(longitude,latitude,pft_climate_mask'); %pft_domi'
    for i = 1 : length(sparse_grid)
        m_text(sparse_grid(i,1)-0.25,sparse_grid(i,2)-0.25,'o');
    end
    colormap('jet');
    shading flat;
    m_grid('box','on','tickdir','in');
    title(sprintf('%d # of representative sites, %.1f%% spatial coverage',sum(repre_site),sum(sum(pft_climate_mask>0))/length(MAT_MAP_PFT)*100),'FontWeight','bold','FontSize',14);
    set(gca, ...
        'Box'         , 'on'     , ...
        'TickDir'     , 'in'     , ...
        'TickLength'  , [.01 .02] , ...
        'XMinorTick'  , 'off'      , ...
        'YMinorTick'  , 'on'      , ...
        'YGrid'       , 'off'      , ...
        'XColor'      , [.01 .01 .01], ...
        'YColor'      , [.01 .01 .01], ...
        'LineWidth'   , 2         );
end

fid = fopen('global_sparse_grid.txt','w');
for i  = 1 : length(sparse_grid)
    fprintf(fid,'%.2f\t%.2f\t%d\t%d\t%d\t%d\n',sparse_grid(i,1),sparse_grid(i,2),sparse_grid(i,3),sparse_grid(i,4),sparse_grid(i,5),sparse_grid(i,6));
end

nsites = length(sparse_grid);

%% iLAMB data for each representative climate-biome gridcell
%{
% split data for representative climate-biome gridcell
repre_site = int32(pft_count_pct);
repre_site(4) = repre_site(4)*2; % BET Tro
repre_site(6) = repre_site(6)*2; % BDT Tro
repre_site(9) = 1; % BES Temp
repre_site(11) = repre_site(11) + 1; % BDS Bor
repre_site(10) = repre_site(10); % BDS Temp
repre_site(14) = repre_site(14); % C4 grass
ensemble_case = 3;
repre_site = repre_site*ensemble_case;
repre_MAT = NaN(15,max(repre_site));
repre_MAP = NaN(15,max(repre_site));
repre_lon = NaN(15,max(repre_site));
repre_lat = NaN(15,max(repre_site));
for i = 1:15
    tmp = MAT_MAP_PFT(MAT_MAP_PFT(:,3) == i,1:2);
    centers = {(-20:1:30),(0:100:4000)};
    [n,c] = hist3(tmp,'Ctrs', centers);
    [a b] = find(n==(max(max(n))));
    if length(a) > 1
        a = a(1);
    end
    if length(b) >  1
        b = b(1);
    end
    repre_MAT(i,1) = c{1,1}(a);
    repre_MAP(i,1) = c{1,2}(b);
    if repre_site(i) > 1
        for j = 1 : repre_site(i) - 1
            n(a,b) = 0;
            [a b] = find(n==(max(max(n))));
            if length(a) > 1
                a = a(1);
            end
            if length(b) >  1
                b = b(1);
            end
            repre_MAT(i,1+j) = c{1,1}(a);
            repre_MAP(i,1+j) = c{1,2}(b);
        end
    end
end
repre_site_cumsum = cumsum(repre_site);
repre_site_cumsum = [0; repre_site_cumsum];
pft_mask = NaN(720,360);
pft_climate_mask = NaN(720,360);
for i = 1 : 15
    for j = repre_site_cumsum(i)+1 : repre_site_cumsum(i+1)
        for l = 1 : 720
            for ll = 1 : 360
                if abs(MAT(l,ll)-repre_MAT(i,j-repre_site_cumsum(i)))<=0.5 && abs(MAP(l,ll)-repre_MAP(i,j-repre_site_cumsum(i)))<=50
                    if (i == pft_domi(l,ll)) % multiple pft can share the same MAT/MAP climate
                        pft_mask(l,ll) = pft_domi(l,ll);
                        pft_climate_mask(l,ll) = j;
                    end
                end
            end
        end
    end
end
%imagesc(pft_climate_mask);

% gpp 82-08, resp 82-08, lai 82-10, et 80-11
start_year = [1982,1982,1982,1980,1,1];
end_year = [2008,2008,2010,2011,1,1];
var_name = {'gpp','nee','lai','et','biomass','soilc'};
for var = 1 : 6
    
    % read in data
    data = NaN(720,360,end_year(var)-start_year(var) + 1);
    myindex = 1;
    if var <=4 % temporally resolved data, gp, nee, lai, et
        for i = start_year(var) : end_year(var)
            for j = 1 : 12
                data(:,:,myindex) = ncread(sprintf('/Users/QZhu/Desktop/data/ILAMB_data/%s/%.4d/%s_0.5x0.5_%.4d%.2d.nc',var_name{var},i,var_name{var},i,j),var_name{var});
                myindex = myindex + 1;
            end
        end
    else % not temporally resolved data, biomass, soilc
        data(:,:) = ncread(sprintf('/Users/QZhu/Desktop/data/ILAMB_data/%s/%s_0.5x0.5.nc',var_name{var},var_name{var}),var_name{var});
    end
    if var == 1 || var == 2 % gpp nee kgC/m2/s
        data = data * 3600*24*1000; % g/m2/day
    elseif var == 4 % et kg/m2/s [mm/s]
        data = data * 3600*24; % mm/mon
    end
    [a, b, c] = size(data);
    data_subpft_mean = NaN(15,max(repre_site),c);
    data_subpft_std = NaN(15,max(repre_site),c);
    data_mean = NaN(15,c);
    data_std = NaN(15,c);
    for time_index = 1 : c
        for i = 1 : 15
            for j = repre_site_cumsum(i)+1 : repre_site_cumsum(i+1)
                tmp = data(:,:,time_index);
                data_subpft_mean(i,j-repre_site_cumsum(i),time_index) = nanmean(tmp((pft_climate_mask(:,:) == j) & (pft_mask(:,:) == i)));
                data_subpft_std(i,j-repre_site_cumsum(i),time_index) = nanstd(tmp((pft_climate_mask(:,:) == j) & (pft_mask(:,:) == i)));
            end
            data_mean(i,time_index) = nanmean(tmp(pft_mask(:,:) == i));
            data_std(i,time_index) = nanstd(tmp(pft_mask(:,:) == i));
        end
    end

    save('pft_plot_data.mat', 'data_subpft_mean');
    save('pft_plot_data.mat', 'data_subpft_std','-append');
    save('pft_plot_data.mat', 'data_mean','-append');
    save('pft_plot_data.mat', 'data_std','-append');
    movefile('pft_plot_data.mat', [var_name{var} '_pft_plot_data.mat']);

end

% plot
pft_name = {'NET Temp','NET Bor','NDT Bor','BET Tro','BET Temp','BDT Tro','BDT Temp','BDT Bor','BES Temp','BDS Temp','BDS Bor','C3 arc','C3','C4','CropR','CropI'};
%gpp
load gpp_pft_plot_data.mat;
for i = 1 : 3
    figure('pos',[200 200 1800 1200]);
    for j = 1 : 5
        for myindex = (i-1)*5 + j
            subplot(5,1,j);errorbar(data_mean(myindex,:)',data_std(myindex,:)');hold;
            title([pft_name{myindex} '  GPP gC/m2/day'],'FontSize',18,'FontWeight','bold');
        end
    end
end
%nee
load nee_pft_plot_data.mat;
for i = 1 : 3
    figure('pos',[200 200 1800 1200]);
    for j = 1 : 5
        for myindex = (i-1)*5 + j
            subplot(5,1,j);errorbar(data_mean(myindex,:)',data_std(myindex,:)');hold;
            title([pft_name{myindex} '  NEE gC/m2/day'],'FontSize',18,'FontWeight','bold');
        end
    end
end
%et
load et_pft_plot_data.mat;
for i = 1 : 3
    figure('pos',[200 200 1800 1200]);
    for j = 1 : 5
        for myindex = (i-1)*5 + j
            subplot(5,1,j);errorbar(data_mean(myindex,:)',data_std(myindex,:)');hold;
            title([pft_name{myindex} '  ET mm/mon'],'FontSize',18,'FontWeight','bold');
        end
    end
end
%lai
load lai_pft_plot_data.mat;
for i = 1 : 3
    figure('pos',[200 200 1800 1200]);
    for j = 1 : 5
        for myindex = (i-1)*5 + j
            subplot(5,1,j);errorbar(data_mean(myindex,:)',data_std(myindex,:)');hold;
            title([pft_name{myindex} '  LAI m2/m2'],'FontSize',18,'FontWeight','bold');
        end
    end
end
%soilc
load soilc_pft_plot_data.mat;
figure;barwitherr(data_std,data_mean);hold;xlim([0 17]);
ylabel('SoilC KgC/m2','FontSize',18,'FontWeight','bold');
xticklabel_rotate([1:16],45,pft_name,'FontSize',18,'FontWeight','bold');
%for i = 1 : 15
%    figure;barwitherr([0 nanstd(data_subpft_std(i,:))],[data_std(i) nanmean(data_subpft_std(i,:))]);hold;
%    ylabel(pft_name{i},'FontSize',18,'FontWeight','bold');
%end
%biomass
load biomass_pft_plot_data.mat;
figure;barwitherr(data_std,data_mean);hold;xlim([0 17]);
ylabel('Biomass KgC/m2','FontSize',18,'FontWeight','bold');
xticklabel_rotate([1:16],45,pft_name,'FontSize',18,'FontWeight','bold');
%}

end