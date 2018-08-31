function step3_generate_surface_domain(gridfile,ensemblefile,nparam)

% input, 1. 'global_sparse_grid.txt': sites information
% input, 2. ALM_param.txt: parameter ensembles
% input, 3, nvar: number of variable want to tune
% output, surface and domain files, ready to be used by ELM-ECA-CNP

ngrids = numel(textread(gridfile,'%1c%*[^\n]')); % determine how many parameter to tune%% user lon lat info

%% read in site infomation
fid = fopen(gridfile,'r');
C = textscan(fid,'%.2f\t%.2f\t%d\t%d\t%d\t%d\t%d\t%d\n',ngrids);
sparse_grid = [double(C{1,1}) double(C{1,2}) double(C{1,5}) double(C{1,6})];
fclose(fid);

lonc = sparse_grid(:,1);
latc = sparse_grid(:,2);

%% read in parameter ensemble information and assign to gridcell level 
nensemble = numel(textread(ensemblefile,'%1c%*[^\n]')); % determine how many parameter to tune%% user lon lat info
nensemble = nensemble - 1; % remove one header line

fid = fopen(ensemblefile,'r');
parname = textscan(fid,'%s',nparam);
parname = char(parname{1,1});

tmp = textscan(fid,'%f');
param_val = reshape(tmp{1,1},nparam,nensemble)';
for i = 1 : nensemble
    param_val_gridcell((i-1)*ngrids+1:i*ngrids,:) = repmat(param_val(i,:),ngrids,1);
end
fclose(fid);

%% create large ensemble, all in nsites
lonc = repmat(lonc,nensemble,1);
latc = repmat(latc,nensemble,1);
nsites = ngrids * nensemble;

%% convert (-180 - 180) to CLM lon (0-360) if necessary
for i = 1 : nsites
    if lonc(i,1) < 0
        lonc(i,1) = 360 + lonc(i,1);
    end
end

%% calculate vertex
latv = zeros(nsites,4);
lonv = zeros(nsites,4);
dlon = 1;
dlat = 1;
for i = 1:nsites
    lonv(i,:) = [lonc(i)-dlon/2 lonc(i)+dlon/2 lonc(i)+dlon/2 lonc(i)-dlon/2];
    latv(i,:) = [latc(i)-dlat/2 latc(i)-dlat/2 latc(i)+dlat/2 latc(i)+dlat/2];
end

%% create domain file
c = clock;
c_year = mod(c(1),100);
c_month = c(2);
c_day = c(3);
fname_out = sprintf('domain_sparse_grids_global_c%.2d%.2d%.2d.nc',c_year,c_month,c_day)
ncid_inp = netcdf.open('domain.lnd.360x720.130305.nc','NC_NOWRITE');
%ncid_inp = netcdf.open('domain.lnd.fv1.9x2.5_gx1v6.090206.nc','NC_NOWRITE');
ncid_out = netcdf.create(fname_out,'NC_CLOBBER');
info_inp = ncinfo('domain.lnd.360x720.130305.nc');
%info_inp = ncinfo('domain.lnd.fv1.9x2.5_gx1v6.090206.nc');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid_inp);
nv = size(lonv,2);
ni = size(lonv,1);
nj = 1;

for ii = 1:ndims
    [dimname, ndim] = netcdf.inqDim(ncid_inp,ii-1);
    switch dimname
        case 'ni'
            ndim = ni;
        case 'nj'
            ndim = nj;
        case 'nv'
            ndim = nv;
    end
    dimid(ii) = netcdf.defDim(ncid_out,dimname,ndim);
end

for ivar = 1:nvars
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncid_inp,ivar-1);
    varid(ivar) = netcdf.defVar(ncid_out,varname,xtype,dimids);
    varnames{ivar} = varname;
    
    for iatt = 1:natts
        attname = netcdf.inqAttName(ncid_inp,ivar-1,iatt-1);
        attvalue = netcdf.getAtt(ncid_inp,ivar-1,attname);
        
        netcdf.putAtt(ncid_out,ivar-1,attname,attvalue);
    end
    
end

varid = netcdf.getConstant('GLOBAL');
[~,user_name]=system('echo $USER');
netcdf.putAtt(ncid_out,varid,'Created_by' ,user_name(1:end-1));
netcdf.putAtt(ncid_out,varid,'Created_on' ,datestr(now,'ddd mmm dd HH:MM:SS yyyy '))

netcdf.endDef(ncid_out);

%R = 6371000;% m
for ivar = 1:nvars
    
    data = netcdf.getVar(ncid_inp,ivar-1);
    [varname,vartype,vardimids,varnatts] = netcdf.inqVar(ncid_inp,ivar-1);
    
    switch varname
        case 'xc'
            data = lonc;
        case 'yc'
            data = latc;
        case 'xv'
            data = lonv;
        case 'yv'
            data = latv;
        case 'mask'
            data = ones(length(lonc),1);
        case 'frac'
            data = ones(length(lonc),1);
        case 'area'
            data = ones(length(lonc),1);
    end
    netcdf.putVar(ncid_out,ivar-1,data);
end
netcdf.close(ncid_inp);
netcdf.close(ncid_out);

%% create surfdata file
fname_out = sprintf('surfdata_sparse_grids_global_c%.2d%.2d%.2d.nc',c_year,c_month,c_day)
ncid_inp = netcdf.open('surfdata_360x720cru_simyr2000_c170913.nc','NC_NOWRITE');
%ncid_inp = netcdf.open('surfdata_1.9x2.5_simyr2000_c170818.nc','NC_NOWRITE');
ncid_out = netcdf.create(fname_out,'NC_CLOBBER');
info_inp = ncinfo('surfdata_360x720cru_simyr2000_c170913.nc');
%info_inp = ncinfo('surfdata_1.9x2.5_simyr2000_c170818.nc');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid_inp);

dimid(1:ndims) = -1;
lonlat_found = 0;

for idim = 1:ndims
    [dimname, dimlen] = netcdf.inqDim(ncid_inp,idim-1);
    %disp(['Inp: Dimension name:' dimname])
    
    switch dimname
        case {'lsmlon','lsmlat'}
            if (strcmp(dimname,'lsmlat'))
                lat_dimid = idim;
            else
                lon_dimid = idim;
            end
            
            if (lonlat_found == 0)
                lonlat_found = 1;
                dimname = 'gridcell';
                dimlen = length(lonc);
                %disp(['Out: Dimension name:' dimname])
                dimid(idim) = netcdf.defDim(ncid_out,dimname,dimlen);
            end
        case 'time'
            %disp(['Out: Dimension name:' dimname])
            dimid(idim) = netcdf.defDim(ncid_out,dimname,netcdf.getConstant('NC_UNLIMITED'));
        otherwise
            %disp(['Out: Dimension name:' dimname])
            for ii=1:length(info_inp.Dimensions)
                if (strcmp(info_inp.Dimensions(ii).Name,dimname) == 1)
                    [dimname, dimlen] = netcdf.inqDim(ncid_inp,ii-1);
                end
            end
            dimid(idim) = netcdf.defDim(ncid_out,dimname,dimlen);
    end
end

%PCT_NAT_var = 0;
for ivar = 1:nvars
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncid_inp,ivar-1);
    %disp(['varname : ' varname ' ' num2str(dimids)])
    if(isempty(dimids)==0)
        if(dimids(1) == 0 && dimids(2) == 1)
            dimids_new =  [0 dimids(3:end)-1];
            dimids = dimids_new;
        else
            dimids = dimids - 1;
        end
    end
    %if ~strcmp(varname,'PCT_NAT_PFT')
    varid(ivar) = netcdf.defVar(ncid_out,varname,xtype,dimids);
    varnames{ivar} = varname;
    %disp([num2str(ivar) ') varname : ' varname ' ' num2str(dimids)])
    
    for iatt = 1:natts
        attname = netcdf.inqAttName(ncid_inp,ivar-1,iatt-1);
        attvalue = netcdf.getAtt(ncid_inp,ivar-1,attname);
        
        netcdf.putAtt(ncid_out,ivar-1,attname,attvalue);
    end
    %else
    %    PCT_NAT_var = 1;
    %end
end
varid = netcdf.getConstant('GLOBAL');

[~,user_name]=system('echo $USER');
netcdf.putAtt(ncid_out,varid,'Created_by' ,user_name(1:end-1));
netcdf.putAtt(ncid_out,varid,'Created_on' ,datestr(now,'ddd mmm dd HH:MM:SS yyyy '));
netcdf.endDef(ncid_out);

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Find the nearest neighbor index for (lonc,latc) within global
% dataset
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% Get Lat/Lon for global 2D grid.
for ivar = 1:length(varnames)
    if(strcmp(varnames{ivar},'LATIXY'))
        latixy = netcdf.getVar(ncid_inp,ivar-1);
    end
    if(strcmp(varnames{ivar},'LONGXY'))
        longxy = netcdf.getVar(ncid_inp,ivar-1);
    end
end

% allocate memoery
ii_idx = zeros(size(lonc));
jj_idx = zeros(size(lonc));

% find the nearest CLM land gridcell index
pftmask = ncread('surfdata_360x720cru_simyr2000_c170913.nc','PFTDATA_MASK');
latixy(pftmask==0) = -9999;
longxy(pftmask==0) = -9999;

for ii=1:size(lonc,1)
    for jj=1:size(lonc,2)
        dist = (longxy - lonc(ii,jj)).^2 + (latixy - latc(ii,jj)).^2;
        [nearest_cell_i_idx, nearest_cell_j_idx] = find( dist == min(min(dist)));
        if (length(nearest_cell_i_idx) > 1)
            disp(['  WARNING: Site with (lat,lon) = (' sprintf('%f',latc(ii,jj)) ...
                sprintf(',%f',lonc(ii,jj)) ') has more than one cells ' ...
                'that are equidistant.' char(10) ...
                '           Picking the first closest grid cell.']);
            for kk = 1:length(nearest_cell_i_idx)
                disp(sprintf('\t\tPossible grid cells: %f %f', ...
                    latixy(nearest_cell_i_idx(kk),nearest_cell_j_idx(kk)), ...
                    longxy(nearest_cell_i_idx(kk),nearest_cell_j_idx(kk))));
            end
        end
        ii_idx(ii,jj) = nearest_cell_i_idx(1);
        jj_idx(ii,jj) = nearest_cell_j_idx(1);
    end
end

mypft = zeros(nsites,17);
sitepft = repmat(sparse_grid(:,3)',1,nensemble);
for i = 1 : nsites
    mypft(i,sitepft(i)+1) = 100;
end

for ivar = 1:nvars
    
    %disp(varnames{ivar})
    [varname,vartype,vardimids,varnatts]=netcdf.inqVar(ncid_inp,ivar-1);
    data = netcdf.getVar(ncid_inp,ivar-1);
    switch varname
        case {'LATIXY'}
            netcdf.putVar(ncid_out,ivar-1,latc);
        case {'LONGXY'}
            netcdf.putVar(ncid_out,ivar-1,lonc);
        case {'PCT_NAT_PFT'}
            netcdf.putVar(ncid_out,ivar-1,mypft);
        otherwise
            
            switch varname
                case {'PCT_URBAN','PCT_CROP','PCT_WETLAND','PCT_LAKE','PCT_GLACIER',...
                        'PCT_GLC_MEC'}
                    data = data*0;
                case {'PCT_NATVEG'}
                    data = data*0 + 100;
            end
            
            switch varname
                case {'PFTDATA_MASK', 'LANDFRAC_PFT'}
                    data = data*0 + 1;
            end
            
            switch length(vardimids)
                case 0
                    netcdf.putVar(ncid_out,ivar-1,data);
                case 1
                    data = 0;
                    netcdf.putVar(ncid_out,ivar-1,0,length(data),data);
                case 2
                    if (min(vardimids) == 0)
                        data_2d = zeros(size(lonc));
                        for ii=1:size(lonc,1)
                            for jj=1:size(lonc,2)
                                data_2d(ii,jj) = data(ii_idx(ii,jj),jj_idx(ii,jj));
                            end
                        end
                        
                        % (lon,lat) --> % (gridcell)
                        vardimids_new =  [0 vardimids(3:end)-1];
                        vardimids = vardimids_new;
                        dims = size(data_2d);
                        if (length(dims)>2)
                            dims_new = [dims(1)*dims(2) dims(3:end)];
                        else
                            dims_new = [dims(1)*dims(2) 1];
                        end
                        data_2d_new = reshape(data_2d,dims_new);
                        data_2d = data_2d_new;
                        
                        netcdf.putVar(ncid_out,ivar-1,data_2d);
                    else
                        netcdf.putVar(ncid_out,ivar-1,data);
                    end
                case 3
                    if (min(vardimids) == 0)
                        nx = size(lonc,1);
                        ny = size(lonc,2);
                        nz = size(data,3);
                        data_3d = zeros(nx,ny,nz);
                        for ii = 1:nx
                            for jj = 1:ny
                                for kk = 1:nz
                                    data_3d(ii,jj,kk) = data(ii_idx(ii,jj),jj_idx(ii,jj),kk);
                                end
                            end
                        end
                        
                        % (lon,lat,:) --> % (gridcell,:)
                        vardimids_new =  [0 vardimids(3:end)-1];
                        vardimids = vardimids_new;
                        dims = size(data_3d);
                        if (length(dims)>2)
                            dims_new = [dims(1)*dims(2) dims(3:end)];
                        else
                            dims_new = [dims(1)*dims(2) 1];
                        end
                        data_3d_new = reshape(data_3d,dims_new);
                        data_3d = data_3d_new;
                        
                        netcdf.putVar(ncid_out,ivar-1,data_3d);
                    else
                        netcdf.putVar(ncid_out,ivar-1,data);
                    end
                case 4
                    if (min(vardimids) == 0)
                        nx = size(lonc,1);
                        ny = size(lonc,2);
                        nz = size(data,3);
                        na = size(data,4);
                        data_4d = zeros(nx,ny,nz,na);
                        for ii = 1:nx
                            for jj = 1:ny
                                for kk = 1:nz
                                    for ll = 1:na
                                        data_4d(ii,jj,kk,ll) = data(ii_idx(ii,jj),jj_idx(ii,jj),kk,ll);
                                    end
                                end
                            end
                        end
                        
                        % (lon,lat,:) --> % (gridcell,:)
                        vardimids_new =  [0 vardimids(3:end)-1];
                        vardimids = vardimids_new;
                        dims = size(data_4d);
                        if (length(dims)>2)
                            dims_new = [dims(1)*dims(2) dims(3:end)];
                        else
                            dims_new = [dims(1)*dims(2) 1];
                        end
                        data_4d_new = reshape(data_4d,dims_new);
                        data_4d = data_4d_new;
                        
                        netcdf.putVar(ncid_out,ivar-1,zeros(length(size(data_4d)),1)',size(data_4d),data_4d);
                    else
                        netcdf.putVar(ncid_out,ivar-1,data);
                    end
                otherwise
                    disp('error')
            end
    end
end
netcdf.close(ncid_inp);
netcdf.close(ncid_out);

%% add new varialbes: parameters value for each sites
for i = 1 : nparam
    par_name = strtrim(parname(i,:));
    nccreate(fname_out,par_name,'Dimensions', {'gridcell',nsites});
    ncwrite(fname_out,par_name,param_val_gridcell(:,i)');
end

end