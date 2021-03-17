function[HDF5]=LoadlittleData2(ThisFilename)
%% Load entire HDF5 file in a struct, and extract the dropshots, remember to update the eventcodes if they change
% Tim Brandt van Driel 2016
HDF5.info=h5info(ThisFilename);
file = H5F.open (ThisFilename, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');


% no groups
if numel(HDF5.info.Groups)==0
    try
        
        for gg=1:numel(HDF5.info.Datasets) % for each group load all datasets
            try
                %n1=HDF5.info.Groups(ii).Name;
                n2=HDF5.info.Datasets(gg).Name;
                %if num2str(n2(1))<10
                %    n2=['_' n2]
                %end
                
                dataset=[num2str(n2)];
                dset = H5D.open (file, dataset);
                dataspace = H5D.get_space(dset);
                tmp= H5D.read(dset, 'H5ML_DEFAULT', 'H5S_ALL', dataspace, 'H5P_DEFAULT');
                H5D.close (dset);
                H5S.close (dataspace);
                
                eval(['HDF5.' num2str(n2) '=tmp;']);
                
            catch
                if strcmp(n2,'91')
                    HDF5.DropL=tmp;
                    %disp([num2str(n1) '/' num2str(n2) ' saved as DropL' ])
                    %elseif strcmp(n2,'96')
                    %    HDF5.DropX=tmp;
                elseif strcmp(n2,'162')
                    HDF5.DropX=tmp;
                    %disp([num2str(n1) '/' num2str(n2) ' saved as DropX' ])
                    
                else
                    %disp([num2str(n1) '/' num2str(n2) ' failed' ])
                end
            end
        end
        % add one deeper after this
        
        
        
        disp('no groups loaded all directly')
    catch
        disp('no groups and failed to load directly')
    end
    
end

% for each HDF5(ii).Datasets
for ii=1:numel(HDF5.info.Groups)
    %tmp=HDF5.info.Groups(ii).Datasets.Name;
    for gg=1:numel(HDF5.info.Groups(ii).Datasets) % for each group load all datasets
       try
        n1=HDF5.info.Groups(ii).Name;
        n2=HDF5.info.Groups(ii).Datasets(gg).Name;
        %if num2str(n2(1))<10
        %    n2=['_' n2]
        %end
        
        dataset=[num2str(n1) '/' num2str(n2)];
        dset = H5D.open (file, dataset);
        dataspace = H5D.get_space(dset);
        tmp= H5D.read(dset, 'H5ML_DEFAULT', 'H5S_ALL', dataspace, 'H5P_DEFAULT');
        H5D.close (dset);
        H5S.close (dataspace);
        
        eval(['HDF5.' num2str(n1(2:end)) '.' num2str(n2) '=tmp;']);
        
       catch
           if strcmp(n2,'91')
               HDF5.DropL=tmp;
               %disp([num2str(n1) '/' num2str(n2) ' saved as DropL' ])
               %elseif strcmp(n2,'96')
               %    HDF5.DropX=tmp;
           elseif strcmp(n2,'162')
               HDF5.DropX=tmp;
               %disp([num2str(n1) '/' num2str(n2) ' saved as DropX' ])
                        
           else
               %disp([num2str(n1) '/' num2str(n2) ' failed' ])
           end
       end
    end
    
    %for each group load all groups - if group go one deeper
    if ~isempty(isempty(HDF5.info.Groups(ii).Groups))
        %gg=1:numel(HDF5.info.Groups(ii).Datasets)
        for gg=1:numel(HDF5.info.Groups(ii).Groups)
            
            
            for hh=1:numel(HDF5.info.Groups(ii).Groups(gg).Datasets)
                try
                    n1=HDF5.info.Groups(ii).Groups(gg).Name;
                    
                    n1str=HDF5.info.Groups(ii).Groups(gg).Name;
                    n1str(n1str=='/')='.';
                    
                    n2=HDF5.info.Groups(ii).Groups(gg).Datasets(hh).Name;
                    
                    %n1=HDF5.info.Groups(ii).Name;
                    %n2=HDF5.info.Groups(ii).Datasets(gg).Name;
                    %if num2str(n2(1))<10
                    %    n2=['_' n2]
                    %end
                    
                    dataset=[num2str(n1) '/' num2str(n2)];
                    dset = H5D.open (file, dataset);
                    dataspace = H5D.get_space(dset);
                    tmp= H5D.read(dset, 'H5ML_DEFAULT', 'H5S_ALL', dataspace, 'H5P_DEFAULT');
                    H5D.close (dset);
                    H5S.close (dataspace);
                    
                    %eval(['HDF5.' num2str(n1(2:end)) '.' num2str(n2) '=tmp;']);
                    eval(['HDF5' num2str(n1str) '.' num2str(n2) '=tmp;']);
                    
                end
            end
        end
   end
    
end
clear tmp
H5F.close(file)