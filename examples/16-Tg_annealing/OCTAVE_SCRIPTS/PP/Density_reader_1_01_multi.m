function [ groupeddata,temps ] = Density_reader_1_01_multi(foldername,ftype)
%Loads density/temperature data from various file formats
%output is a data array with a format that is readable by Tg analysis
%scripts


if strcmp(ftype,'LAMMPScustom') == 1
    direc='';
    prefix=foldername;
    fpathoffset=2;
    dfile='densities.txt';
    
    
    %Search routine for all data folders
    %basically, the goal of the next few lines is
    %to find all the folders associated with the
    %raw data.
    if (isunix+ismac) == 1         %i.e. if this is a unix
        dum=what;
        olddirec=dum.path;
        slashsym='/';
        direc=strcat(prefix,direc);
        %cd(strcat(olddirec,'/',direc));
        cd(direc)
        dum=dir;
        numfolders=max(size(dum))-2;
        for jjjj=1:numfolders
%            fpaths{jjjj,:}= dum(jjjj+fpathoffset).name;
            fpaths{jjjj}= dum(jjjj+fpathoffset).name;
        end
        cd(olddirec);
        indexcorrector=0;
    end
    if ispc == 1          
        indexcorrector=2;
        slashsym='\';
        direc=strcat(prefix,direc,'\');
        fpaths=ls(direc);           %this command looks for all the folders in the directory, saves them in an array
        dummy=size(fpaths);         %finds the number of folders + 2
        numfolders=dummy(1)-2;      %number of folders in the directory containing all the data
    end
    
    
    for kkll=1:numfolders
        jjjj=kkll+indexcorrector;
        
        %Now we load the density vs temperature data
        if ispc == 1
            filename=strcat(direc,fpaths(jjjj,:),'\',dfile)
        end
        if (ismac+isunix) == 1
%            filename=char(strcat(direc,fpaths(jjjj,:),'/',dfile))
            filename=char(strcat(direc,fpaths{jjjj},'/',dfile))            
        end
        data=Density_reader_1_01C(filename);     %opens the density versus temp data
        numtemps=numel(data(:,1));
        data=flipud(data); %reorganizing the data a bit to make it easier to handle
        
        if kkll == 1
            groupeddata=zeros(numtemps,numfolders);
            groupeddata(:,1)=data(:,1);
            temps=data(:,2);
        end

        if sum(temps == data(:,2)) ~= numel(temps)
            disp('Error: datasets have different numbers of temperatures')
            disp(filename)
            die
        end
        
        if kkll> 1
            groupeddata(:,kkll) = data(:,1);
        end
        
        
    end
   
end





end

