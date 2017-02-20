clear 
clc
close all

%% ====================================================================== %%
%% Create Grid
showgridflag                = 0; % set to 1 to see grids
nreps                       = 100; 
includeendcasesflag         = 0; % set to 1 to replace 1st and 2nd K grids with parallel and series

nColumns                    = 25;
nRows                       = nColumns;
nGrids                      = nColumns*nColumns;
nPoints                     = [1:nGrids]';
Kgrid                       = reshape(nPoints,[nRows nColumns]);

kHigh                       = 1;
kLow                        = 0.01;  
fractionhighK               = [0.25, 0.5, 0.75]; % the example has to be 25%,50% and 75%
fractionhighK               = fractionhighK';
numberhighK                 = floor(nGrids*fractionhighK);
numberlowK                  = nGrids-numberhighK;
fractionhighK               = fractionhighK';


figure('Color',[1 1 1]);
for j=1:size(fractionhighK,2);
    for iii=1:nreps
        
        Keff_parallel(j)=fractionhighK(:,j).*kHigh+(1-fractionhighK(:,j)).*kLow;%parallel system alpha =1
        Keff_series(j)=1/((fractionhighK(:,j)./kHigh)+((1-fractionhighK(:,j))./kLow));%series in layer alpha =-1
        Keff_geome(j)=((fractionhighK(:,j).*kHigh^0.5)+((1-fractionhighK(:,j)).*kLow^0.5)).^(1/0.5);% geometric random distr alpha = 0.5
       
%         ReSamples                   = (ceil(rand(nGrids,nResampling).*nGrids));
%         ReSamples                   = reshape(ReSamples(:,1),[nRows nColumns]);
        
        %% Assign K values to the matrix
        %     Kgrid(1:25,1:25)=1;
        %     Kgrid(find(ReSamples<=nGrids-numberhighK))=0.01; %or 0.5 in the case
        
        % put the following into the loop to replace your matrix definition
        tempKgrid=rand(nColumns,nRows); % make a grid of random numbers between 0 and 1
        tempKgridunwrap=reshape(tempKgrid,1,nColumns*nRows); % unwrap grid into a vector
        temp1=sort(tempKgridunwrap); % sort the vector
        tempKgrid2(1:nColumns, 1:nRows)=kHigh; % first fill new K grid with all high K
        tempKgrid2(find(tempKgrid<0.0000001+temp1(numberlowK(j))))=kLow; % substitute lowest random values with low K
        Kgrid=tempKgrid2;
        
        if includeendcasesflag==1; % set flag to one to have case 1 parallel and 2 series
            if iii==1; % make a parallel example
                tempKgrid2(1:nColumns, 1:nRows)=kLow; % first fill new K grid with all high K
                nkHighrows=floor(numberhighK/nColumns);
                nkHighextracolumns=numberhighK-nkHighrows*nColumns;
                tempKgrid2(1:nkHighrows,:)=kHigh;
                tempKgrid2(nkHighrows+1,1:nkHighextracolumns)=kHigh;
                Kgrid=tempKgrid2;
            end;
            
            if iii==2; % make a series example
                tempKgrid2(1:nColumns, 1:nRows)=kLow; % first fill new K grid with all high K
                nkHighcolumns=floor(numberhighK/nRows);
                nkHighextrarows=numberhighK-nkHighcolumns*nRows;
                tempKgrid2(:,1:nkHighcolumns)=kHigh;
                tempKgrid2(1:nkHighextrarows,nkHighcolumns+1)=kHigh;
                Kgrid=tempKgrid2;
            end;
        end;
        
        %% Plot K Grid
        if showgridflag==1
            figure('Color',[1 1 1]);
            imagesc((Kgrid));  %instead of pcolor
            view([0 -90])
            title('K-Grid','FontWeight','bold','FontSize',14,...
                'FontName','Times New Roman');
            set(gca,'FontName','Times New Roman');
        end;
        
        %% Create Modflow Input FileholdKgrid
        fileID = fopen('PJ1.bcf6','w');
        Header1 = [-1 1.00E+30 0 0 0 0];
        Header2 = 10;
        Header3 = 'CONSTANT 1.0 TRPY';
        Header4 = 'INTERNAL 1 (FREE) 1 LAYER 1';
        dlmwrite('PJ1.bcf6',Header1,'delimiter','\t','precision',4)
        dlmwrite('PJ1.bcf6',Header2,'-append','delimiter','\t','precision',4)
        dlmwrite('PJ1.bcf6',Header3,'-append','delimiter','','precision',4)
        dlmwrite('PJ1.bcf6',Header4,'-append','delimiter','','precision',4)
        dlmwrite('PJ1.bcf6',Kgrid,'-append','delimiter','\t','precision',4)
        fclose(fileID);
        
        %% Run MODFLOW from MATLAB
        %     addpath(genpath('C:\Users\TUCSON\Desktop\thesis\MODFLOW EXAMPLE'))
        %     cd 'C:\Users\TUCSON\Desktop\thesis\MODFLOW EXAMPLE'
        %     !mf.exe
        %     system('mf.exe PJ1 &');
        !mf.exe PJ1
        %!~/modflow/bin/mf2005
        disp('MODFLOW Runs completed...')
        
        
        %% Import MODFLOW Results
        
        % Import flow
        Flow = Import_Modflow_Result('PJ1.lst');
        Flow
        
        % Write a subroutine to read in the Head value
        Heads = importHead('PJ1.lst', 322, 396);
        HeadsResh  = reshape(Heads',36,25);
        Head = (HeadsResh([2:23 25:27],:))';
        
        
        %Calculate K effective
        Keff = (nColumns-1)*Flow/(nRows*100);
        holdKeff(j,iii)=Keff;
        holdKgrid(j,iii,1:nColumns,1:nRows)=Kgrid;
        holdFlow(j,iii)=Flow;
        holdHead(j,iii,:,:)=Head;
        
    end

        %% Plot Histogram
        subplot(3,1,j) %3
        % hist(holdKeff);
        h = histfit(holdKeff(j,:),15,'normal');
        yL = get(gca,'YLim');
        line([Keff_series(j) Keff_series(j)],yL,'Color','b');
        line([Keff_geome(j) Keff_geome(j)],yL,'Color','g');
        line([Keff_parallel(j) Keff_parallel(j)],yL,'Color','y');
        set(gca, 'Xlim',[Keff_series(j) Keff_parallel(j)]);
        xlabel('Keff');
        ylabel('Frequency');
        axis([0 1 0 20])
        title('Histogram & Normal Density','FontWeight','bold','FontSize',14,...
        'FontName','Times New Roman');
        
   %%  K effective 
    
    Keff_low(j,:)  = nonzeros(sort(holdKeff(j,:)));
    Keff_high(j,:) = nonzeros(sort(holdKeff(j,:),'descend'));
      
end
four_Klow = Keff_low(:,1:4);
four_Khigh = Keff_high(:,1:4);


%% Figures Fourth K low (25%,50%,75% fraction K high)
for j=1:size(fractionhighK,2)
    
    for k = 1:4;  
        
       figure(1+j)

       subplot(2,2,k)                                       % find grid with lowest Keff
       LowestKrep(j,k)=find(holdKeff(j,:)==four_Klow(j,k)); 
       hh(j) = imagesc((squeeze(holdKgrid(j,LowestKrep(j,k),:,:))));
       view([0 90])
       title('Grid for Low Keff','FontWeight','bold','FontSize',12,...
       'FontName','Times New Roman');
       set(gca,'FontName','Times New Roman');
       myColorMap = jet(2);
       colormap(myColorMap);
       L = line(nan(2), nan(2),'LineStyle','none'); % 'nan' creates 'invisible' data
       set(L, {'MarkerEdgeColor'}, num2cell(myColorMap, 2),...
       {'MarkerFaceColor'},num2cell(myColorMap, 2),... % setting the markers to filled squares
       'Marker','s'); 
       text(17,-1,[num2str(four_Klow(j,k))]);
    end
       legend(L, {'kLow','kHigh'})
end
       

%% Figures Fourth K high (25%,50%,75% fraction K high)
for j=1:size(fractionhighK,2)
    
    for k = 1:4;
       figure(4+j) 

       subplot(2,2,k)                                       % find grid with highest Keff
       HighestKrep(j,k)=find(holdKeff(j,:)==four_Khigh(j,k)); 
       hh(j) = imagesc((squeeze(holdKgrid(j,HighestKrep(j,k),:,:))));
       view([0 90])
       title('Grid for High Keff','FontWeight','bold','FontSize',12,...
       'FontName','Times New Roman');
       set(gca,'FontName','Times New Roman');
       myColorMap = jet(2);
       colormap(myColorMap);
       L = line(nan(2), nan(2),'LineStyle','none'); % 'nan' creates 'invisible' data
       set(L, {'MarkerEdgeColor'}, num2cell(myColorMap, 2),...
       {'MarkerFaceColor'},num2cell(myColorMap, 2),... % setting the markers to filled squares
       'Marker','s'); 
       text(17,-1,[num2str(four_Khigh(j,k))]);  
    end
       legend(L, {'kLow','kHigh'})
end


%% Figures 8 panel plots for Grid of Khigh and kLow distribution

figure('Color',[1 1 1]);
for k = 1:4;
        for j=1;
            
        subplot(6,4,k)                                  % find grid with lowest Keff
        
        LowestKrep1=find(holdKeff(j,:)==four_Klow(j,k)); 
        imagesc((squeeze(holdKgrid(j,LowestKrep1,:,:)))); %lowestKrep1
        view([0 90])
        title('Grid for Low Keff','FontWeight','bold','FontSize',10,...
        'FontName','Times New Roman');
        set(gca,'FontName','Times New Roman');
        myColorMap = jet(2);
        colormap(myColorMap);
        L = line(nan(2), nan(2),'LineStyle','none'); % 'nan' creates 'invisible' data
        set(L, {'MarkerEdgeColor'}, num2cell(myColorMap, 2),...
        {'MarkerFaceColor'},num2cell(myColorMap, 2),... % setting the markers to filled squares
        'Marker','s');
        text(20,-3,[num2str(four_Klow(j,k))]); % consider low k effective value
%       legend(L, {'kLow','kHigh'})
        
        subplot(6,4,k+4)                                % find grid with highest Keff
        
        HighestKrep1=find(holdKeff(j,:)==four_Khigh(j,k)); 
        imagesc((squeeze(holdKgrid(j,HighestKrep1,:,:))));
        view([0 90])
        title('Grid for High Keff','FontWeight','bold','FontSize',10,...
        'FontName','Times New Roman');
        set(gca,'FontName','Times New Roman');
        myColorMap = jet(2);
        colormap(myColorMap);
        L = line(nan(2), nan(2),'LineStyle','none'); % 'nan' creates 'invisible' data
        set(L, {'MarkerEdgeColor'}, num2cell(myColorMap, 2),...
        {'MarkerFaceColor'},num2cell(myColorMap, 2),... % setting the markers to filled squares
        'Marker','s');
        text(20,-3,[num2str(four_Khigh(j,k))]); % consider high k effective value
%       legend(L, {'kLow','kHigh'})
    end
        for j=2;
        
        subplot(6,4,k+8)                               % find grid with lowest Keff
        
        LowestKrep1=find(holdKeff(j,:)==four_Klow(j,k)); 
        imagesc((squeeze(holdKgrid(j,LowestKrep1,:,:))));
        view([0 90])
        title('Grid for Low Keff','FontWeight','bold','FontSize',10,...
        'FontName','Times New Roman');
        set(gca,'FontName','Times New Roman');
        myColorMap = jet(2);
        colormap(myColorMap);
        L = line(nan(2), nan(2),'LineStyle','none'); % 'nan' creates 'invisible' data
        set(L, {'MarkerEdgeColor'}, num2cell(myColorMap, 2),...
        {'MarkerFaceColor'},num2cell(myColorMap, 2),... % setting the markers to filled squares
        'Marker','s');
        text(20,-3,[num2str(four_Klow(j,k))]); % consider low k effective value
%       legend(L, {'kLow','kHigh'})
        
        subplot(6,4,k+12)                                 % find grid with highest Keff
        
        HighestKrep1=find(holdKeff(j,:)==four_Khigh(j,k)); 
        imagesc((squeeze(holdKgrid(j,HighestKrep1,:,:))));
        view([0 90])
        title('Grid for High Keff','FontWeight','bold','FontSize',10,...
        'FontName','Times New Roman');
        set(gca,'FontName','Times New Roman');
        myColorMap = jet(2);
        colormap(myColorMap);
        L = line(nan(2), nan(2),'LineStyle','none'); % 'nan' creates 'invisible' data
        set(L, {'MarkerEdgeColor'}, num2cell(myColorMap, 2),...
        {'MarkerFaceColor'},num2cell(myColorMap, 2),... % setting the markers to filled squares
        'Marker','s');
        text(20,-3,[num2str(four_Khigh(j,k))]); % consider high k effective value
%       legend(L, {'kLow','kHigh'})
    end
        for j=3;

        subplot(6,4,k+16)                               % find grid with lowest Keff
        
        LowestKrep1=find(holdKeff(j,:)==four_Klow(j,k)); 
        imagesc((squeeze(holdKgrid(j,LowestKrep1,:,:))));
        view([0 90])
        title('Grid for Low Keff','FontWeight','bold','FontSize',10,...
        'FontName','Times New Roman');
        set(gca,'FontName','Times New Roman');
        myColorMap = jet(2);
        colormap(myColorMap);
        L = line(nan(2), nan(2),'LineStyle','none'); % 'nan' creates 'invisible' data
        set(L, {'MarkerEdgeColor'}, num2cell(myColorMap, 2),...
        {'MarkerFaceColor'},num2cell(myColorMap, 2),... % setting the markers to filled squares
        'Marker','s');
        text(20,-3,[num2str(four_Klow(j,k))]); % consider low k effective value
%       legend(L, {'kLow','kHigh'})
            
        subplot(6,4,k+20)                                 % find grid with highest Keff

        HighestKrep1=find(holdKeff(j,:)==four_Khigh(j,k)); 
        imagesc((squeeze(holdKgrid(j,HighestKrep1,:,:))));
        view([0 90])
        title('Grid for High Keff','FontWeight','bold','FontSize',10,...
        'FontName','Times New Roman');
        set(gca,'FontName','Times New Roman');
        myColorMap = jet(2);
        colormap(myColorMap);
        L = line(nan(2), nan(2),'LineStyle','none'); % 'nan' creates 'invisible' data
        set(L, {'MarkerEdgeColor'}, num2cell(myColorMap, 2),...
        {'MarkerFaceColor'},num2cell(myColorMap, 2),... % setting the markers to filled squares
        'Marker','s');
        text(20,-3,[num2str(four_Khigh(j,k))]); % consider low k effective value
        end
         
end
        legend(L, {'kLow','kHigh'})


%% Flux computed based on Head of Kgrid for The Lowest and Highest K effective

for j=1:size(fractionhighK,2);
    for k=1:nreps;
                    [FX(j,k,:,:),FY(j,k,:,:)]   = gradient(squeeze(holdHead(j,k,:,:)));
                    qx(j,k,:,:)                 = (FX(j,k,:,:)).*((holdKgrid(j,k,:,:)));
                    qy(j,k,:,:)                 = (FY(j,k,:,:)).*((holdKgrid(j,k,:,:)));
                    qtotal(j,k,:,:)             = sqrt(qx(j,k,:,:).^2+qy(j,k,:,:).^2);  
    end
end

%% Here we extract the 4 max and 4 min flux (q)

for j=1:size(fractionhighK,2);
    for k=1:nreps;
        for  m = 1:4;
                    qtotal_Min(j,m,:,:)       = qtotal(j,LowestKrep(j,m),:,:);                  
                    qtotal_Max(j,m,:,:)       = qtotal(j,HighestKrep(j,m),:,:);
        end
    end
end

%% look at which level of velocity we should use for velocity binning

figure(200)
c = 1;
for jj=1:size(qtotal,1);
    for ii=1:size(qtotal,2);
        
        chooselevel=jj;     %chooselevel = 3 conditions 25%HighK, 50%HighK, and 75%HighK
        chooserealiz=ii;    %chooserealiz = number of repetitions 
        qtotalrealiz=squeeze(qtotal(chooselevel,chooserealiz,:,:));
        unwrapqtotal=reshape(qtotalrealiz,1,size(qtotalrealiz,1)*size(qtotalrealiz,2));

%       figure(100)
%       hist(unwrapqtotal);

        unwrapqtotalfinal(c,:)=sort(squeeze(unwrapqtotal));

        plot(unwrapqtotalfinal(c,:)); %'b.'
        hold on; 
        c = c+1;

    end
end
hold off

dataqtotal = squeeze(unwrapqtotalfinal');

qtotal1 = dataqtotal(:,[1:nreps]);         % qtotal data for first n reps (25% highK)
qtotal2 = dataqtotal(:,nreps+1:nreps*2);   % qtotal data for the second n reps (50% highK)
qtotal3 = dataqtotal(:,nreps*2+1:nreps*3); % qtotal data for the third n reps (75% highK)

% Plot figure of level represents low, intermediate, and high velocity
% overqtotal data. Use the plot that we developed together to decide how to define three classes of velocity ... 
%very high, very low, and intermediate. Then use red, yellow, green in your velocity field to
%show the slow, medium, and fast regions for the low and high effective K values.

%% plot low, intermediate and high velocity for 25% highK
% j = 1
figure(201) 
plot(1:492,qtotal1(1:492,:),'r.');  
hold on
plot(493:600,qtotal1(493:600,:),'y.')
plot(601:625,qtotal1(601:end,:),'g.');
%% j = 2
% plot low, intermediate and high velocity for 50% highK
plot(1:312,qtotal2(1:312,:),'r.');
plot(313:600,qtotal2(313:600,:),'y.')
plot(601:625,qtotal2(601:end,:),'g.'); 

%% j = 3
% plot low, intermediate and high velocity for 75% highK
plot(1:155,qtotal3(1:155,:),'r.');
plot(156:600,qtotal3(156:600,:),'y.')
plot(601:625,qtotal3(601:end,:),'g.');
hold off;
legend;
set(legend, 'Ytick'); %// 3 yticks, each "in the middle" of one color
set(legend, 'YTickLabel', {'low','medium','high'});
box on


%% qtotal_Min to look at which velocity we should consider the limit of low,
% medium, and high velocity for cases with 

figure(300)
c = 1;
for jj=1:size(qtotal_Min,1);
    for ii=1:size(qtotal_Min,2);
        
        chooselevel=jj;     %chooselevel = 3 conditions 25%HighK, 50%HighK, and 75%HighK
        chooserealiz=ii;    %chooserealiz = number of repetitions (10)
        qtotal_Minrealiz=squeeze(qtotal_Min(chooselevel,chooserealiz,:,:));
        unwrapqtotal_Min=reshape(qtotal_Minrealiz,1,size(qtotal_Minrealiz,1)*size(qtotal_Minrealiz,2));

%       figure(100)
%       hist(unwrapqtotal);

        unwrapqtotal_Minfinal(c,:)=sort(squeeze(unwrapqtotal_Min));

        plot(unwrapqtotal_Minfinal(c,:)); %'k.'
        hold on; 
        c = c+1;

    end
end
hold off

dataqtotal_Min = squeeze(unwrapqtotal_Minfinal');

qtotal_Min1 = dataqtotal_Min(:,[1:4]);     % qtotal data for the first 4 qtotal_Min (25% highK)
qtotal_Min2 = dataqtotal_Min(:,4+1:4*2);   % qtotal data for the second 4 qtotal_Min (50% highK)
qtotal_Min3 = dataqtotal_Min(:,4*2+1:4*3); % qtotal data for the third 4 qtotal_Min (75% highK)

% Stablish Limit for Low and High Fluxes

%Low Fluxes

%limit fluxes for qtotal_Min
qtotal_Min1_lim1 = qtotal_Min1(492,:); % first limit qtotal Min 25% high K
% qtotal_Min1_lim1 = qtotal_Min1_lim1';  % trasponse of first limit
qtotal_Min1_lim2 = qtotal_Min1(600,:); % second limit qtotal Min 25% high K
% qtotal_Min1_lim2 = qtotal_Min1_lim2';  % trasponse of second limit

qtotal_Min2_lim1 = qtotal_Min2(315,:); % first limit qtotal Min 50% high K
% qtotal_Min2_lim1 = qtotal_Min2_lim1';  % trasponse of first limit
qtotal_Min2_lim2 = qtotal_Min2(600,:); % second limit qtotal Min 50% high K
% qtotal_Min2_lim2 = qtotal_Min2_lim2';  % trasponse of second limit

qtotal_Min3_lim1 = qtotal_Min3(157,:); % first limit qtotal Min 75% high K 
% qtotal_Min3_lim1 = qtotal_Min3_lim1';  % trasponse of first limit
qtotal_Min3_lim2 = qtotal_Min3(600,:); % second limit qtotal Min 75% high K 
% qtotal_Min3_lim2 = qtotal_Min3_lim2';  % trasponse of second limit

% qtotal_Max to look at which velocity we should consider the limit of low,
% medium, and high velocity for cases with 


%Plot fluxes 4 qtotal_Min

for j=1:size(fractionhighK,2)
    for k = 1:4;
        figure(7+j) %8
        subplot(2,2,k)
        if j==1;
            plotFluxes(squeeze(qtotal_Min(j,k,:,:)),1,4);
            view([0 90])
            text(18,-0.7,[num2str(four_Klow(j,k))]); % consider low k effective value         
        end
        if j==2;
            plotFluxes(squeeze(qtotal_Min(j,k,:,:)),1,4);
            view([0 90])
            text(18,-0.7,[num2str(four_Klow(j,k))]); % consider low k effective value         
        end    
        if j==3;
            plotFluxes(squeeze(qtotal_Min(j,k,:,:)),1,4);
            view([0 90])
            text(18,-0.7,[num2str(four_Klow(j,k))]); % consider low k effective value         
        end
        
    end
end

%   legend(hcb, {'Low','kMedium','kHigh'})


%% qtotal_Max to look at which velocity we should consider the limit of low,
% medium, and high velocity for cases with 

figure(400)
c = 1;
for jj=1:size(qtotal_Max,1);
    for ii=1:size(qtotal_Max,2);
        
        chooselevel=jj;     %chooselevel = 3 conditions 25%HighK, 50%HighK, and 75%HighK
        chooserealiz=ii;    %chooserealiz = number of repetitions (10)
        qtotal_Maxrealiz=squeeze(qtotal_Max(chooselevel,chooserealiz,:,:));
        unwrapqtotal_Max=reshape(qtotal_Maxrealiz,1,size(qtotal_Maxrealiz,1)*size(qtotal_Maxrealiz,2));

%       figure(100)
%       hist(unwrapqtotal);

        unwrapqtotal_Maxfinal(c,:)=sort(squeeze(unwrapqtotal_Max));

        plot(unwrapqtotal_Maxfinal(c,:)); %'k.'
        hold on; 
        c = c+1;

    end
end
hold off

dataqtotal_Max = squeeze(unwrapqtotal_Maxfinal');

qtotal_Max1 = dataqtotal_Max(:,[1:4]);     % qtotal data for first 4 qtotal_Max (25% highK)
qtotal_Max2 = dataqtotal_Max(:,4+1:4*2);   % qtotal data for the second 4 qtotal_Max (50% highK)
qtotal_Max3 = dataqtotal_Max(:,4*2+1:4*3); % qtotal data for the third 4 qtotal_Max (75% highK)

%% High Fluxes

%limit fluxes for qtotal_Min
qtotal_Max1_lim1 = qtotal_Max1(492,:); % first limit qtotal Max 25% high K
% qtotal_Max1_lim1 = qtotal_Max1_lim1';  % trasponse of first limit
qtotal_Max1_lim2 = qtotal_Max1(600,:); % second limit qtotal Max 25% high K
% qtotal_Max1_lim2 = qtotal_Max1_lim2';  % trasponse of second limit

qtotal_Max2_lim1 = qtotal_Max2(315,:); % first limit qtotal Max 50% high K
% qtotal_Max2_lim1 = qtotal_Max2_lim1';  % trasponse of first limit
qtotal_Max2_lim2 = qtotal_Max2(600,:); % second limit qtotal Max 50% high K
% qtotal_Max2_lim2 = qtotal_Max2_lim2';  % trasponse of second limit

qtotal_Max3_lim1 = qtotal_Max2(157,:); % first limit qtotal Max 75% high K
% qtotal_Max3_lim1 = qtotal_Max3_lim1';  % trasponse of second limit
qtotal_Max3_lim2 = qtotal_Max2(600,:); % second limit qtotal Max 75% high K 
% qtotal_Max3_lim2 = qtotal_Max3_lim2';  % trasponse of second limit

%plot fluxes 4 qtotal_Max

for j=1:size(fractionhighK,2)
    for k = 1:4;
        figure(10+j) %8
        subplot(2,2,k)
        if j==1;          % 25% high K
            plotFluxes(squeeze(qtotal_Max(j,k,:,:)),1,4);
            view([0 90])
            text(18,-0.7,[num2str(four_Khigh(j,k))]); % consider low k effective value         
        end
        if j==2;          % 50% high K
            plotFluxes(squeeze(qtotal_Max(j,k,:,:)),1,4);
            view([0 90])
            text(18,-0.7,[num2str(four_Khigh(j,k))]); % consider low k effective value         
        end    
        if j==3;          % 75% high K
            plotFluxes(squeeze(qtotal_Max(j,k,:,:)),1,4);
            view([0 90])
            text(18,-0.7,[num2str(four_Khigh(j,k))]); % consider low k effective value         
        end
        
    end
end

%% Plot comparison between distribution of qtotal, qtotal_Min and qtotal_Max

figure(500)
plot(dataqtotal,'r');
hold on
plot(dataqtotal_Min,'k');
plot(dataqtotal_Max,'b');
hold off

        
%% 8 panels for qflux

figure('Color',[1 1 1]);
for k = 1:4;
        for j=1;
            
        subplot(6,4,k)                             % find grid with qflux with lowest Keff
        plotFluxes(squeeze(qtotal_Min(j,k,:,:)),1,4)
        view([0 90])
        text(25,-5,[num2str(four_Klow(j,k))])
        title('Flux for Lowest K','FontWeight','bold','FontSize',12,...
       'FontName','Times New Roman');
   
        subplot(6,4,k+4)                            % find grid with qflux with highest Keff
        plotFluxes(squeeze(qtotal_Max(j,k,:,:)),1,4)
        view([0 90])
        text(25,-5,[num2str(four_Khigh(j,k))])
        title('Flux for Highest K','FontWeight','bold','FontSize',12,...
       'FontName','Times New Roman');
   
    end
        for j=2;
        
        subplot(6,4,k+8)                              % find grid with qflux with lowest Keff
        plotFluxes(squeeze(qtotal_Min(j,k,:,:)),1,4)
        view([0 90])
        text(25,-5,[num2str(four_Klow(j,k))])
        title('Flux for Lowest K','FontWeight','bold','FontSize',12,...
       'FontName','Times New Roman');
   
        
        subplot(6,4,k+12)                             % find grid with qflux with highest Keff
        plotFluxes(squeeze(qtotal_Max(j,k,:,:)),1,4)
        view([0 90])
        text(25,-5,[num2str(four_Khigh(j,k))])
        title('Flux for Highest K','FontWeight','bold','FontSize',12,...
       'FontName','Times New Roman');

    end
        for j=3;

        subplot(6,4,k+16)                         % find grid with qflux with lowest Keff
        plotFluxes(squeeze(qtotal_Min(j,k,:,:)),1,4)
        view([0 90])
        text(25,-5,[num2str(four_Klow(j,k))])
        title('Flux for Lowest K','FontWeight','bold','FontSize',12,...
       'FontName','Times New Roman');
            
        subplot(6,4,k+20)                         % find grid with qflux with highest Keff
        plotFluxes(squeeze(qtotal_Max(j,k,:,:)),1,4)
        view([0 90])
        text(25,-5,[num2str(four_Khigh(j,k))])
        title('Flux for Highest K','FontWeight','bold','FontSize',12,...
       'FontName','Times New Roman');
   
    end
        
end

% hcb = colorbar;
% set(hcb, 'Ytick'); %// 4 yticks, each "in the middle" of one color
% set(hcb, 'YTickLabel', {'low','medium','high'});
% %Red, yellow, green in your velocity field to show the slow, medium, and fast regions


%% check the number of high K in each grid

firstsum=squeeze(sum(holdKgrid,4));
sumKvals=squeeze(sum(firstsum,3));
aveKvals=sumKvals/nRows/nColumns;
normaveKvals=(aveKvals-kLow)/(kHigh-kLow);

%% Generate a color distribution of low, intermedia, and high velocity vector

tempvar=reshape(qtotalrealiz,1,size(qtotalrealiz,1)*size(qtotalrealiz,2)); % turn qs into a vector
% hist(tempvar,4);
% [A B]=hist(tempvar,3)

% use kmeans clustering to define three groups
opts = statset('Display','final');
[cidx, ctrs] = kmeans(tempvar', 3, 'Distance','city','Replicates',5,...
'Options',opts);

tempvar(2,:)=cidx;
tempvar=sort(tempvar',1);
xvals=1:size(tempvar,1);

figure(550)
colorstr='ryg';

for ii=1:size(tempvar,1)
    plotstr=strcat('.',+colorstr(tempvar(ii,2)));
    plot(xvals(ii),tempvar(ii,1),plotstr);
    hold on;
end;

%find limit of the tempvar sort distribution of the 3 colors plot

% Make first limit color for example with tempvar

Limit1 = find(tempvar(:,2)==1);  % find red initial value
L1 = tempvar(:,1);
index = find(Limit1==1);   % May be multiple indexes, possibly
yDesired1 = L1(Limit1);
Limit_initial = yDesired1(end,1);

% Make second limit color for example with tempvar

Limit2 = find(tempvar(:,2)==2); % find green initial value
L2 = tempvar(:,1); % the same values that L1 but written as L2
index = find(Limit1==2);   % May be multiple indexes, possibly
yDesired2 = L2(Limit2);
Limit_final = yDesired2(end,1);


%% Find a correct plot for low, intermedia, and high velocity

figure(555)
plotFluxes(qtotalrealiz,Limit_initial,Limit_final)  %Limit_initial,Limit_final
view([0 90])
% colorbar('YTicklabel',...
% {'Low Velocity','Intermediate Velocity','High Velocity'})
% % hcb = colorbar;
% set(hcb, 'Ytick'); %// 4 yticks, each "in the middle" of one color
% set(hcb, 'YTickLabel', {'low','medium','high'});


%% Create plot K vs Velocity
%velocity

for i = 1
   for j = 1
   velocity_25 = qtotal(i,j,:,:);  %the first velocity grid for 25% high K
   velocity_25 = squeeze(velocity_25);
   velocity_25 = reshape(velocity_25,1,size(velocity_25,1)*size(velocity_25,2));
   velocity_25 = velocity_25';
   end 
end

for i = 2
   for j = 1
   velocity_50 = qtotal(i,j,:,:);  %the first velocity grid for 50% high K
   velocity_50 = squeeze(velocity_50);
   velocity_50 = reshape(velocity_50,1,size(velocity_50,1)*size(velocity_50,2));
   velocity_50 = velocity_50';
   end
end

for i =3
   for j = 1 
   velocity_75 = qtotal(i,j,:,:);  %the first velocity grid for 75% high K
   velocity_75 = squeeze(velocity_75);
   velocity_75 = reshape(velocity_75,1,size(velocity_75,1)*size(velocity_75,2));
   velocity_75 = velocity_75';
   end
end

% Kvalues
for i = 1
   for j = 1
   Kvalues_25 = holdKgrid(i,j,:,:);
   Kvalues_25 = squeeze(Kvalues_25);
   Kvalues_25 = reshape(Kvalues_25,1,size(Kvalues_25,1)*size(Kvalues_25,2));
   Kvalues_25 = Kvalues_25';
   end
end

for i =2
   for j = 1
   Kvalues_50 = holdKgrid(i,j,:,:);
   Kvalues_50 = squeeze(Kvalues_50);
   Kvalues_50 = reshape(Kvalues_50,1,size(Kvalues_50,1)*size(Kvalues_50,2));
   Kvalues_50 = Kvalues_50';
   end
   end
   
for i = 3
   for j = 1
   Kvalues_75 = holdKgrid(i,j,:,:);
   Kvalues_75 = squeeze(Kvalues_75);
   Kvalues_75 = reshape(Kvalues_75,1,size(Kvalues_75,1)*size(Kvalues_75,2));
   Kvalues_75 = Kvalues_75';
   end
end

%Plot velocity vs Kvalues as dot

figure(777)
subplot(3,1,1)
plot(Kvalues_25,velocity_25,'.r');
hold on
ylabel('velocity');
title('Velocity vs Kvalues');
legend('velocity at 25% highK',4);
subplot(3,1,2)
plot(Kvalues_50,velocity_50,'.b');
ylabel('velocity');
title('Velocity vs Kvalues');
legend('velocity at 50% highK',4);
hold on
subplot(3,1,3)
plot(Kvalues_75,velocity_75,'.g');
ylabel('velocity');
title('Velocity vs Kvalues');
legend('velocity at 75% highK',4);
hold off

%histogram of v showing the results for both k values in the same axes

figure('Color',[1 1 1]);

subplot(3,1,1)
RandColor = rand(2,3);
X_25 = velocity_25(Kvalues_25==0.01);
XX_25 = velocity_25(Kvalues_25==1);
hist(X_25, 50);
h = findobj(gca,'Type','patch');
set(h,'FaceColor', RandColor(1,:),'EdgeColor',RandColor(1,:));
alpha(0.3);
hold on
hist(XX_25, 20);
h = findobj(gca,'Type','patch');
set(h,'FaceColor', RandColor(2,:),'EdgeColor',RandColor(2,:));
alpha(0.3);
xlabel('velocity of Klow and Khigh');
ylabel('Frequency');
title('Frequency of Velocity ffor Klow and Khigh (25%)');
legend('velocity at 25% highK',4);
hold on;

subplot(3,1,2)
RandColor = rand(2,3);
X_50 = velocity_50(Kvalues_50==0.01);
XX_50 = velocity_50(Kvalues_50==1);
hist(X_50, 50);
h = findobj(gca,'Type','patch');
set(h,'FaceColor', RandColor(1,:),'EdgeColor',RandColor(1,:));
alpha(0.3);
hold on
hist(XX_50, 20);
h = findobj(gca,'Type','patch');
set(h,'FaceColor', RandColor(2,:),'EdgeColor',RandColor(2,:));
alpha(0.3);
xlabel('velocity of Klow and Khigh');
ylabel('Frequency');
title('Frequency of Velocity for Klow and Khigh (50%)');
legend('velocity at 50% highK',4);
hold on;

subplot(3,1,3)
RandColor = rand(2,3);
X_75 = velocity_75(Kvalues_75==0.01);
XX_75 = velocity_75(Kvalues_75==1);
hist(X_75, 50);
h = findobj(gca,'Type','patch');
set(h,'FaceColor', RandColor(1,:),'EdgeColor',RandColor(1,:));
alpha(0.3);
hold on
hist(XX_75, 20);
h = findobj(gca,'Type','patch');
set(h,'FaceColor', RandColor(2,:),'EdgeColor',RandColor(2,:));
alpha(0.3);
xlabel('Velocity of Klow and Khigh');
ylabel('Frequency');
title('Frequency of Velocity for Klow and Khigh (75%)');
legend('velocity at 75% highK',4);
hold off;

%  myhist(XX,XX2)


hist(velocity_25,20);
hold all
hist(Kvalues_25,100);
set(get(gca,'child'),'FaceColor','none','EdgeColor','r');
xlabel('Kvalues_25');
ylabel('velocity_25');
axis([0 2 0 1])
title('Histogram of velocity for each Kvalue','FontWeight','bold','FontSize',14,...
        'FontName','Times New Roman');

    
%% High K over the ratio between Klow and Khigh for all the cases (25%, 50%, 75%)

kHigh                       = 1;
kLow                        = 0.01;  
fractionhighK               = [0.25, 0.5, 0.75]; % the example has to be 25%,50% and 75%
fractionhighK               = fractionhighK';

ratio_high_low              = [kLow./kHigh,kLow./kHigh,kLow./kHigh]; 

HighKplot                   = [fractionhighK, ratio_high_low];

figure(560)
plot(fractionhighK, ratio_high_low,'ro');


%% Generate new limits for qtotal_Min and qtotal_Max

tempvar1=reshape(qtotal_Min1(:,1),1,size(qtotal_Min1(:,1),1)*size(qtotal_Min1(:,1),2)); % turn qs into a vector
% hist(tempvar,3);
% [A B]=hist(tempvar,3)

% use kmeans clustering to define three groups
opts = statset('Display','final');
[cidx, ctrs] = kmeans(tempvar1', 3, 'Distance','city','Replicates',5,...
'Options',opts);

tempvar1(2,:)=cidx;
tempvar1=sort(tempvar1',1);
xvals=1:size(tempvar1,1);

figure(600)
colorstr='ryg';

for ii=1:size(tempvar1,1)
    plotstr=strcat('.',+colorstr(tempvar1(ii,2)));
    plot(xvals(ii),tempvar1(ii,1),plotstr);
    hold on;
end;

%find limit of the tempvar sort distribution of the 3 colors plot

% Make first limit color for example with tempvar

Limit1_ = find(tempvar1(:,2)==1);  % find red initial value
L1_= tempvar1(:,1);
index_ = find(Limit1_==1);   % May be multiple indexes, possibly
yDesired1_ = L1_(Limit1_);
Limit_initial_ = yDesired1_(end,1);

% Make second limit color for example with tempvar

Limit2_ = find(tempvar1(:,2)==2); % find green initial value
L2_ = tempvar1(:,1); % the same values that L1 but written as L2
index_ = find(Limit1==2);   % May be multiple indexes, possibly
yDesired2_ = L2_(Limit2_);
Limit_final_ = yDesired2_(end,1);
