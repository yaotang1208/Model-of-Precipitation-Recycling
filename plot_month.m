%==============================================================%
%==============Plot Monthly Average Recycling Ratio============%
%==============================================================%

    clear
    clc
    load RR_ncep03.out            
    load U1_ncep03.out             
    load V1_ncep03.out            
    load E1_ncep03.out            
    load W1_ncep03.out            

    IY=41;                      % Central US IY=33, JX=41
    JX=61;
    n0 = 0;
    ntime = 239; 

  t1 = n0+1;             % start time of display 
  t2 = n0+ntime;      % end time (4xdaily) 
  dn = ntime; 

  UM = zeros(IY,JX);
  VM = zeros(IY,JX);
  EM = zeros(IY,JX);
  WM = zeros(IY,JX);
  RM = zeros(IY,JX);

  figure, orient landscape

  for t=t1:1:t2    
    y1=(t-1)*IY+1;
    y2=t*IY;
    UM=UM+U1_ncep03(y1:1:y2,:)/dn;
    VM=VM+V1_ncep03(y1:1:y2,:)/dn;
    EM=EM+E1_ncep03(y1:1:y2,:)/dn;
    WM=WM+W1_ncep03(y1:1:y2,:)/dn;
    RM=RM+RR_ncep03(y1:1:y2,:)/dn;
        
  end
  
      rvar = num2str(mean(mean(100*RM))); 
%       xlon = (360-105):0.25:(360-95);
      xlon = -100:0.25:-85;
      ylat = 32.5:0.25:42.5;
  
     hr=[1 2 3 4 5 6 7 8 9 10];
    [csr, hr]=contourf(xlon,ylat,100*RM);
    clabel(csr,hr), 
%     grid, 
    colorbar,
    caxis([0 100])
    hold on
    quiver(xlon,ylat,UM,VM), 
    xlabel('longitude')
    ylabel('latitude')
    axis('square')
    title(strcat('Monthly Recycling Ratio (%) ,June 1998, <\rho>_A=',rvar,'%'))
    hold off

%     subplot(223),
%     he=[100 200 300 400 500 600];
%     clabel(contour(xlon_ncep03-360,ylat_ncep03,EM.*WM,he));
%     hold on
%     grid
%     xlabel('longitude')
%     ylabel('latitude')
%     title(strcat('Evaporation (W m^{-2})'))
%     contour(xlon_ncep03-360,ylat_ncep03,land_ncep03,[1],'k');
%     hold off
% 
%     subplot(224),
%     clabel(contour(xlon_ncep03-360,ylat_ncep03,WM));
%     hold on
%     grid
%     xlabel('longitude')
%     ylabel('latitude')
%     title(strcat('Precipitable Water (kg m^{-2})'))
%     contour(xlon_ncep03-360,ylat_ncep03,land_ncep03,[1],'k');
%     hold off
