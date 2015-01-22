 function [best_strehl,best_kW]=res_batch()
 %close all
 clear;
 
Nb_courbes = 1;
diametre=[8 16 32 40];    

lambda = 1.654; 
pixsize = 0.3412; % pixsize a Shannon (= y_wfs(1:).lambda/(y_tel.diam/y_wfs(1:).nxsub)/2)

% bruit_detecteur_erms = [...
%     1.604000e-01 1.436600e+00 2.025300e+00 2.477900e+00 5.149800e+00 6.848400e+00 ...
%     8.202500e+00 1.646460e+01 2.179190e+01 2.605160e+01 5.212210e+01 6.895460e+01 ];
% 
% 
% bruit_gamma_px2 = [1.000e-04 4.000e-04 7.000e-04  1.000e-03 4.000e-03 7.000e-03 ...
%                    1.000e-02 4.000e-02 7.000e-02  1.000e-01 4.000e-01 7.000e-01 ];

bruit_detecteur_erms = cell(length(diametre),1);
bruit_alpha_arcsec2 = cell(length(diametre),1);
bruit_gamma_px2 = cell(length(diametre),1);
bruit_alpha_rad2 = cell(length(diametre),1);
bruit_phi_rad2 = cell(length(diametre),1);

bruit_detecteur_erms{1} = [ 0.16 1.02 1.66 2.48 4.04 5.77 8.20 13.01 18.41 26.05 41.20 58.28 82.42];% 130.3];

 bruit_detecteur_erms{2} = [ 0.241 0.535 0.986 1.45 2.1 3.36 4.76 6.75 10.68 15.11 21.38];% 33.80];
 
  bruit_detecteur_erms{3} = [ 0.173 0.3213 0.4997 0.8302 1.1924 1.6991 2.6985 3.8219 5.4090];% 8.5562]';
 
 bruit_detecteur_erms{4} = [0.049 0.1805 0.3049 0.5235 0.7589 1.0861 1.7229 2.4515 3.4709];% 5.4918]';



bruit_gamma_px2{1} = [ 1.000e-04 2.5000e-04 5.000e-04 ...
    1.000e-03 2.5000e-03 5.000e-03 ...
    1.000e-02 2.5000e-02 5.000e-02 ...
    1.000e-01 2.5000e-01 5.000e-01 ...
    1.000e-00];% 2.5000e-00];

 bruit_gamma_px2{2} = [5.000e-04 ... 
     1.000e-03 2.5000e-03 5.000e-03 ...
     1.000e-02 2.5000e-02 5.000e-02 ...
     1.000e-01 2.5000e-01 5.000e-01 ...
     1.000e-00];% 2.5000e-00];
 
 
 bruit_gamma_px2{3} = [2.5000e-03 5.000e-03 ...
   1.000e-02 2.5000e-02 5.000e-02 ...
   1.000e-01 2.5000e-01 5.000e-01 ...
   1.000e-00];% 2.5000e-00]';
 
bruit_gamma_px2{4} = [ 2.5000e-03 5.000e-03 ...
    1.000e-02 2.5000e-02 5.000e-02 ...
    1.000e-01 2.5000e-01 5.000e-01 ...
    1.000e-00 ];%2.5000e-00]';

%  bruit_gamma_px2{4} = [ 2.5000e-03            ...
%      1.000e-02            5.000e-02 ...
%                2.5000e-01           ...
%      1.000e-00 ];%2.5000e-00]';


for dd=1:length(diametre)
    bruit_alpha_arcsec2{dd} = bruit_gamma_px2{dd} * pixsize^2;
    bruit_alpha_rad2{dd} = bruit_alpha_arcsec2{dd} * (pi/3600/180)^2;
    bruit_phi_rad2{dd} = bruit_gamma_px2{dd} * pi^2;
end
A1=load('resultats_scripts_kalman.dat');
best_strehl = cell(length(diametre),1);
best_kW = cell(length(diametre),1); 
for dd=1:length(diametre)
    best_strehl{dd} = zeros(length(bruit_detecteur_erms{dd}),1);
    best_kW{dd} = zeros(length(bruit_detecteur_erms{dd}),1);
end


for dd=1:length(diametre)
    for bb=1:length(bruit_detecteur_erms{dd})
        B_tmp1=A1(A1(:,1)==diametre(dd) & A1(:,2)==bruit_detecteur_erms{dd}(bb),:);
        if ~isempty(B_tmp1)
            best_strehl{dd}(bb) = max(B_tmp1(:,4));
            best_kW{dd}(bb) = B_tmp1(B_tmp1(:,4)==max(B_tmp1(:,4)),3);
        end
    end
end

col='kbgrmy';
figure
hold on ; pl3=zeros(length(diametre),1)-1;
for dd=1:length(diametre),if ~isempty(best_strehl{dd}),pl3(dd)=plot(bruit_phi_rad2{dd},best_strehl{dd},['.-' col(dd)]);end,end 
title('bruit de mesure = bruit de detecteur + bruit de photons');
xlabel('bruit de mesure ASO sur le dephasage (rad^2)');ylabel('rapport de Strehl');
set(gca,'XScale','log')
set(pl3,'MarkerSize',20);
grid
legende = cell(1,length(diametre));
for dd=1:length(diametre),legende{dd}=['Kalman ASO diffractif ' num2str(diametre(dd)) 'm'];end
legend(legende,'Location','SouthWest')

figure
hold on ; pl4=zeros(length(diametre),1)-1;
for dd=1:length(diametre),if ~isempty(best_kW{dd}),pl4(dd)=plot(bruit_phi_rad2{dd},best_kW{dd},['.-' col(dd)]);end,end
title('bruit de mesure = bruit de detecteur + bruit de photons');
xlabel('bruit de mesure ASO sur le dephasage (rad^2)');ylabel('kW');
set(gca,'XScale','log');set(gca,'YScale','log')
set(pl4,'MarkerSize',20);
grid
legende = cell(1,length(diametre));
for dd=1:length(diametre),legende{dd}=['Kalman ASO diffractif ' num2str(diametre(dd)) 'm'];end
legend(legende,'Location','SouthWest')

 
