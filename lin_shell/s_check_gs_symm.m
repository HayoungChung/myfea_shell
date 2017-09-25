gs = textread('gpts_sens_test.txt');

gs_ry = gs;
gs_ry(:,2) = -gs(:,2)+20;

gs = sortrows(gs,[2,1,3]);
gs_ry = sortrows(gs_ry,[2,1,3]);

dels = abs((gs - gs_ry)./gs);

figure(10); clf;
scatter(gs(:,1),gs(:,2),40,gs(:,3),'filled')
figure(11); clf;
scatter(gs_ry(:,1),gs_ry(:,2),40,gs_ry(:,3),'filled')
figure(12); clf;
scatter(gs(:,1),gs(:,2),40,dels(:,3),'filled')

figure; %% a y-axis symmtry is broken
plot(gs_ry(:,2),gs_ry(:,3),'o');
hold on
plot(gs(:,2),gs(:,3),'r+');

