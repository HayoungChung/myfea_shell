bs = textread('bpts_sens_test.txt');

bs_ry = bs;
bs_ry(:,2) = -bs(:,2)+20;

bs = sortrows(bs,[2,1,3]);
bs_ry = sortrows(bs_ry,[2,1,3]);

dels = abs((bs - bs_ry)./bs);

figure(10); clf;
scatter(bs(:,1),bs(:,2),40,bs(:,3),'filled')
figure(11); clf;
scatter(bs_ry(:,1),bs_ry(:,2),40,bs_ry(:,3),'filled')
figure(12); clf;
scatter(bs(:,1),bs(:,2),40,dels(:,3),'filled')

figure; %% a y-axis symmtry is broken
plot(bs_ry(:,2),bs_ry(:,3),'o');
hold on
plot(bs(:,2),bs(:,3),'r+');

