%% now just plot voxel of choice with fit
root = load_root;
vdir = 'vRFs';
fitdir = 'vRFfits';
sub = 'AR'; voi = 'V3AB';
vfn = sprintf('%s%s/%s_ridgeVRF_CVBilat_%s.mat',root,vdir,sub,voi);
load(vfn);
ffn = sprintf('%s%s/%s_%s_GridFit_ThreshVRFs_LambdaMinBIC_ALLSESS.mat',...
    root,fitdir,sub,voi);
load(ffn,'bfpar');

basis_size = FWHM * 1/rad2fwhm(1);
xx = chanX * (basis_size/2);
yy = -1*(chanY * (basis_size/2));

%%
vv = 32;

figure;
colormap hot;

ps(1) = subplot(2,1,1);  hold all;
imagesc(xx,yy,vrf_mat{3}(:,:,vv),[-0.6 0.9]);  hold on; axis manual;
cx = rad2fwhm(bfpar{3}(vv,3))*cos(linspace(0,2*pi,1001)) + bfpar{3}(vv,1);
cy = rad2fwhm(bfpar{3}(vv,3))*sin(linspace(0,2*pi,1001)) + (-1*bfpar{3}(vv,2));
plot(bfpar{3}(vv,1),-1*bfpar{3}(vv,2),'ko');
plot(cx,cy,'w--','LineWidth',2);
xlim([-5.5 5.5]); ylim([-3 3]);
axis off; axis equal; clear cx cy

ps(2) = subplot(2,1,2);  hold all;
imagesc(xx,yy,vrf_mat{1}(:,:,vv),[-0.6 0.9]); axis manual;
cx = rad2fwhm(bfpar{1}(vv,3))*cos(linspace(0,2*pi,1001)) + bfpar{1}(vv,1);
cy = rad2fwhm(bfpar{1}(vv,3))*sin(linspace(0,2*pi,1001)) + (-1*bfpar{1}(vv,2));
plot(bfpar{1}(vv,1),-1*bfpar{1}(vv,2),'ko');
plot(cx,cy,'w--','LineWidth',2);
xlim([-5.5 5.5]); ylim([-3 3]);
axis off; axis equal; clear cx cy

clear ps;
suptitle('AR V3AB vox 32 (fix v left)');