clear all;
close all;

I1=load -ascii('cameraman.mat');

I3=load -ascii('cameraman_flou.mat');


I3image=reshape(I3.y,256,256);
figure();
subplot(1,2,1)
imagesc(I1.x);
colormap('gray');
title('image originale');
%subplot(1,3,2)
%imagesc(I2image);
subplot(1,2,2)
imagesc(I3image);
colormap('gray');
title('image dégradée');

figure();

I4=load('flou.mat');
%subplot(1,2,1)
imagesc(I4.H);
colormap('gray');
title('floue oubruit ???');


figure()
err= I3image - real(ifft2(fft2(I1.x).*I4.H));
% subplot(1,2,1)
imagesc(err);
colormap('gray');
% subplot(1,2,2)
figure()
err_reshape=reshape(err,65536,1);
% h=hist(err_reshape,100);
% sigma1=std(err_reshape)
% bar(h); %pic Ã  1/(sigma*(sqrt(2pi)))
% title('histogramme du bruit');
% sigma=1/(96*sqrt(2*pi)) %sigma proche de 0 donc centrÃ©
ecartType_buit=std(err_reshape(:));
esperance_bruit=mean(err_reshape(:));

nbrClasse=100;
hist(err_reshape,nbrClasse);

title({['\bfHistogramme de la reponse frequentielle du filtre pour ',num2str(nbrClasse),' classes'],
    ['Esperance : ',num2str(esperance_bruit)],
    ['Ecart-type : ',num2str(ecartType_buit)]});
set(gcf,'color','white')
print(3,'-dmeta','f3');
 
 

Ilap=load('laplacien.mat');
%% Question 3
% Reconstitution de l'image d'origine
lambda=0.5;
I4H=I4.H; %image du flou
ILG=Ilap.G;

Ireconstruit = real(ifft2((1./(conj(I4H).*I4H + lambda*conj(ILG).*ILG )).*(conj(I4H)) .*fft2(I3image)));

figure();
subplot(1,2,1)
imagesc(Ireconstruit);
colormap('gray');
title({['\bfReconstruction par la methode'],
    ['de Tikhonov pour \lambda=',num2str(lambda)]});

subplot(1,2,2)
imagesc(I3image);
title('\bfImage degradee');


erreur=norm(I1.x-Ireconstruit);
display(['Erreur de reconstuction pour lamba=',num2str(lambda),' : ',num2str(erreur)]);
 
figure()

% Evolution de l'erreur de restauration en fonction du parametre lambda
lambdaMin=0.001;
lambdaMax=0.5;
pasLambda=0.001;
 
lambdaTest=lambdaMin:pasLambda:lambdaMax;
erreur=zeros(1,max(size(lambdaTest)));
for i=1:max(size(lambdaTest))
    IrecontruitTest=real(ifft2((1./(conj(I4H).*I4H + lambdaTest(i)*conj(ILG).*ILG )).*(conj(I4H)) .*fft2(I3image)));
    erreur(i)=norm(I1.x-IrecontruitTest);
end
 
plot(lambdaTest,erreur,'LineSmoothing','on');
title({['\bfEvolution de l''erreur de restauration'],['pour \lambda variant de ',num2str(lambdaMin),' a ',num2str(lambdaMax),' avec un pas de ',num2str(pasLambda)]});
xlabel('\lambda');
ylabel('Erreur de restauration');
set(gcf,'color','white')
grid on
print(5,'-dmeta','f5');


%%%%%%%%%%%%ù
n=100;
yI=fft2(I3image); % image floutÃ©
x_n= I3image;

beta = lambda + 1;
Gamma = 2/(beta + 2);

for k=0:1:n
  x_n=x_n-Gamma.*(conj(I4H).*(I4H.*x_n-yI) +(lambda .*conj(ILG).*ILG.*x_n));
  
  F_n=0.5*norm(I4H.*x_n - yI).^2 +( lambda .* norm(ILG.*x_n).^2);  
    
  err = norm(x_n-I1.x).^2;

  
end

figure()

imagesc(real(ifft2(x_n)));
colormap('gray');
title('Image x_n');


figure()
plot(F_n)
title('F_n');

figure()
plot(err)
title('erreur');
%% Question 4
cameramanCS=load('cameraman_cs.mat');
Decimation=load('decimation.mat');
D=Decimation.D;
CS_R = cameramanCS.y;
D = reshape(D,256,256);

figure();

imagesc(D);
colormap('gray');
title('decimation');

tI=256;

Z=zeros(tI,tI);


ind=find(D==1);
Z(ind)=CS_R;


Z=reshape(Z,tI,tI);


figure();

imagesc(Z);
colormap('gray');
title('Image avec decimÃ©e');

x_nD= I3image;
for k=0:1:n
  x_nD=x_nD-Gamma.*(conj(Z).*(Z.*x_n-yI) +(lambda .*conj(ILG).*ILG.*x_nD));
    
end
figure()

imagesc(real(ifft2(x_nD))),colormap('gray');
title('Image x_nD');