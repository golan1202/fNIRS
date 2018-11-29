% Name: Golan Hammer
% ID: 021953237

function fNIRS_project 

load('220.mat');
channels= OxyBlockData(:,1:24);
Time=OxyBlockData(:,25);

% choose 2 fNIRS channels to compare
numChannel =1;
numChannel2 =3;


% plotting original numChannel
figure;
subplot(4,1,1);
plot(Time,channels(:,numChannel));


title(['original signal ',num2str(numChannel)]);
% Filter Heart rate and breathing from the signals
filteredChannels= zeros(length(Time),24);
for i=1:24
    filteredChannels(:,i)=fHeartBreath(channels(:,i));
end

% plotting numChannel after Filter Heart rate and breathing
subplot(4,1,2);
plot(Time,filteredChannels(:,numChannel));
title('after Filtering Heart rate and breathing');

% Detrand the signals using DCT
[w,f] = frequencyVector(Time);

afterHPF=zeros(length(Time),24);
for i=1:24
    afterHPF(:,i)=myHPF(w,filteredChannels(:,i),0.0001);
end


% plotting numChannel after detrand using DCT
subplot(4,1,3);
plot(Time,afterHPF(:,numChannel));
title('after detrand using DCT');

% AMAF to all channels 
order=10;
windowSize=15;
afterAMAF=zeros(length(Time),24);
for i=1:24
afterAMAF(:,i) = AMAF(afterHPF(:,i),order,windowSize);  
end
subplot(4,1,4);
plot(Time,afterAMAF(:,numChannel));
title('after AMAF');

% show the diffrence between the two channels
figure;
plot(Time,afterAMAF(:,numChannel),Time,afterAMAF(:,numChannel2));
title(['channel ',num2str(numChannel),' & channel ',num2str(numChannel2),' after filtering']);

% Running Correlation on each pair of the filtered signals 
vMovCorr= zeros(length(Time),24);
medMovCorr= zeros(24,24);


k=500;  % window size
for i=1:24
    for t=1:24
        vMovCorr(:,t) = MovCorr(afterAMAF(:,i),afterAMAF(:,t),k);
        medMovCorr(i,t)= median(vMovCorr(:,t),'omitnan');
    end
end

% plotting the median correlation matrix - heatmap
figure;
h = heatmap(medMovCorr);
title('median correlation matrix - heatmap');

% plotting the median correlation matrix - jet color map 
figure;
n = (24*24-22)/2;
rgb = ind2rgb(gray2ind(abs(medMovCorr),n),jet(n));
colormap jet;
image(rgb);
title('median correlation matrix - jet colormap');
end 


% -------------------------------------------------------------------------


% Filter Heart rate and breathing from the signals , using Band Stop
% Filter(BSF)
function filteredSignal = fHeartBreath(channel)

% cleanning heartbeat frequency from 0.9 to 1.6 Hz 
% and cleanning breathing frequency from  0.15 to 0.35 Hz
[b,a] = butter(5,[0.15/5 0.35/5],'stop');
channelf = filter(b,a,channel');

[b,a] = butter(5,[0.9/5 1.6/5],'stop');
channelff = filter(b,a,channelf);

filteredSignal=channelff;
end


function [w,f] = frequencyVector(t)

L = length(t);
if mod(length(t),2)
    x = linspace(-L/2, L/2, L);
else
    x = linspace(-L/2 - 1, L/2, L);
end
f = x / (L * median(diff(t)));
w = 2*pi*f;

end



function sig = myHPF(w, sig, c)
f = fftshift(dct(sig));
inds = find( w > -c & w < c);
f(inds) = 0;
sig = idct(ifftshift(f));
end



% AMAF function 
function sig = AMAF(f,order,windowSize)
f=f';
inflexionPmin=islocalmin(f);
inflexionPmax=islocalmax(f);

inflexion=zeros(1,length(f));
for i=1:length(f)
  if ((inflexionPmin(i)==1)||(inflexionPmax(i)==1))
      inflexion(i)=1;
  end
end

% count the extremum points in the window size 
% devide the signal into sections. Each section size is "windowSize"

Mat=zeros(round(length(f)/windowSize),windowSize);
Matcount=zeros(round(length(f)/windowSize),windowSize);
fpad=f(length(f))*ones(1,windowSize-rem(length(f),windowSize));
Mpad=inflexion(length(f))*ones(1,windowSize-rem(length(f),windowSize));
fnew=[f,fpad];
M=[inflexion,Mpad];
S=zeros(1,round(length(f)/windowSize));
j=1;
for i=1:round(length(f)/windowSize)
       Mat(i,:)=fnew(j:j+windowSize-1);
       Matcount(i,:)=M(j:j+windowSize-1);
       S(i)=sum(Matcount(i,:));
        j=j+windowSize;      
end
% convert S to a new scala between 1 - 2*order+1
windowSizes=rescale(S,1,2*order+1);

%padding each section in symetric way 
Matpad=zeros(round(length(f)/windowSize),windowSize+2*order);
for i=1:round(length(f)/windowSize)
    order=round(windowSizes(i)/2-0.5);
    Matpad(i,1:order)=ones(1,order)*Mat(i,1);
    Matpad(i,1+order:order+windowSize)=Mat(i,:);
    Matpad(i,1+order+windowSize:2*order+windowSize)=ones(1,order)*Mat(i,windowSize);
    
end
%clear the Matrix from elements equal to zero 
colsWithZeros = any(Matpad==0);
R = Matpad(:, ~colsWithZeros);


% moving average filter for each section accordding to neccesary
smoothedSignal=zeros(round(length(f)/windowSize),windowSize);
for i = 1 : length(windowSizes)
    smoothedSignal(i,:) = movmean(R(i,:),round(windowSizes(i)));
end
% returning the Matrix into a vector
Res=zeros(1,length(M));
j=1;
for i=1:round(length(f)/windowSize)   
        Res(j:j+windowSize-1)=smoothedSignal(i,:);
        j=j+windowSize; 
end
Res=Res(1:length(f));
sig=Res';
end



% Running Correlation function
function Cor = MovCorr(Data1,Data2,k)
y = zscore(Data2);
n = size(y,1);

if (n<k)
    Cor = NaN(n,1);
else
    x = zscore(Data1);
    x2 = x.^2;
    y2 = y.^2;
    xy = x .* y;
    A=1;
    B = ones(1,k);
    Stdx = sqrt((filter(B,A,x2) - (filter(B,A,x).^2)*(1/k))/(k-1));
    Stdy = sqrt((filter(B,A,y2) - (filter(B,A,y).^2)*(1/k))/(k-1));
    Cor = (filter(B,A,xy) - filter(B,A,x).*filter(B,A,y)/k)./((k-1)*Stdx.*Stdy);
    Cor(1:(k-1)) = NaN;
end
end