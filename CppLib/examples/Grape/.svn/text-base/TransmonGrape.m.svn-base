function Run()
'Start'
Run1Qubit();


%This function basically calls the Grape algorithm (GetPulse) with a bunch
%of different configuration options.  OptimizeEnvelope optimizes the 
%envelope of the pulse.  trydifferentcutoffs tries some different amplitude
%cutoffs for the controls.  It is also possible to loop through various
%pulse durations (j), initial configuration randomizations (l) (in order to
%avoid local maxima), or amplitude cut-offs (f).
function [f1, f2] = Run1Qubit()
PropagatorFinal3LevelNOT = [ 1 0 0; 0 0 1; 0 1 0];
PropagatorFinal3LevelPauliY = [ 1 0 0; 0 0 -i; 0 i 0];
PropagatorFinal3LevelX90 = [ 1 0 0; 0 sqrt(2)/2 -sqrt(2)*i/2; 0 -sqrt(2)*i/2 sqrt(2)/2];
PropagatorFinal3LevelY90 = [ 1 0 0; 0    0.5 + i/2    -0.5 - i/2; 0    0.5 + i/2   0.5 + i/2];
GateFinal = PropagatorFinal3LevelNOT;

optimizeEnvelope = true;
trydifferentcutoffs = false;
plotfidelities = true;


%Required accuracy of pulse
eps=0.0000001;

Ec=2*pi*383;
Ej=40*Ec;
g0=2*pi*100;
g1=g0*sqrt(2);
g2=g0*sqrt(3);
w01=sqrt(8*Ec*Ej);
w12=w01-Ec;
w23=w12-Ec;
w02=w01+w12;
wr=w01+6282;
d0=w01-wr;
d1=w12-wr;
d2=w23-wr;
lambda0=g0/d0;
lambda1=g1/d1;
n0=g0*g1*(w01-w02);
X0=g0*g0/d0;
X1=g1*g1/d1;
X2=g2*g2/d2;
%wd=0.999*(w01+X0)
wd=1.00*(w01+X0)
dr0=wr +0 -X0 -wd;
dr1=wr +X0 -X1 -wd;
dr2=wr +X1 -X2 -wd;

e0=1/dr0;
e1=1/dr1;
e2=1/dr2;
E1=0; %w01 - wd + X0
E2=w02 - 2*wd + X1;

LE1=w01+X0;
LE2=w02+X1;

%transmon Hamiltonian in Rotating Frame (2*2 controls)
 rwaHamiltonDrift3Level1 = [ E2 0 0; 0 E1 0; 0 0 0];
 rwaHamiltonControlX3Level = [0 lambda1 0; lambda1 0 lambda0; 0 lambda0 0];
 rwaHamiltonControlX3LevelSq = [e2 0 0; 0 e1 0; 0 0 e0];
 rwaHamiltonControlY3Level = [0 -lambda1*i 0; lambda1*i 0 -lambda0*i; 0 lambda0*i 0];
 rwaHamiltonControlY3LevelSq = [e2 0 0; 0 e1 0; 0 0 e0];
 
%transmon Lab-Frame Hamiltonian (1st and 2nd order) for 2 controls 
                      
labHamiltonDrift3Level1 = [ LE2 0 0; 0 LE1 0; 0 0 0];
labHamiltonControlX3Level1 = 2*[0 lambda1 0; lambda1 0 lambda0; 0 lambda0 0];
labHamiltonControlX3Level1Sq = 0*4*[e2 0 0;0 e1 0;0 0 e0];
labHamiltonControlY3Level1 = 2* [0 -lambda1*i 0; lambda1*i 0 -lambda0*i; 0 lambda0*i 0];
labHamiltonControlY3Level1Sq = 0*4*[e2 0 0;0 e1 0;0 0 e0];


rwaDriftHamiltonians = cell(1, 1);
rwaDriftHamiltonians{1} = rwaHamiltonDrift3Level1;
numRWAControls=2;
rwaHamiltonians = cell(1, numRWAControls);
rwaHamiltoniansSq = cell(1, numRWAControls);
rwaHamiltonians{1} = rwaHamiltonControlX3Level;
rwaHamiltoniansSq{1} = rwaHamiltonControlX3LevelSq;
rwaHamiltonians{2} = rwaHamiltonControlY3Level;
rwaHamiltoniansSq{2} = rwaHamiltonControlY3LevelSq;

labframe_Hdrift = cell(1, 1);
labframe_Hdrift{1} = labHamiltonDrift3Level1;
numLabControls=2;
labframeH = cell(1, numLabControls);
labframeH2 = cell(1, numLabControls);
labframeH{1} = labHamiltonControlX3Level1;
labframeH2{1} = labHamiltonControlX3Level1Sq;
labframeH{2} = labHamiltonControlY3Level1;
labframeH2{2} = labHamiltonControlY3Level1Sq;



%N : number of pulse durations to try (to find shortest optimal)
N = 0

if(N==0)
   plotfidelities=false; 
end

%R : number of random initial pulse shapes to try
R=1
global PenaltyStrength;
PenaltyStrength=5;

TimeSlices = 10000;
FidelityPrecision = 1e-16;
topfidelity = -100;
topfidelunpen=-1;
numcutoffs = 1 + trydifferentcutoffs*2;
bestcontrols=zeros(2,TimeSlices);
bestl=-1;

fidelities = zeros(5, N+1);  
pulselengths = zeros(1, N+1);
maxamps=zeros(numcutoffs,4);
maxamps(3,:) = [500 500 500 500];
maxamps(2,:) = [0.7 0.2 0.7 0.2];
maxamps(1,:) = [0.7 0.1 0.7 0.1];


Start = 0.004%*TimeSlices%/200 
End = Start +0.001*N;
Delta = 0.001%/300;

%cycle through possible total-times for pulse   
for j = 0:N

    j    
    time = Start + j * Delta
    
    if(optimizeEnvelope)
        numControls=numRWAControls;    
        driftH = rwaDriftHamiltonians;
        controlH = rwaHamiltonians;
        controlSqH = rwaHamiltoniansSq;
        timesperpixel=1;
    else
         numControls=numLabControls;
        driftH = labframe_Hdrift;
        controlH = labframeH;
        controlSqH = labframeH2; 
        timesperpixel=70;
        Urot = [exp(-i*2*wd*time) 0 0; 0 exp(-i*wd*time) 0; 0 0 1];
        GateFinal = Urot*GateFinal*Urot'
        return;
    end

    ntimes = timesperpixel*TimeSlices;
    times = GenerateTimeSlicesNumber(time, ntimes); 
    dt = times(1);
    %try various random initial conditions for pulse
   for l=1:R  
      l     %Gaussian
    % randcontrols1 = GetRandomControlsPalindromic(-800, 0, 2, TimeSlices);
     %randcontrols2 = GetRandomControlsPalindromic(-50, 50, 2, TimeSlices);
     randcontrols1 = GetRandomSawtooth(-2000, 100, ntimes);
     randcontrols2 = GetRandomSawtooth(-130, 130, ntimes);
     sigma=time/5;
     %ntimes=pixels
     initialControls = zeros(1,ntimes);
     gaussEnvelope = GaussPulse(time, TimeSlices, 2.5, sqrt(pi/2)/lambda0/2);
     
      for k = 1:ntimes
          if k (~mod(k-1,timesperpixel))
            prefac = gaussEnvelope(1, (k-1)/timesperpixel+1);
          end            
            initialControls(1, k) = randcontrols1(1,k);%hermitecontrols(k);%0.55; % % + initialControlstemp(k);\
%            if(numControls>1)
               initialControls(2, k) = randcontrols2(1,k);%hermitecontrols(k);%-0.3;%initialControlstemp2(k); % + initialControlstemp(k);     
%            end
            if (l==1)
                if optimizeEnvelope
                   initialControls(1, k) = prefac; % * sin(wd*dt*k);
                else 
                   initialControls(1, k) = prefac * sin(wd*dt*k); 
                end
              % initialControls(1, k) = prefac;
               initialControls(2, k) = 0; %gaussEnvelope(1,k);
            end
            %if(numControls==4)
%               initialControls(3, k) = randcontrols1(1,k);%0.54;%initialControlstemp2(k); %picontrols(l,k); % + initialControlstemp(k);
%               initialControls(4, k) = randcontrols2(1,k);%0.57;%initialControlstemp2(k); %picontrols(l,k); % + initialControlstemp(k);
%            end
      end
      pulselengths(1, j+1) = time; 
% 
  controls = initialControls;
%  bestcontrols = controls;
  propagators = GetPropagators(driftH, controlH, controlSqH,  controls, times);
  
%  bestpropagators=propagators;  
  
% [p1, p2, p3] = GetLevelPopulation1Qubit(propagators,0);
%  figure;
   [p1, p2, p3] = GetLevelPopulation1Qubit(propagators,1);
%   subplot(2,1,1);
%   PlotPopulation(times, p1, p2, p3);
%   subplot(2,1,2);
%figure();
%   PlotnControls(initialControls, numControls, times);   
[forward, backward] = GetPropagation(propagators);
%forward{1,5}
  fidel = FidelityRelativePhaseInvariantSquare(GateFinal, forward)
  
% fidelities(1, j+1) = fidel;
% continue;
    
     %compare different amplitude cut-offs for controls
    for f=1:numcutoffs
tic        
controls = GetPulse(times, initialControls, driftH, controlH, controlSqH,  GateFinal, 1-eps, FidelityPrecision, ...
          @GradientRelativePhaseInvariant, @FidelityRelativePhaseInvariantSquare, @PenaltyVoid, @PenaltyVoidGradient, maxamps(f,:));
toc
        propagators = GetPropagators(driftH, controlH, controlSqH,  controls, times);
        [forward, backward] = GetPropagation(propagators);
        fidel = FidelityRelativePhaseInvariantSquare(GateFinal, forward)  +  PenaltyPopulation3Level(GateFinal, times, forward);
         if(fidel>fidelities(f, j+1))   
            fidelities(f, j+1) = FidelityRelativePhaseInvariantSquare(GateFinal, forward);
            bestl=l;
            if(~trydifferentcutoffs)
                bestpropagators=propagators;
                bestcontrols = controls;
                bestinitialcontrols = initialControls;
                topfidelity = fidel;
                topfidelunpen = FidelityRelativePhaseInvariantSquare(GateFinal, forward);
            end
         end           
    end      
  end
   
   csvwrite(strcat('labconXNOT',num2str(TimeSlices),'.txt'),bestcontrols(1,:));
   csvwrite(strcat('labconYNOT',num2str(TimeSlices),'.txt'),bestcontrols(2,:));
   bestl
   topfidelity
   topfidelunpen
   TimeSlices=TimeSlices+1;
end

if(trydifferentcutoffs)
figure
 plot(pulselengths(1,:), fidelities(1,:), pulselengths(1,:), fidelities(2,:), pulselengths(1,:), fidelities(3,:));  
 xlabel('t (\mus)')
 ylabel('fidelity');
elseif(plotfidelities)
 figure
 plot(pulselengths(1,:), fidelities(1,:));  
 xlabel('t (\mus)')
 ylabel('fidelity');
else
'Top Fidelity'
[p1, p2, p3] = GetLevelPopulation1Qubit(bestpropagators,0);
figure;
subplot(2,2,1);
PlotPopulation(times, p1, p2, p3);
[p1, p2, p3] = GetLevelPopulation1Qubit(bestpropagators,1);
subplot(2,2,2);
PlotPopulation(times, p1, p2, p3);
subplot(2,2,3);
PlotnControls(bestcontrols, numControls, times);
%PlotControls1(bestcontrols(2,:), times)
subplot(2,2,4);
%PlotnControls(bestinitialcontrols, numControls, times);
PlotPulseSpectrum(bestcontrols, numControls, ntimes, time);
topfidelity
topfidelunpen
end 


return; %Testing consistency of Rotating Wave Approx
    %Reconstruct lab-frame Hamiltonian
    if(optimizeEnvelope)
        labcontrols =  zeros(numLabControls, ntimes);
            for k = 1:ntimes
            labcontrols(1,k) = controls(1,k) + picontrols(1,k) + controls(3,k);
            end
            for k = 1:ntimes
            labcontrols(2,k) = controls(2,k) + picontrols(1,k) + controls(4,k);
            end
        PlotControls2(labcontrols, times);   
        propagators = GetPropagators(labframe_Hdrift, labframeH, labcontrols, times);
        [p1, p2, p3] = GetLevelPopulation1Qubit(propagators);
        PlotPopulation(times, p1, p2, p3);

        %if fidelity4 > 0.9999
        %    'found fault tolerant'
        %end
        %dlmwrite('f4_8.txt', fidelities4, 'precision', 16);
    else    
    %Construct RWA wave from lab frame H
        rwcontrols =  zeros(numRWAControls, ntimes);
            for k = 1:ntimes
                rwcontrols(1,k) = controls(1,k) - c_u*picontrols(1,k);
                rwcontrols(3,k) = controls(1,k) - c_u*picontrols(1,k);
                rwcontrols(2,k) = controls(2,k) - c_phi*picontrols(1,k);
                rwcontrols(4,k) = controls(2,k) - c_phi*picontrols(1,k);
            end
        PlotControls4(rwcontrols, times);   
        figure;       
        PlotPulseSpectrum(rwcontrols, 4, ntimes, time); 
        propagators = GetPropagators(rwaDriftHamiltonians, rwaHamiltonians, rwcontrols, times);
        [p1, p2, p3] = GetLevelPopulation1Qubit(propagators);
        PlotPopulation(times, p1, p2, p3);        
    
    %can this pulse be re-optimized in new RWA frame?
%     controls = GetPulse(times, rwcontrols, rwaDriftHamiltonians, rwaHamiltonians, GateFinal, ...
%         1-eps, FidelityPrecision, @GradientRelativePhaseInvariant, @FidelityRelativePhaseInvariantSquare, @PenaltyVoid);
%      PlotControls4(controls, times); 
%      
%         propagators = GetPropagators(rwaDriftHamiltonians, rwaHamiltonians, controls, times);
%         [p1, p2, p3] = GetLevelPopulation1Qubit(propagators);        
%         PlotPopulation(times, p1, p2, p3); 
%        figure;
%        PlotPulseSpectrum(controls, numControls, ntimes, time)
%        

    end

%returns a random smooth palindromic pulse with amplitude between (from,to)
function controls = GetRandomControlsPalindromic(from, to, numberOfControls, numberOfTimes)
numberOfTimesHalf = numberOfTimes/2;
initialXHalf = 1:15:numberOfTimesHalf;
initialYHalf = rand(numberOfControls, size(initialXHalf, 2));
initialXHalfSize = size(initialXHalf, 2);
initialX = zeros(1, initialXHalfSize*2);
initialX(1, 1:initialXHalfSize) = initialXHalf(1, :);
for k = 1:initialXHalfSize
    initialX(1, 2*initialXHalfSize - k + 1) = 2*numberOfTimesHalf - initialXHalf(1, k);
end
initialY = zeros(1, initialXHalfSize*2);
initialY(1, 1:initialXHalfSize) = initialYHalf(1, :);
for k = 1:initialXHalfSize
    initialY(1, 2*initialXHalfSize - k + 1) = initialYHalf(1, k);
end
controls = spline(initialX, initialY, 1:numberOfTimes);
controls = from + (to - from)*controls;   

function controls = GetRandomSawtooth(from, to, numberOfTimes)
width = to-from;
controls = zeros(1, numberOfTimes);
for k = 1:numberOfTimes
    controls(1,k)=width*rand()+from;
end

%returns Rabi pulse at frequency d01 but capped at the max amplitudes
function controls = CappedInitialControls(amp_u, amp_phi, time, ntimes, d01)
times = GenerateTimeSlicesNumber(time, ntimes);
period = round(2*ntimes/time/d01);
t = zeros(1, size(times, 2));
prefac = 1;
controls = zeros(2, size(times, 2));
for a = 1:size(times, 2)
     if(mod(a,period)==0) 
         prefac = -prefac;
     end
     controls(2,a) = prefac.*amp_u;
     controls(1,a) = amp_phi;
end    

function controls = GaussPulse(time, ntimes, alpha, amplitude)
N = ntimes;
controls = zeros(1, N);

Start = -time/2;   
End = time/2;
Delta = (End - Start) / N;

for k = 1:N
    t = Start + k * Delta - Delta/2;
    controls(1, k) = ((2*alpha*amplitude)/time)*exp(-2*alpha^2*(t/time)^2);
end

function controls = GaussPulseDerivativeAmplitude(time, times, alpha, amplitude)
N = times;
controls = zeros(1, N);

Start = -time/2;
End = time/2;
Delta = (End - Start) / N;

for k = 1:N
    t = Start + k * Delta - Delta/2;
    controls(1, k) = (2*alpha/time)*exp(-2*alpha^2*(t/time)^2);
end

function controls = GaussPulseDerivativeAlpha(time, times, alpha, amplitude)
N = times;
controls = zeros(1, N);

Start = -time/2;
End = time/2;
Delta = (End - Start) / N;

for k = 1:N
    t = Start + k * Delta - Delta/2;
    controls(1, k) = (2*amplitude/time)*exp(-2*alpha^2*(t/time)^2) + (-2*2*alpha*(t/time)^2)*((2*alpha*amplitude)/time)*exp(-2*alpha^2*(t/time)^2);
end


%This is essentially the GRAPE algorithm.  The inputs are the Hamiltonians,
%the method for calculating the gradient, and the penalty function. The
%steepest ascent algorithm climbs to the targetFidelity, or when the
%controls come nearest to the intended propagatorFinal
function returnedControls = GetPulse(times, initialControls, hamiltonDrift, hamiltonControl, hamiltonControl2, propagatorFinal, targetFidelity, fidelityPrecision, gradientMethod, fidelityMethod, penaltyMethod, penaltyGradientMethod, maxampls)
N = size(times, 2);
numberOfControls = size(hamiltonControl, 2);
controls = initialControls;

tracker = 0;
length = 0;
propagators = GetPropagators(hamiltonDrift, hamiltonControl, hamiltonControl2, controls, times);
[propagationForward, propagationBackward] = GetPropagation(propagators);
currentFidelityUnpenalized = fidelityMethod(propagatorFinal, propagationForward)
%penalty = penaltyMethod(controls, maxampls);
%penalty2= PenaltyPopulation3Level(propagatorFinal, times, propagationForward);
%currentFidelityUnpenalized = currentFidelity - penalty
currentFidelity = currentFidelityUnpenalized %+ penalty2
ntriespengrad=0;

while 1
    tracker = tracker + 1;
    if ~mod(tracker, 100)
        tracker
        newFidelity
        newFidelityUnpenalized
    end
    %penalty = penaltyMethod(hamiltonDrift, hamiltonControl, propagatorFinal, controls, times);
    %penalty = PenaltySharpRiseAndFall2(controls, times, propagators);
    %penalty = PenaltyPopulation3Level(hamiltonDrift, hamiltonControl, propagatorFinal, controls, times);
    gradient = gradientMethod(hamiltonControl, hamiltonControl2, propagatorFinal, controls, times, propagationForward, propagationBackward);
   % penaltyGradient_2 = PenaltyPopulation3LevelGradient(hamiltonControl, hamiltonControl2, propagatorFinal, controls, times, propagators, propagationForward);
    %gradient = gradient + (penaltyGradient_2);
   % penaltyGradient = penaltyGradientMethod(controls, gradient, maxampls);   
    %gradient = gradient + penaltyGradient;
    while 1
        
        for penaltyStrength = 0:ntriespengrad
        
            modifiedGradient = (2^length) *2000*(gradient);%+(ntriespengrad-penaltyStrength+1)*penaltyGradient_2);
            propagators = GetPropagators(hamiltonDrift, hamiltonControl, hamiltonControl2, controls + modifiedGradient, times);
            [propagationForward, propagationBackward] = GetPropagation(propagators);
            newFidelityUnpenalized = fidelityMethod(propagatorFinal, propagationForward) ;
            %penalty = penaltyMethod(controls + modifiedGradient, maxampls);
            penalty2= 0; %PenaltyPopulation3Level(propagatorFinal, times, propagationForward);
            %newFidelity = newFidelityUnpenalized - penalty;
            newFidelity = newFidelityUnpenalized + penalty2;

            if newFidelity >1
               'Grape failed'
                returnedControls = controls;
                currentFidelity
                tracker
                return; 
            end
            if newFidelity >= targetFidelity
                returnedControls = controls + modifiedGradient;
                'Fidelity'
                tracker
                newFidelity
                newFidelityUnpenalized
                return;
            end

            if newFidelity > (currentFidelity + fidelityPrecision)
                break;
            end        
        
        end 
        
        if abs(newFidelity - currentFidelity) < fidelityPrecision
            returnedControls = controls;
            'Precision'
            fidelityPrecision
            currentFidelity
            newFidelityUnpenalized
            length
            if(length>-5)
                'bad'
            end
            tracker
            return;
        end
        
        if newFidelity > currentFidelity
            length = length + 1;
            controls = controls + modifiedGradient;
            currentFidelity=newFidelity;    
            break;
        end
        
        length = length - 1;
                
        if length < -52
            'Grape failed'
            returnedControls = controls;
            currentFidelity
            tracker
            return;
        end
        if length > 52
            length = 52;
        end
        
    end
end
returnedControls = controls;

function gradient = GradientBlock(Perturbed, Unperturbed)
gradient = -2*real(trace(Perturbed)*trace(Unperturbed));

%Calculates fidelity to 3-Level Gate while ignoring phase of level 3
function fidelity = FidelityRelativePhaseInvariantSquare(propagatorFinal, propagationForward)
N = size(propagationForward, 2);
dimension = size(propagatorFinal, 1);
M = propagatorFinal'*propagationForward{N};
block2Level = M(2:3, 2:3);
overlap1 = abs(trace(block2Level));
%overlap2 = abs(M(1, 1));
overlap2 = 0;
fidelity = (overlap1^2 + overlap2^2)/4.0;

%Calculates the gradient f the 3-Level propagators while ignoring the phase
%aquired by the third level (since it will not be populated)
%Calculates the gradient f the 3-Level propagators while ignoring the phase
%aquired by the third level (since it will not be populated)
function gradient = GradientRelativePhaseInvariant(hamiltonControl, hamiltonControl2, propagatorFinal, controls, times, propagationForward, propagationBackward)
N = size(times, 2);
numberOfControls = size(hamiltonControl, 2);
gradient = zeros(numberOfControls, N);
for j = 1:N
    for k = 1:numberOfControls
        Perturbed = i*times(1, j)*(hamiltonControl{k}+2*controls(k,j)*hamiltonControl2{k})*propagationForward{j};
        Unperturbed = propagationForward{j};
        Perturbed = propagatorFinal' * propagationBackward{j} * Perturbed;
        Unperturbed = propagatorFinal' * propagationBackward{j} * Unperturbed; 
        UnperturbedDagger = Unperturbed';

        Perturbed2Level = Perturbed(2:3, 2:3);
        Unperturbed2Level = UnperturbedDagger(2:3, 2:3);
        grad = GradientBlock(Perturbed2Level, Unperturbed2Level);
        gradient(k, j) = grad;
    end
end

%Dummy penalty funtions
function penalty = PenaltyVoid(controls, times, propagators, maxampls)
penalty = 0;

function penalty = PenaltyVoidGradient(controls, gradient, times, propagators, maxampls)
penalty = 0;

%Penalty function forbidding amplitudes larger than maxampls
function penalty = PenaltyLargeAmplitude(controls, maxampls)
numControls = size(controls, 1);
numberOfTimes = size(controls, 2);
penalty = 0;

for j = 1:numControls 
for k = 1:numberOfTimes
    penalty = penalty + (controls(j, k)>maxampls(j)) + (controls(j, k)<-maxampls(j)); %((controls(j, k)/maxampls(j)).^514);    
end
end

%gradient of max-amplitude penalty function
function penaltygradient = PenaltyGradientLargeAmplitude(controls, gradient, maxampls)
newcontrols = controls + gradient;
numControls = size(controls, 1);
numberOfTimes = size(controls, 2);
penaltygradient = zeros(numControls, numberOfTimes);

for j = 1:numControls 
for k = 1:numberOfTimes
    penaltygradient(j, k) = (newcontrols(j, k)>maxampls(j))*(maxampls(j)-newcontrols(j, k) +(controls(j,k)-maxampls(j))/1.3) + ...
                            (newcontrols(j, k)<-maxampls(j))*(-maxampls(j)-newcontrols(j, k) +(-controls(j,k)-maxampls(j))/1.3);
                    %514*((controls(j, k)/maxampls(j)).^513)/maxampls(j);    
end
end

%Penalized population in the third level during pulse duration (gradient)
function penalty = PenaltyPopulation3LevelGradient(hamiltonControl, hamiltonControl2, propagatorFinal, controls, times, propagators, forward)
numControls = size(controls, 1);
numberOfTimes = size(controls, 2);
gamma = 500;
penalty = zeros(numControls, numberOfTimes);
for j = 1:numControls 
for k = 1:numberOfTimes-1
    one = hamiltonControl{j} * forward{k} + 2*hamiltonControl2{j} * forward{k}*controls(j,k);  
for l = k:numberOfTimes-1  
    Perturbed = i * times(1, k) * one(1, 1);
    two = forward{l};
    Unperturbed = two(1, 1);
    penalty(j, k) = penalty(j, k) + GradientBlock(Perturbed, Unperturbed');
    one= propagators{l+1}*one;
end
    penalty(j, k) = gamma * penalty(j, k)*times(1,k);
end
end

%Penalized population in the third level during pulse duration 
function penalty = PenaltyPopulation3Level(propagatorFinal, times, forward)
numberOfTimes = size(forward, 2);
gamma = 500;
penalty = 0;
for k = 1:numberOfTimes
    one = forward{k};
    two = forward{k};
    Perturbed = one(1, 1);
    Unperturbed = two(1, 1);
    penalty = penalty +abs(Perturbed * Unperturbed')*times(1,k);
end
penalty=gamma *(penalty - numberOfTimes*times(1,1));

%Given Hamiltonian, calculates propagator = e^iHt, for all time steps
function propagators = GetPropagators(hamiltonDrift, hamiltonControl, hamiltonControl2, controls, times)
N = size(controls, 2);       %number of times
propagators = cell(1, N); 
dimension = size(hamiltonDrift{1}, 1);
Identity = eye(dimension);
for g = 1:N
    hamilton = 0*Identity;
    numberOfDrifts = size(hamiltonDrift, 2);
    for k = 1:numberOfDrifts
        hamilton = hamilton + hamiltonDrift{k};
    end
    numberOfControls = size(hamiltonControl, 2);
    for k = 1:numberOfControls
        hamilton = hamilton + controls(k, g) * hamiltonControl{k} + controls(k, g) * controls(k, g) * hamiltonControl2{k};
    end
    propagators{1, g} = expm(-i * times(1, g) * hamilton);
end

%Splits Gate into 2 Gates: before and after t-th step, for each step t
function [forward, backward] = GetPropagation(propagators)
N = size(propagators, 2);
dimension = size(propagators{1}, 1);
forward = cell(1, N);
backward = cell(1, N);
Identity = eye(dimension);
propagationForward = Identity;
for j = 1:N
    propagationForward = propagators{j} * propagationForward;
    forward{1, j} = propagationForward;
end
propagationBackward = Identity;
for j = 1:N
    backward{1, N-j+1} = propagationBackward;
    propagationBackward = propagationBackward * propagators{N-j+1};
end

%splits time into equal intervals
function times = GenerateTimeSlicesNumber(time, N)
trotterizeTime = time/N;
times = trotterizeTime * ones(1, N);

%obtains qutrit populations over span of pulse
function [population1, population2, population3] = GetLevelPopulation1Qubit(propagators, initiallevel)
N = size(propagators, 2);
population1 = zeros(1, N+1);
population2 = zeros(1, N+1);
population3 = zeros(1, N+1);
propagator = eye(3);
level1 = [ 0 0 1 ];
level2 = [ 0 1 0 ];
level3 = [ 1 0 0 ];
startlevel = (initiallevel)*level2 + (~initiallevel)*level1;
for k = 1:N+1
    step = propagator * startlevel';
    population1(1, k) = abs(level1 * step)^2;
    population2(1, k) = abs(level2 * step)^2;
    population3(1, k) = abs(level3 * step)^2;
    if k ~= N+1
        propagator = propagators{k} * propagator;
    end
end

%
%Plotting functions for displaying pulseinformation graphically
%
function PlotPopulation(times, population1, population2, population3)
px = zeros(1, size(population1, 2)); 
for a = 1:size(px, 2)
    if a == 1
        px(1, a) = 0;
    else
        px(1, a) = px(1, a-1) + times(1, a-1);
    end
end
%figure;
plot(px, population1, px, population2, px, population3, px, 1 - population1 - population2 - population3)
xlabel('t (\mus)')
ylabel('Population');

function PlotnControls(controls, ncontrols, times)
if(ncontrols==4)
    PlotControls4(controls, times);
else if (ncontrols==2)
       PlotControls2(controls, times); 
    else
      PlotControls1(controls, times); 
    end
end

function PlotPulseSpectrum(controls, ncontrols, ntimes, totalTime)

nfreqs = ntimes;    
fs = 2*pi*ntimes/totalTime;
freqs = [(0.5-nfreqs/2):(nfreqs/2-0.5)]*fs/nfreqs;
for i=1:ncontrols
   %subplot(ncontrols,1,i);
   FT = fftshift(abs(fft(controls(i,:), nfreqs)));
   %subplot(2,2,3+i);
   plot(freqs, FT);
   xlabel('f (MHz)');
   xlim([-fs/16, fs/16]);
   %if(fs/3<90) 
   %set(gca,'XTick',[-fs/8:7.5:90]);
   %end
end


function PlotControls1(controls, times)
cx = zeros(1, size(times, 2)); 
for a = 1:size(cx, 2)
    if a == 1
        cx(1, a) = times(1, a);
    else
        cx(1, a) = cx(1, a-1) + times(1, a);
    end
end
%figure;
plot(cx, controls(1,:))
xlabel('t (\mus)');
ylabel('\lambda');

function PlotControls2(controls, times)
cx = zeros(1, size(times, 2)); 
for a = 1:size(cx, 2)
    if a == 1
        cx(1, a) = times(1, a);
    else
        cx(1, a) = cx(1, a-1) + times(1, a);
    end
end
%figure;
plot(cx, controls(1,:), '-',  cx, controls(2,:), '--')
xlabel('t (\mus)');
ylabel('\lambda ');

function PlotControls4(controls, times)
cx = zeros(1, size(times, 2)); 
for a = 1:size(cx, 2)
    if a == 1
        cx(1, a) = times(1, a);
    else
        cx(1, a) = cx(1, a-1) + times(1, a);
    end
end
%figure;
plot(cx, controls(1,:), '-', cx, controls(2,:), '--', cx, controls(3,:), '.', cx, controls(4,:), '*')
xlabel('t (\mus)');
ylabel('\lambda');


