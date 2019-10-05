clear
clc

foil = 'foils/SG6040.txt'; %Airfoil being used
foil2 = 'foils/SG6043.txt'; %Second airfoil
airfoilCrossoverPoint = 0.50; %Change to airfoil 2 at 50% of the blade

alphaRange = 0:0.05:15; %AoA range to test
windspeed = 12; %Windspeed on turbine
TSR = 6.3; %Target TSR
Radius = 0.35; %Radius of turbine swept area, 35 cm
radiusShapeStart = 0.03; %3cm from the root start following c~1/r
endTaperAmt = 0.01; %1cm from the end of the wing planform will be truncated
angularVel = (TSR * windspeed) / Radius; %Angular vel of turbine
pts = 20;

%Using blade area = (R*Vwind)/(nB*w) from handout
bladeArea = (Radius * windspeed) / (2*angularVel);
%Using chord=k/r (A = k*ln(r1/r0) + k) determine k
planformStart = radiusShapeStart; %r0
planformEnd = Radius - endTaperAmt; %r1
areaConst = bladeArea / (log(planformEnd / planformStart) + 1); %k
chordMax = areaConst / planformStart;

%Generate windDist - a structure with points along radius with the wind
%angle and speed at each
windDist.r = 0:(Radius/pts):Radius;
for i=1:length(windDist.r)
    [angle,speed] = getWindAngleAndSpeed(windDist.r(i),angularVel,windspeed);
    windDist.angle(i) = angle;
    windDist.speed(i) = speed;
    %TODO Calculate chord lengths
    if windDist.r(i) <= planformStart
        windDist.chord(i) = chordMax;
        continue;
    end
    windDist.chord(i) = areaConst / windDist.r(i); %5 cm for now
end

disp("Area const(k): "+areaConst);
plot(windDist.r, windDist.chord);
axis equal;
waitforbuttonpress;

ReLow = (1.225 * windDist.speed(1) * windDist.chord(1)) / (1.82*10^-5); %Using viscosity 1.82e-5 and density 1.225
midIndex = ceil(airfoilCrossoverPoint * length(windDist.speed));
ReMid = (1.225 * windDist.speed(midIndex) * windDist.chord(midIndex)) / (1.82*10^-5);
ReHigh = (1.225 * windDist.speed(length(windDist.speed)) * windDist.chord(length(windDist.speed))) / (1.82*10^-5);
[alphaOptimumRoot,alphaStallRoot] = getOptimumAlpha(foil,alphaRange,ReLow,windspeed/343);
[alphaOptimumMidFoil1,alphaStallFoil1Mid] = getOptimumAlpha(foil,alphaRange,ReMid,windspeed/343);
[alphaOptimumMidFoil2,alphaStallFoil2Mid] = getOptimumAlpha(foil2,alphaRange,ReMid,windspeed/343);
alphaOptimumMid = (alphaOptimumMidFoil1 + alphaOptimumMidFoil2) / 2.0;
[alphaOptimumTip,alphaStallTip] = getOptimumAlpha(foil2,alphaRange,ReHigh,windDist.speed(length(windDist))/343);
disp("Root: Re "+ReLow+" Alpha "+alphaOptimumRoot);
disp("Mid: Re "+ReMid+" Alpha "+alphaOptimumMid);
disp("Tip: Re "+ReHigh+" Alpha "+alphaOptimumTip);

for i=1:length(windDist.r)
    %Compute Re and M for this speed
    Re = (1.225 * windDist.speed(i) * windDist.chord(i)) / (1.82*10^-5); %Using viscosity 1.82e-5 and density 1.225
    %Linearly interpolate to get optimum alpha
    if(i<=midIndex)
        windDist.alpha(i) = alphaOptimumRoot + (alphaOptimumMid - alphaOptimumRoot)*((Re-ReLow)/(ReMid-ReLow));
    else 
        windDist.alpha(i) = alphaOptimumMid + (alphaOptimumTip - alphaOptimumMid)*((Re-ReMid)/(ReHigh-ReMid));
    end
    windDist.beta(i) = windDist.angle(i) - windDist.alpha(i);
end

function [alphaOptimum,alphaStall] = getOptimumAlpha(foil,alphaTestRange,Re,M)
    %Generate polar for given airfoil
    polar = xfoil(foil,alphaTestRange,Re,M,'oper iter 500');
%     figure();
%     plot(polar.alpha,polar.CL);
%     title("Cl vs Alpha");
%     xlabel("Alpha");
%     ylabel("Cl");
    ClCd = polar.CL ./ polar.CD;
    plot(polar.alpha,ClCd);
    title("Cl/Cd Vs Alpha");
    xlabel("Alpha");
    ylabel("Cl/Cd");

    disp("Re: "+Re+", M: "+M);
    %Find optimum-ish alpha for this Re ('ish' since not using interpolation)
    [ClCdMax, maxClCdIndex] = max(ClCd);
    alphaOptimumDeg = polar.alpha(maxClCdIndex);
    disp("Optimum alpha: "+alphaOptimumDeg+"deg (Cl/Cd="+ClCdMax+")");
    alphaOptimum = deg2rad(alphaOptimumDeg);

    %Find stall angle
    [ClMax, maxClIndex] = max(polar.CL);
    alphaStallDeg = polar.alpha(maxClIndex);
    disp("Stall alpha: "+alphaStallDeg+"deg (ClMax="+ClMax+")");
    alphaStall = deg2rad(alphaOptimumDeg);
    
    disp("Press key to continue");
    waitforbuttonpress;
end

%Get local wind angle and speed at point along simple 2D blade. Angle is
%-ve for 
function [angle,speed] = getWindAngleAndSpeed(r,angularVel,windspeed)
    vel = getLocalWindVel(r,angularVel,windspeed);
    angle = atan(vel(2) / vel(1));
    speed = sqrt(vel(1)^2 + vel(2)^2);
end

%Get local wind velocity at point along simple 2D blade
function [totalVel] = getLocalWindVel(r, angularVel, windspeed)
    wakeInductionFac = 1/3; %Assume wake induction factor 1/3
    windVel = [0,windspeed*(1-wakeInductionFac)];
    apparentWindVel = [angularVel*r,0]; %Apparent wind is positive
    totalVel = windVel + apparentWindVel;
end