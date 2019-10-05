clear
clc

disp(pwd);

%Mechanical properties of material
yieldStress = 26000000; %In WEAKEST direction (for safety)
E = 2180000; %Young's modulus in WEAKEST direction
flexuralE = 1760000000;%1760000000; %Mathematically same as E, but for bending according to wikipedia, value from material datasheet

%Specify airfoil at root
foils.startPoint(1) = 0;
foils.file(1) = {'foils/SG6040.txt'};
foils.startPoint(2) = 0.60;%0.40; %0.35
foils.file(2) = {'foils/SG6040.txt'};
foils.file(3) = {'foils/FX63-137.txt'};
foils.startPoint(3) = 0.80;
%Specify airfoil from 50% outwards
foils.startPoint(4) = 0.95;%0.90; %0.80
foils.file(4) = {'foils/SG6043.txt'};
numFoilPts = 100; %Number of airfoil sections to test in x foil to interpolate between, higher is more accurate but slower

materialDensity = 1040; %ABS M30 density
alphaRange = 0:0.05:15; %AoA range to test in degrees
windspeed = 12; %Windspeed on turbine
TSR = 6.3; %Target TSR
Radius = 0.30; %Radius of turbine swept area, 35 cm
radiusShapeStart = 0.065;%0.06;%0.03; %3cm from the root start following c~1/r
endTaperAmt = 0.01; %1cm from the end of the wing planform will be truncated
angularVel = (TSR * windspeed) / Radius; %Angular vel of turbine
rpm = angularVel * 9.55;
disp("RPM: "+rpm);
pts = 60;
minAngleFromStallDeg = 1;
hubDiameter = 0.03; %chord length that blade blends into, 4cm
bladeStartRadius = 0.04 / 2; %first d/2 is hub, diameter is 4.5 cm

%Using blade area = (R*Vwind)/(nB*w) from handout
bladeArea = (Radius * windspeed) / (2*angularVel);
%Using chord=k/r (A = k*ln(r1/r0) + k) determine k
planformStart = radiusShapeStart+0.03; %r0 - add 0.02 to try and increase chord lengths a tad to reduce flex
planformEnd = Radius; %r1
areaConst = bladeArea / (log(planformEnd / planformStart));%{ + 1}%); %k
areaConst = areaConst * 1.06; %Add 6% so that chords are a bit better structurally
chordMax = areaConst / planformStart;
minAllowedChord = 0.022; %2.2 cm min chord

%Generate windDist - a structure with points along radius with the wind
%angle and speed at each
windDist.r = 0:(Radius/pts):Radius;
for i=1:length(windDist.r)
    [angle,speed] = getWindAngleAndSpeed(windDist.r(i),angularVel,windspeed);
    windDist.angle(i) = angle;
    windDist.speed(i) = speed;
    %Calculate chord lengths
    if windDist.r(i) <= planformStart
        if(windDist.r(i) <= bladeStartRadius)
            windDist.chord(i) = hubDiameter;
        else
            windDist.chord(i) = hubDiameter + (chordMax-hubDiameter)*((windDist.r(i)-bladeStartRadius)/planformStart-bladeStartRadius);
        end
        continue;
    end
    windDist.chord(i) = areaConst / windDist.r(i);
    if windDist.chord(i) < minAllowedChord
        windDist.chord(i) = minAllowedChord;
    end
end

disp("Area const(k): "+areaConst);
plot(windDist.r, windDist.chord);
title('Chord length vs Radius');
axis equal;
%waitforbuttonpress;

ReMin = 9999999;
ReMax = 0;
for i=1:length(windDist.r)
    %Compute Re and M for this speed
    Re = (1.225 * windDist.speed(i) * windDist.chord(i)) / (1.82*10^-5); %Using viscosity 1.82e-5 and density 1.225
    ReDist(i) = Re;
    if(Re < ReMin)
        ReMin = Re;
    end
    if(Re > ReMax)
        ReMax = Re;
    end
end

figure();
plot(windDist.r,ReDist);
title('Re vs Radius');

ReTestStep = 10000;
ReTestMin = floor(ReMin / (1000 * 10.0)) * (1000 * 10);
ReTestMax = ceil(ReMax / (1000 * 10.0)) * (1000 * 10);
disp("Min Re: "+ReMin+", Max Re: "+ReMax+"; Test min Re: "+ReTestMin+", Test max Re: "+ReTestMax+", Re Test step: "+ReTestStep+" - Tune test step if jump is too large");
%waitforbuttonpress;

disp("Calculating blended foil coords for structural analysis...");
foilData = getInterpFoilCoords(Radius,pts+1,foils);
approxMass = 0;
for ri=1:length(windDist.r)
   r = windDist.r(ri);
   foilAtPt = foilData.foil(ri);
   structure.radius(ri) = r;
   structure.foil(ri) = foilAtPt;
   %Remove random line from trailing edge that xfoil seems to generate
    for fx=1:length(foilAtPt.x)
        %Find the end coordinate for the airfoil in the file (1,0)
       if fx~=1 && foilAtPt.x(fx) == 1 && foilAtPt.y(fx) == 0
          %End coordinate for the airfoil
          endIndex = fx;
          foilX = foilAtPt.x(1:endIndex);
          foilY = foilAtPt.y(1:endIndex);
       end
    end
    foilX = foilX.*windDist.chord(ri);
    foilY = foilY.*windDist.chord(ri);
    structure.area(ri) = polyarea(foilX,foilY);
    structure.sectionThickness = Radius / pts;
    structure.approxVol(ri) = structure.area(ri)*structure.sectionThickness;
    if(r>=bladeStartRadius)
       approxMass = approxMass + structure.approxVol(ri)*materialDensity; 
    end
end
disp("Approx mass: "+approxMass+"Kg");

%Get airfoil properties along the span
ReToTest = ReTestMin:ReTestStep:ReTestMax; %Approx Re across the airfoil
foilDataMap = getInterpPolars(Radius,numFoilPts,foils,alphaRange,ReToTest,windspeed/343);

objAoaFig = figure();
ClAlphaFig = figure();

for i=1:length(windDist.r)
    %Compute Re and M for this speed
    Re = (1.225 * windDist.speed(i) * windDist.chord(i)) / (1.82*10^-5); %Using viscosity 1.82e-5 and density 1.225
    %Get airfoil data at this point
    airfPolarAtPoint = getFoilPolarFromInterp(windDist.r(i),foilDataMap,Re);
    windDist.airfPolar(i) = airfPolarAtPoint;
    curMax = -Inf; %Current value of objective function we are optimising
    curAlphaDeg = NaN;
    curAlphaIndex = 0;
    secondaryAlphaDeg = NaN;
    secondaryAlphaIndex = 0;
    alphaStallDeg = rad2deg(airfPolarAtPoint.alphaStallRad);
    obj = zeros(1,length(alphaRange));
    CLAtPt = 0;
    CDAtPt = 0;
    for ai=1:length(alphaRange) %Crudely iterate over alpha as this is sufficient
       alphaTest = alphaRange(ai);
       %Polar indexes for alpha are the same as for alphaRange
       %Want max obj=Clsin(O)-Dcos(O) where O is the wind angle, 90deg at
       %root and 0 at tip (Although program uses radians)
       obj(ai) = objectiveFunction(windDist.r(i),Radius,airfPolarAtPoint.CL(ai),airfPolarAtPoint.CD(ai),windDist.angle(i));
       if(alphaStallDeg - minAngleFromStallDeg < alphaTest) %If airfoil is near stalling
           break; %Stop looking
       end
       if(obj(ai) > curMax)
           CLAtPt = airfPolarAtPoint.CL(ai);
           CDAtPt = airfPolarAtPoint.CD(ai);
           curMax = obj(ai);
           %Function seems to generally increase with aoa, until it doesn't
           if curAlphaIndex < length(alphaRange)
               secondaryAlphaIndex = curAlphaIndex+1;
               secondaryAlphaDeg = alphaRange(secondaryAlphaIndex);
           end
           curAlphaDeg = alphaTest;
           curAlphaIndex = ai;
       end
    end
    
    figure(objAoaFig);
    plot(alphaRange,obj);
    title("Obj vs AoA, r/R= "+(windDist.r(i)/Radius)+", wind angle= "+windDist.angle(i));
    figure(ClAlphaFig);
    plot(alphaRange,airfPolarAtPoint.CL,'b');
    hold on
    plot(alphaRange,airfPolarAtPoint.CD,'r');
    hold off
    title("Cl,Cd vs AoA, r/R= "+(windDist.r(i)/Radius)+", wind angle= "+windDist.angle(i));
    
    if(~isnan(secondaryAlphaDeg)) %If found two points with high val of function
        %Iterate between curAlphaDeg and secondaryAlphaDeg and linearly
        %interpolate to find best alpha (Assumes best val of function between
        %best and second best valued points)
        if (secondaryAlphaDeg > curAlphaDeg)
            swap = secondaryAlphaDeg;
            secondaryAlphaDeg = curAlphaDeg;
            curAlphaDeg = swap;
            swap = secondaryAlphaIndex;
            secondaryAlphaIndex = curAlphaIndex;
            curAlphaIndex = swap;
        end
        alphasToTest = secondaryAlphaDeg:0.005:curAlphaDeg; %Alphas to test over
        CL1 = airfPolarAtPoint.CL(secondaryAlphaIndex);
        CL2 = airfPolarAtPoint.CL(curAlphaIndex);
        CD1 = airfPolarAtPoint.CD(secondaryAlphaIndex);
        CD2 = airfPolarAtPoint.CD(curAlphaIndex);
        for jj=1:length(alphasToTest) %Linearly interpolate at each alpha and eval obj
           CLTest = CL1 + ((CL2-CL1) * ((alphasToTest(jj)-secondaryAlphaDeg)/(curAlphaDeg-secondaryAlphaDeg)));
           CDTest = CD1 + ((CD2-CD1) * ((alphasToTest(jj)-secondaryAlphaDeg)/(curAlphaDeg-secondaryAlphaDeg)));
           obj2 = objectiveFunction(windDist.r(i),Radius,CLTest,CDTest,windDist.angle(i));
           if(obj2 > curMax) %Update our best angle
               CLAtPt = CLTest;
               CDAtPt = CDTest;
               curMax = obj2;
               curAlphaDeg = alphasToTest(jj);
           end
        end
    end
    disp("Final opt alpha "+curAlphaDeg+", Stall "+alphaStallDeg);
    %pause(0.05);
    windDist.CLAtOpt(i) = CLAtPt;
    windDist.CDAtOpt(i) = CDAtPt;
    windDist.alpha(i) = deg2rad(curAlphaDeg);
    windDist.beta(i) = windDist.angle(i) - windDist.alpha(i);
end

%Check at 6m/s that everything is okay
windspeed = 6;
for i=1:length(windDist.r)
    [angle,speed] = getWindAngleAndSpeed(windDist.r(i),angularVel,windspeed);
    windDist2.angle(i) = angle;
    windDist2.speed(i) = speed;
end
ReMin2 = 9999999;
ReMax2 = 0;
for i=1:length(windDist.r)
    %Compute Re and M for this speed
    Re = (1.225 * windDist2.speed(i) * windDist.chord(i)) / (1.82*10^-5); %Using viscosity 1.82e-5 and density 1.225
    ReDist2(i) = Re;
    if(Re < ReMin2)
        ReMin2 = Re;
    end
    if(Re > ReMax2)
        ReMax2 = Re;
    end
end
for i=1:length(windDist.r)
    %Compute Re and M for this speed
    Re = (1.225 * windDist2.speed(i) * windDist.chord(i)) / (1.82*10^-5); %Using viscosity 1.82e-5 and density 1.225
    %Get airfoil data at this point
    airfPolarAtPoint = getFoilPolarFromInterp(windDist.r(i),foilDataMap,Re);
    curMax = -Inf; %Current value of objective function we are optimising
    curAlphaDeg = NaN;
    curAlphaIndex = 0;
    secondaryAlphaDeg = NaN;
    secondaryAlphaIndex = 0;
    alphaStallDeg = rad2deg(airfPolarAtPoint.alphaStallRad);
    aoa = windDist.alpha(i);
    if aoa > alphaStallDeg
       disp("Blade STALLS at 6m/s at r="+windDist.r(i)+", stall at "+alphaStallDeg+", AoA is "+aoa); 
    elseif aoa > alphaStallDeg-minAngleFromStallDeg
       disp("Blade NEAR STALL at 6m/s at r="+windDist.r(i)+", stall at "+alphaStallDeg+", AoA is "+aoa); 
    end
end

figure();
plot(windDist.r,rad2deg(windDist.alpha));
title('AoA vs Radius');

figure();
plot(windDist.r,rad2deg(windDist.beta));
title('Twist angle vs Radius');

%Structures
windspeed = 12;
rotFoilPlot = figure();
for ri=1:length(windDist.r)
   r = windDist.r(ri);
   foilAtPt = foilData.foil(ri);
   %Remove random line from trailing edge that xfoil seems to generate
    for fx=1:length(foilAtPt.x)
        %Find the end coordinate for the airfoil in the file (1,0)
       if fx~=1 && foilAtPt.x(fx) == 1 && foilAtPt.y(fx) == 0
          %End coordinate for the airfoil
          endIndex = fx;
          foilX = foilAtPt.x(1:endIndex);
          foilY = foilAtPt.y(1:endIndex);
       end
    end
    %Rotate by alpha
    v = [transpose(foilX);transpose(foilY)]; %Merge x and y data into a vector
    v = v.*windDist.chord(ri); %make the correct scale for the chord length
    theta = windDist.beta(ri); %Angle to rotate by
    R = [cos(theta) -sin(theta); sin(theta) cos(theta)]; %Rotation matrix
    v2 = R*v; %Apply rotation, about origin
    foilX = v2(1,:); %Read out x from rotated coords
    foilY = v2(2,:); %Read out y from rotated coords
    [ geom, iner, ~ ] = polygeom( foilX, foilY );
    x_cen = geom(2);
    y_cen = geom(3);
    %foilX = foilX - x_cen;
    %foilY = foilY - y_cen;
    structure.foilX(ri) = {foilX};
    structure.foilY(ri) = {foilY};
    structure.IxxCentroid(ri) = iner(4);
    structure.IyyCentroid(ri) = iner(5);
    %Align coords to centroid
    foilX = foilX - x_cen;
    foilY = foilY - y_cen;
    structure.critYDist(ri) = max(abs(min(foilY)), abs(max(foilY))); %Distance further away from centroid in -y direction
    structure.critXDist(ri) = max(abs(min(foilX)), abs(max(foilX))); %Distance further away from centroid in -y direction
    figure(rotFoilPlot);
    plot(foilX,foilY);
    axis equal;
    pause(0.05);
end

%Export airfoils as pts files
for ri=1:length(structure.radius)
    fName = "bladePts\bladeCoords"+ri+".pts";
    file = fopen(fName,'w');
    foilX = cell2mat(structure.foilX(ri));
    foilY = cell2mat(structure.foilY(ri));
    for fx=1:length(foilX)
       fprintf(file,['%.6f %.6f 0',newline],foilX(fx),foilY(fx)); 
    end
    fclose(file);
end

% disp("Calculaing equiibrium rotation rate...");
% angularVelEquil = 0;
% windspeed = 12;
% turningMomentClosestToZero = 9999999999;
% angularVels = angularVel:5:500;
% for avi=1:length(angularVels) %Test angular velocities from designed to 500
%     angularVel = angularVels(avi);
%     %Update wind angles and speeds at each point
%     for wi=1:length(windDist.r)
%         [angle,speed] = getWindAngleAndSpeed(windDist.r(wi),angularVel,windspeed);
%         windDist3(avi).angle(wi) = angle;
%         windDist3(avi).speed(wi) = speed;
%     end
%     for zi=1:length(windDist.r)
%         %Compute Re and M for this speed
%         Re = (1.225 * windDist3(avi).speed(zi) * windDist.chord(zi)) / (1.82*10^-5); %Using viscosity 1.82e-5 and density 1.225
%         %Get airfoil data at this point
%         loadingTest(avi).airfPolar(zi) = getFoilPolarFromInterp(windDist.r(zi),foilDataMap,Re);
%         aoa = rad2deg(windDist3(avi).angle(zi)-windDist.beta(zi)); %Convert to degrees as airf data is in degrees
%         for alphaIndex=1:length(loadingTest(avi).airfPolar(zi).alpha)
%            alphaTest = loadingTest(avi).airfPolar(zi).alpha(alphaIndex);
%            if alphaIndex < length(loadingTest(avi).airfPolar(zi).alpha) %Not the last elem
%                alphaIndex2 = alphaIndex+1;
%                alphaTest2 = loadingTest(avi).airfPolar(zi).alpha(alphaIndex+1);
%            else
%                alphaIndex2 = alphaIndex;
%                alphaTest2 = alphaTest;
%            end
%            if(alphaTest <= aoa && alphaTest2>=aoa)
%               %Linearly interpolate to get Cd and Cl
%               CL1 = loadingTest(avi).airfPolar(zi).CL(alphaIndex);
%               CL2 = loadingTest(avi).airfPolar(zi).CL(alphaIndex2);
%               CL = CL1 + (CL2-CL1)*((aoa-alphaTest)/(alphaTest2-alphaTest));
%               CD1 = loadingTest(avi).airfPolar(zi).CD(alphaIndex);
%               CD2 = loadingTest(avi).airfPolar(zi).CD(alphaIndex2);
%               CD = CD1 + (CD2-CD1)*((aoa-alphaTest)/(alphaTest2-alphaTest));
%               loadingTest(avi).CL(zi) = CL;
%               loadingTest(avi).CD(zi) = CD;
%               A = windDist.chord(zi)*structure.sectionThickness; %Ref area
%                 L = 0.5 * 1.225 * windDist3(avi).speed(zi)^2 * A * CL;
%                 D = 0.5 * 1.225 * windDist3(avi).speed(zi)^2 * A * CD;
%                 %windDist.angle is the wind angle, which lift and drag are defined
%                 %relative to
%                 loadingTest(avi).rotForce(zi) = L*sin(windDist3(avi).angle(zi)) - D*cos(windDist3(avi).angle(zi));
%                break; %Stop iterating over alphas
%            end
%         end
%     end
%     %Have CL,CD,L,D,rotForce at each point, calculate
%     %turning moment
%     Mtestrot = 0; %M in plane due to force rotating blade
%    for ri2=1:length(windDist.r)
%       Mtestrot = Mtestrot + (loadingTest(avi).rotForce(ri2) * abs(windDist.r(ri2)));
%    end
%    loadingTestData.RPM(avi) = (angularVel / 0.104719755);
%    loadingTestData.M(avi) = Mtestrot;
%    disp("RPM: "+(angularVel / 0.104719755)+", M="+Mtestrot);
%    if abs(Mtestrot) < turningMomentClosestToZero
%        turningMomentClosestToZero = Mtestrot;
%        angularVelEquil = angularVel;
%    end
% end
% 
% maxAngularVel = angularVelEquil;
% maxRpm = angularVelEquil / 0.104719755;
% 
% disp("Blade at equilibrium at "+maxRpm+"rpm +- 50");

maxRpm = 3000;
maxAngularVel = maxRpm * 0.104719755;
%TODO Structural analysis
%Centrifugal forces
centriF = 0;
disp("Centrifugal stresses:");
for ri=length(structure.radius):-1:1
    sectionArea = structure.area(ri);
    dm = sectionArea * structure.sectionThickness * materialDensity;
    centriF = centriF + dm*maxAngularVel^2*structure.radius(ri);
    loading.radius(ri) = structure.radius(ri);
    loading.centrifugalForce(ri) = centriF;
    loading.centrifugalStress(ri) = centriF / sectionArea;
    stressPercent = (loading.centrifugalStress(ri)*100) / yieldStress;
    disp("r="+loading.radius(ri)+", stress="+stressPercent+"%, force="+centriF);
end

%Forces due to lift and drag
for ri=1:length(structure.radius)
    %Calculate force at each point
    %Force rotating around is max when at opt. tsr and angular vel
    %Force out of plane assume to be max at same time for ease (even if not
    %quite true)
    CL = windDist.CLAtOpt(ri);
    CD = windDist.CDAtOpt(ri);
    A = windDist.chord(ri)*structure.sectionThickness; %Ref area
    L = 0.5 * 1.225 * windDist.speed(ri)^2 * A * CL;
    D = 0.5 * 1.225 * windDist.speed(ri)^2 * A * CD;
    %windDist.angle is the wind angle, which lift and drag are defined
    %relative to
    loading.rotForce(ri) = L*sin(windDist.angle(ri)) - D*cos(windDist.angle(ri));
    loading.transverseForce(ri) = L*cos(windDist.angle(ri)) + D*sin(windDist.angle(ri));
end

%Calculate analysis of beam at each point due to rot and transverse forces
for ri=1:length(structure.radius)
   %Bending moment at point is sum of forces further along the beam mult. by each's distance from this point
   Mrot = 0; %M in plane due to force rotating blade
   Mtrans = 0; %M in plane perp. to force rotating blade
   for ri2=ri+1:length(structure.radius)
      Mrot = Mrot + loading.rotForce(ri2) * abs(loading.radius(ri2) - loading.radius(ri));
      Mtrans = Mtrans + loading.transverseForce(ri2) * abs(loading.radius(ri2) - loading.radius(ri));
   end
   loading.Munit(ri) = 1*abs(loading.radius(length(structure.radius)) - loading.radius(ri));
   loading.Mrot(ri) = Mrot;
   loading.Mtrans(ri) = Mtrans;
   loading.MOverIRot(ri) = loading.Mrot(ri) / structure.IyyCentroid(ri);
   loading.OneOverEIRot(ri) = 1 / (flexuralE*structure.IyyCentroid(ri));
   loading.MOverITransverse(ri) = loading.Mrot(ri) / structure.IxxCentroid(ri);
   loading.OneOverEITransverse(ri) = 1 / (flexuralE*structure.IxxCentroid(ri));
   loading.maxRotStress(ri) = loading.MOverIRot(ri)*structure.critXDist(ri);
   loading.maxTransverseStress(ri) = loading.MOverITransverse(ri)*structure.critYDist(ri);
   maxAllowableStressWithoutWarn = yieldStress * 0.8;
   rotStressPercent = loading.maxRotStress(ri)*100 / yieldStress;
   transStressPercent = loading.maxTransverseStress(ri)*100 / yieldStress;
   if(loading.maxRotStress(ri) > maxAllowableStressWithoutWarn)
      disp("Stress due to rot. at r="+structure.radius(ri)+" is "+loading.maxRotStress(ri)+" ("+rotStressPercent+"%)"); 
   end
   if(loading.maxTransverseStress(ri) > maxAllowableStressWithoutWarn)
      disp("Stress out of rot. plane at r="+structure.radius(ri)+" is "+loading.maxTransverseStress(ri)+" ("+transStressPercent+"%)"); 
   end
end

%Tip deflections:
for ri=1:length(structure.radius)
    %Rotations at each point
    %dTheta = (M/EI) * dz
    loading.defTrans(ri) = 0;
    loading.defRot(ri) = 0;
    for ri2=1:ri
        loading.defTrans(ri) = loading.defTrans(ri)+loading.Mtrans(ri2)*loading.Munit(ri2)*loading.OneOverEITransverse(ri2)*structure.sectionThickness;
        loading.defRot(ri) = loading.defRot(ri)+loading.Mrot(ri2)*loading.Munit(ri2)*loading.OneOverEIRot(ri2)*structure.sectionThickness;
    end
end
tipDeflectionTrans = loading.defTrans(length(structure.radius));
tipDeflectionRot = loading.defRot(length(structure.radius));
% for ri=1:length(structure.radius)
%     %Theta = dx/dz; dx = Theta * dz; x = sum of dx
%     tipDeflectionTrans = tipDeflectionTrans + loading.rotTrans(ri)*structure.sectionThickness;
%     tipDeflectionRot = tipDeflectionRot + loading.rotRot(ri)*structure.sectionThickness;
% end
disp("Tip deflection in plane of rotation: "+tipDeflectionRot+"m");
disp("Tip deflection out of plane of rotation: "+tipDeflectionTrans+"m");

%Function to maximise at each point along airfoil
function [x] = objectiveFunction(r,Radius,Cl,Cd,windAngleRad)
    x = Cl*sin(windAngleRad) - Cd*cos(windAngleRad);
    %x = Cl / Cd;
end

%Get the polar of an airfoil at a given radius along the span
function polar = getFoilPolarFromInterp(r,foilDataMap,Re)
    for Re1Index=1:length(foilDataMap.Re) %foilDataMap.Re is in ascending order
        Re2Index = Re1Index; %Second Re polar to blend with is this one
        if (foilDataMap.Re(Re1Index) > Re) %If the datapoint's Re is larger than here
           Re1Index = Re1Index - 1; %Use last datapoint as the first polar to blend
           if Re1Index < 1
               Re1Index = 1;
           end
           break; %And end loop
        end
        if Re1Index == length(foilDataMap.Re) %Reached end of Re vals, and not large enough
            Re1Index = Re1Index - 1; %Use last datapoint as the first polar to blend
            %Note extrapolation quality not guaranteed
        end
    end
    
    %Blending is necessary
    try
        blendAmt = (Re-foilDataMap.Re(Re1Index)) / (foilDataMap.Re(Re2Index) - foilDataMap.Re(Re1Index));
    catch
        blendAmt = 0; %Single airfoil, no need to blend
    end
    if(isnan(blendAmt))
        blendAmt = 0; %Single airfoil
    end
    polar1 = getFoilPolarFromSingleReInterp(r,foilDataMap.foilData(Re1Index));
    polar2 = getFoilPolarFromSingleReInterp(r,foilDataMap.foilData(Re2Index));
    %Generate blended polar
    polar.name="BlendedAirfoil";
    alphas = foilDataMap.foilData(Re1Index).alphasUsed;
    %Unfortunately xfoil omits alphas where it doesn't converge, so
    %linearly interpolating is a bit messy
    for ai=1:length(alphas) %For each alpha value
       polar.alpha(ai) = alphas(ai); %get the correct val of alpha
       pd1Index = 1; %Index of where polar1's data for this alpha is at
       pd2Index = 1; %Index of where polar2's data for this alpha is at
       
       %This basically makes non-existing alphas use the values from the
       %ones that existed before it. May produce very bad results if huge
       %gaps in alpha data present
       for d=0:ai-1 %For each alpha value from this to the ones before it
           pd1Index = ai-d;
           if(ai-d > length(polar1.alpha)) %If this alpha is outside of range successfully simulated
               continue; %Skip to next iteration where we try lower alpha
           end
           alphaP1 = polar1.alpha(ai-d);
           if alphaP1 <= polar.alpha(ai) %if this alpha is correct
               break; %We have our index, so end the for loop
           end
       end
       for d2=0:ai-1 %For each alpha value from this to the ones before it
           pd2Index = ai-d2;
           if(ai-d2 > length(polar2.alpha)) %If this alpha is outside of range successfully simulated
               continue; %Skip to next iteration where we try lower alpha
           end
           alphaP2 = polar2.alpha(ai-d2);
           if alphaP2 <= polar.alpha(ai) %if this alpha is correct
               break; %We have our index, so end the for loop
           end
       end
       polar.CL(ai) = polar1.CL(pd1Index) + blendAmt * (polar2.CL(pd2Index) - polar1.CL(pd1Index)); %Linearly interpolate
       polar.CD(ai) = polar1.CD(pd1Index) + blendAmt * (polar2.CD(pd2Index) - polar1.CD(pd1Index)); %Linearly interpolate
    end
    %Find stall angle
    [ClMax, maxClIndex] = max(polar.CL);
    alphaStallDeg = polar.alpha(maxClIndex);
    polar.alphaStallRad = deg2rad(alphaStallDeg);
    polar.Re = polar1.Re + blendAmt * (polar2.Re - polar1.Re); %Linearly interpolate
    polar.Ncrit = polar1.Ncrit + blendAmt * (polar2.Ncrit - polar1.Ncrit); %Linearly interpolate
end

%Get the polar of an airfoil at a given radius along the span
function polar = getFoilPolarFromSingleReInterp(r,foilData)
    for polarIndex1=1:length(foilData.testRad) %For each foil data point ordered by test radius
        polarIndex2 = polarIndex1; %Second polar to blend with is this one
        if (foilData.testRad(polarIndex1) > r) %If the datapoint's radius is larger than here
           polarIndex1 = polarIndex1 - 1; %Use last datapoint as the first polar to blend
           break; %And end loop
        end
    end
    
%     if(polarIndex1 == polarIndex2)
%        %No blended necessary
%        polar = foilData.polar(polarIndex1);
%        return;
%     end
    %Blending is necessary
    try
        blendAmt = (r-foilData.testRad(polarIndex1)) / (foilData.testRad(polarIndex2) - foilData.testRad(polarIndex1));
    catch
        blendAmt = 0; %Single airfoil, no need to blend
    end
    if(isnan(blendAmt))
        blendAmt = 0; %Single airfoil
    end
    polar1 = foilData.polar(polarIndex1);
    polar2 = foilData.polar(polarIndex2);
    %Generate blended polar
    polar.name="BlendedAirfoil";
    alphas = foilData.alphasUsed;
    %Unfortunately xfoil omits alphas where it doesn't converge, so
    %linearly interpolating is a bit messy
    for ai=1:length(alphas) %For each alpha value
       polar.alpha(ai) = alphas(ai); %get the correct val of alpha
       pd1Index = 1; %Index of where polar1's data for this alpha is at
       pd2Index = 1; %Index of where polar2's data for this alpha is at
       
       %This basically makes non-existing alphas use the values from the
       %ones that existed before it. May produce very bad results if huge
       %gaps in alpha data present
       for d=0:ai-1 %For each alpha value from this to the ones before it
           pd1Index = ai-d;
           if(ai-d > length(polar1.alpha)) %If this alpha is outside of range successfully simulated
               continue; %Skip to next iteration where we try lower alpha
           end
           alphaP1 = polar1.alpha(ai-d);
           if alphaP1 <= polar.alpha(ai) %if this alpha is correct
               break; %We have our index, so end the for loop
           end
       end
       for d2=0:ai-1 %For each alpha value from this to the ones before it
           pd2Index = ai-d2;
           if(ai-d2 > length(polar2.alpha)) %If this alpha is outside of range successfully simulated
               continue; %Skip to next iteration where we try lower alpha
           end
           alphaP2 = polar2.alpha(ai-d2);
           if alphaP2 <= polar.alpha(ai) %if this alpha is correct
               break; %We have our index, so end the for loop
           end
       end
       polar.CL(ai) = polar1.CL(pd1Index) + blendAmt * (polar2.CL(pd2Index) - polar1.CL(pd1Index)); %Linearly interpolate
       polar.CD(ai) = polar1.CD(pd1Index) + blendAmt * (polar2.CD(pd2Index) - polar1.CD(pd1Index)); %Linearly interpolate
    end
    %Find stall angle
    [ClMax, maxClIndex] = max(polar.CL);
    alphaStallDeg = polar.alpha(maxClIndex);
    polar.alphaStallRad = deg2rad(alphaStallDeg);
    polar.Re = polar1.Re + blendAmt * (polar2.Re - polar1.Re); %Linearly interpolate
    polar.Ncrit = polar1.Ncrit + blendAmt * (polar2.Ncrit - polar1.Ncrit); %Linearly interpolate
end

function foilDataMap = getInterpPolars(Radius,numFoilPts,foils,alphas,ReValues,M)
    if isempty(gcp('nocreate'))
        parpoolMsg = msgbox('Initializing parallel pool... (To calculate all airf polars)');
        parpool();
        close(parpoolMsg);
    end
    for ri=1:length(ReValues)
       foilDataMap.Re(ri) = ReValues(ri);
    end
    foilDatas = {length(foilDataMap.Re)}; %Required format for parfor to work
    rootDir = pwd;
    spmd %Put each worker in it's own directory
        cd(rootDir);
        workerDir = [rootDir '\' sprintf('worker%d', labindex)];
        if(~exist(workerDir,'dir'))
            mkdir(workerDir);
        end
        copyfile('xfoil.exe',[workerDir '\xfoil.exe']);
        cd(workerDir);
    end
    parfor ri=1:length(ReValues)
        %Wrap as cell so parfor works
        foilDatas(ri) = {getInterpPolarsForOneRe(Radius,numFoilPts,foils,alphas,ReValues(ri),M)};
    end
    spmd
        cd ..
    end
    cd(rootDir);
    for ri=1:length(ReValues)
       foilDataMap.foilData(ri) = cell2mat(foilDatas(ri)); 
    end
end

%Get a structure containing foil polars from blended foils at different
%points along the radius
function foilData = getInterpPolarsForOneRe(Radius,numFoilPts,foils,alphas,Re,M)
    foilData.alphasUsed = alphas;
    foilData.ReUsed = Re;
    foilData.MUsed = M;
    foilRadiiStep = Radius / (numFoilPts-1); %Increment of radius to calculate 
    for foilPt=1:numFoilPts
       testRad = foilRadiiStep*(foilPt-1);
       foilData.testRad(foilPt) = testRad;
       for j=1:length(foils.startPoint) %For each foil section
          if length(foils.startPoint) > 1 && foils.startPoint(j)*Radius > testRad %If there is more than one section and the section starts further out than us 
              j = j-1; %Set the j (index of the airfoil section to use) to be the prev one in the list
              break; %End this inner for loop
          end
       end
       foilSection1 = char(foils.file(j));
       foilData.desc(foilPt) = "Re"+Re+" M"+M;
       if j < length(foils.startPoint) %If there is another foil section after this one
           foilSection2 = char(foils.file(j+1));
           f1StartRad = foils.startPoint(j)*Radius; %Radius where foil1 starts
           f2StartRad = foils.startPoint(j+1)*Radius; %Radius where foil2 starts
           %Linearly interpolate between airfoils
           blendProportion = (testRad - f1StartRad) / (f2StartRad-f1StartRad);
           fs1Str = convertCharsToStrings(foilSection1);
           fs2Str = convertCharsToStrings(foilSection2);
           if fs1Str == fs2Str
               blendProportion = 0;
           end
           foilData.desc(foilPt) = foilData.desc(foilPt)+ foilSection1+" blend "+foilSection2+" "+blendProportion;
           disp(foilData.desc(foilPt));
       else %This is the only foil section for this region
           foilData.desc(foilPt) = foilData.desc(foilPt)+ foilSection1;
           disp(foilData.desc(foilPt));
       end
       
       %If already calculated this exact case, skip it
       if foilPt > 1 && foilData.desc(foilPt-1) == foilData.desc(foilPt)
           foilData.polar(foilPt) = foilData.polar(foilPt-1);
       else %Need to calculate this case
           cacheFileName = char("..\foilpolars\"+strrep(strrep(foilData.desc(foilPt), '/', ''),'.','')+".mat");
           if exist(cacheFileName, 'file') == 2
              %Load foil polar from cache to save computation time 
              S = load(cacheFileName,'varToSave');
              foilData.polar(foilPt) = S.varToSave;
           else
               writeToCache = true;
               try
               if j < length(foils.startPoint) %If there is another foil section after this one
                    foilData.polar(foilPt) = getBlendedPolar(foilSection1, foilSection2, foilData.desc(foilPt), blendProportion, alphas,Re,M);
               else
                    foilData.polar(foilPt) = xfoil(foilSection1,foilData.desc(foilPt),alphas,Re,M);
               end
               catch err
                   warning('Failed to generate polar for airfoil pt: '+foilData.desc(foilPt)); 
                   writeToCache = false;
                   if(foilPt > 1)
                        foilData.polar(foilPt) = foilData.polar(foilPt-1);
                   end
                   %Try and set a close enough polar for this point
               end
               if (writeToCache)
                   %Save foil data polar to cache
                   varToSave = foilData.polar(foilPt);
                   save(cacheFileName,'varToSave','-mat');
               end
           end
       end
    end
end

%Get a structure containing foil coords from blended foils at different
%points along the radius
function foilData = getInterpFoilCoords(Radius,numFoilPts,foils)
    foilRadiiStep = Radius / (numFoilPts-1); %Increment of radius to calculate 
    for foilPt=1:numFoilPts
       testRad = foilRadiiStep*(foilPt-1);
       foilData.testRad(foilPt) = testRad;
       for j=1:length(foils.startPoint) %For each foil section
          if length(foils.startPoint) > 1 && foils.startPoint(j)*Radius > testRad %If there is more than one section and the section starts further out than us 
              j = j-1; %Set the j (index of the airfoil section to use) to be the prev one in the list
              break; %End this inner for loop
          end
       end
       foilSection1 = char(foils.file(j));
       foilData.desc(foilPt) = " ";
       if j < length(foils.startPoint) %If there is another foil section after this one
           foilSection2 = char(foils.file(j+1));
           f1StartRad = foils.startPoint(j)*Radius; %Radius where foil1 starts
           f2StartRad = foils.startPoint(j+1)*Radius; %Radius where foil2 starts
           %Linearly interpolate between airfoils
           blendProportion = (testRad - f1StartRad) / (f2StartRad-f1StartRad);
           foilData.desc(foilPt) = foilData.desc(foilPt)+ foilSection1+" blend "+foilSection2+" "+blendProportion;
           disp(foilData.desc(foilPt));
       else %This is the only foil section for this region
           foilData.desc(foilPt) = foilData.desc(foilPt)+ foilSection1;
           disp(foilData.desc(foilPt));
       end
       
       %If already calculated this exact case, skip it
       if foilPt > 1 && foilData.desc(foilPt-1) == foilData.desc(foilPt)
           foilData.foil(foilPt) = foilData.foil(foilPt-1);
       else %Need to calculate this case
               if j < length(foils.startPoint) %If there is another foil section after this one
                    foilData.foil(foilPt) = getBlendedFoil(foilSection1, foilSection2, foilData.desc(foilPt), blendProportion, 0,45000,0.03); %Arbitrary alpha, Re, M as not used
               else
                    foilData.foil(foilPt) = getBlendedFoil(foilSection1,foilSection1,foilData.desc(foilPt),0,0,45000,0.03); %Arbitrary alpha, Re, M as not used
               end
       end
    end
end

function [foil] = getBlendedFoil(coord1, coord2, debugName, blendAmt,alpha,Re,Mach,varargin)
    %Execute commands in x foil to blend airfoils
    args(1) = "INTE";
    args(2) = "F";
    args(3) = coord1;
    args(4) = "F";
    args(5) = coord2;
    args(6) = blendAmt;
    args(7) = "BlendedAirfoil";
    args(8) = "PANE";
    for ii = 1:length(varargin)
        args(8+ii) = varargin(ii);
    end
    [~,foil] = xfoil(coord1,debugName,alpha,Re,Mach,args);
end

function [pol] = getBlendedPolar(coord1, coord2, debugName, blendAmt,alpha,Re,Mach,varargin)
    %Execute commands in x foil to blend airfoils
    args(1) = "INTE";
    args(2) = "F";
    args(3) = ['..\',coord1];
    args(4) = "F";
    args(5) = ['..\',coord2];
    args(6) = blendAmt;
    args(7) = "BlendedAirfoil";
    args(8) = "PANE";
    for ii = 1:length(varargin)
        args(8+ii) = varargin(ii);
    end
    [pol] = xfoil(coord1,debugName,alpha,Re,Mach,args);
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