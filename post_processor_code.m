clear all
close all

%Last updated 10/24/25 to include measuring L-junctions from 0-180 degrees
%and X- and T-junctions from 0-90 degrees. This code post processes the
%results of running "junction_formation_code.m" 5 independent times and
%saving the runs as "run1_...", "run2_...", etc.

%Code written by Brittany Bannish

%1st independent run
load run1_junction_anglefull.dat
juncangleA=run1_junction_anglefull;

load run1_components.mat
componentsA=components;

load run1_junction_type.dat
junctypeA=run1_junction_type;

load run1_Ljunc.dat
load run1_parameters.mat

load run1_junction_time.dat
junctimeA=run1_junction_time;

load run1_junction_points.dat
juncptsA=run1_junction_points;


%2nd
load run2_junction_anglefull.dat
juncangleB=run2_junction_anglefull;

load run2_components.mat
componentsB=components;

load run2_junction_type.dat
junctypeB=run2_junction_type;

load run2_Ljunc.dat
load run2_parameters.mat

load run2_junction_time.dat
junctimeB=run2_junction_time;

load run2_junction_points.dat
juncptsB=run2_junction_points;


%3rd
load run3_junction_anglefull.dat
juncangleC=run3_junction_anglefull;

load run3_components.mat
componentsC=components;

load run3_junction_type.dat
junctypeC=run3_junction_type;

load run3_Ljunc.dat
load run3_parameters.mat

load run3_junction_time.dat
junctimeC=run3_junction_time;

load run3_junction_points.dat
juncptsC=run3_junction_points;


%4th
load run4_junction_anglefull.dat
juncangleD=run4_junction_anglefull;

load run4_components.mat
componentsD=components;

load run4_junction_type.dat
junctypeD=run4_junction_type;

load run4_Ljunc.dat
load run4_parameters.mat

load run4_junction_time.dat
junctimeD=run4_junction_time;

load run4_junction_points.dat
juncptsD=run4_junction_points;

%5th
load run5_junction_anglefull.dat
juncangleE=run5_junction_anglefull;

load run5_components.mat
componentsE=components;

load run5_junction_type.dat
junctypeE=run5_junction_type;

load run5_Ljunc.dat
load run5_parameters.mat

load run5_junction_time.dat
junctimeE=run5_junction_time;

load run5_junction_points.dat
juncptsE=run5_junction_points;


LangleA=0;
LangleB=0;
LangleC=0;
LangleD=0;
LangleE=0;

countLangA=0;
countTangA=0;
countXangA=0;

countLangB=0;
countTangB=0;
countXangB=0;

countLangC=0;
countTangC=0;
countXangC=0;

countLangD=0;
countTangD=0;
countXangD=0;

countLangE=0;
countTangE=0;
countXangE=0;

for i=1:length(junctypeA)
    if(junctypeA(i)==1) %L-junctions
        countLangA=countLangA+1;
        LangleA(countLangA,1)=juncangleA(i);
    elseif(junctypeA(i)==2) %T-junctions
        countTangA=countTangA+1;
        TanglefullA(countTangA)=juncangleA(i);
    elseif(junctypeA(i)==3) %X-junctions
        countXangA=countXangA+1;
        XanglefullA(countXangA)=juncangleA(i);
    end
end


for i=1:length(junctypeB)
    if(junctypeB(i)==1) %L-junctions
        countLangB=countLangB+1;
        LangleB(countLangB,1)=juncangleB(i);
    elseif(junctypeB(i)==2) %T-junctions
        countTangB=countTangB+1;
        TanglefullB(countTangB)=juncangleB(i);
    elseif(junctypeB(i)==3) %X-junctions
        countXangB=countXangB+1;
        XanglefullB(countXangB)=juncangleB(i);
    end
end

for i=1:length(junctypeC)
    if(junctypeC(i)==1) %L-junctions
        countLangC=countLangC+1;
        LangleC(countLangC,1)=juncangleC(i);
    elseif(junctypeC(i)==2) %T-junctions
        countTangC=countTangC+1;
        TanglefullC(countTangC)=juncangleC(i);
    elseif(junctypeC(i)==3) %X-junctions
        countXangC=countXangC+1;
        XanglefullC(countXangC)=juncangleC(i);
    end
end

for i=1:length(junctypeD)
    if(junctypeD(i)==1) %L-junctions
        countLangD=countLangD+1;
        LangleD(countLangD,1)=juncangleD(i);
    elseif(junctypeD(i)==2) %T-junctions
        countTangD=countTangD+1;
        TanglefullD(countTangD)=juncangleD(i);
    elseif(junctypeD(i)==3) %X-junctions
        countXangD=countXangD+1;
        XanglefullD(countXangD)=juncangleD(i);
    end
end

for i=1:length(junctypeE)
    if(junctypeE(i)==1) %L-junctions
        countLangE=countLangE+1;
        LangleE(countLangE,1)=juncangleE(i);
    elseif(junctypeE(i)==2) %T-junctions
        countTangE=countTangE+1;
        TanglefullE(countTangE)=juncangleE(i);
    elseif(junctypeE(i)==3) %X-junctions
        countXangE=countXangE+1;
        XanglefullE(countXangE)=juncangleE(i);
    end
end

TangleA=min(TanglefullA,180-TanglefullA)';
XangleA=min(XanglefullA,180-XanglefullA)';

TangleB=min(TanglefullB,180-TanglefullB)';
XangleB=min(XanglefullB,180-XanglefullB)';

TangleC=min(TanglefullC,180-TanglefullC)';
XangleC=min(XanglefullC,180-XanglefullC)';

TangleD=min(TanglefullD,180-TanglefullD)';
XangleD=min(XanglefullD,180-XanglefullD)';

TangleE=min(TanglefullE,180-TanglefullE)';
XangleE=min(XanglefullE,180-XanglefullE)';

filename='model_data.xlsx'

headers = {'X_angle', 'T_angle', 'L_angle'};

% Run A - Sheet 1
n = max([length(XangleA), length(TangleA), length(LangleA)]);
XangleA(end+1:n) = NaN; TangleA(end+1:n) = NaN; LangleA(end+1:n) = NaN;
M = [XangleA(:), TangleA(:), LangleA(:)];
writecell(headers, filename, 'Sheet', 'Run1', 'Range', 'A1');
writematrix(M, filename, 'Sheet', 'Run1', 'Range', 'A2');

% Run B - Sheet 2
n = max([length(XangleB), length(TangleB), length(LangleB)]);
XangleB(end+1:n) = NaN; TangleB(end+1:n) = NaN; LangleB(end+1:n) = NaN;
M = [XangleB(:), TangleB(:), LangleB(:)];
writecell(headers, filename, 'Sheet', 'Run2', 'Range', 'A1');
writematrix(M, filename, 'Sheet', 'Run2', 'Range', 'A2');

% Run C - Sheet 3
n = max([length(XangleC), length(TangleC), length(LangleC)]);
XangleC(end+1:n) = NaN; TangleC(end+1:n) = NaN; LangleC(end+1:n) = NaN;
M = [XangleC(:), TangleC(:), LangleC(:)];
writecell(headers, filename, 'Sheet', 'Run3', 'Range', 'A1');
writematrix(M, filename, 'Sheet', 'Run3', 'Range', 'A2');

% Run D - Sheet 4
n = max([length(XangleD), length(TangleD), length(LangleD)]);
XangleD(end+1:n) = NaN; TangleD(end+1:n) = NaN; LangleD(end+1:n) = NaN;
M = [XangleD(:), TangleD(:), LangleD(:)];
writecell(headers, filename, 'Sheet', 'Run4', 'Range', 'A1');
writematrix(M, filename, 'Sheet', 'Run4', 'Range', 'A2');

% Run E - Sheet 5
n = max([length(XangleE), length(TangleE), length(LangleE)]);
XangleE(end+1:n) = NaN; TangleE(end+1:n) = NaN; LangleE(end+1:n) = NaN;
M = [XangleE(:), TangleE(:), LangleE(:)];
writecell(headers, filename, 'Sheet', 'Run5', 'Range', 'A1');
writematrix(M, filename, 'Sheet', 'Run5', 'Range', 'A2');