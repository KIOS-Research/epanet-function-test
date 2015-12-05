function [ output_args ] = compareEPANET(oldversion, newversion)
%CompareEPANET Compares two EPANET DLLs
%   Example: compareEPANET(', newversion)
    if isdeployed % Stand-alone mode.
        [status, result] = system('path');
        currentDir = char(regexpi(result, 'Path=(.*?);', 'tokens', 'once'));
    else % MATLAB mode.
        addpath(genpath(pwd))
        currentDir = pwd;
    end
    if ~exist('oldversion')
        oldversion='epanet2';
        newversion='epanet2dev21';
    end
	inpname={'Net1.inp','Net2.inp','Net3.inp'};
    fileID = fopen('report.txt','w');
    fprintf(fileID,['Comparison ', oldversion,' with ',newversion,'\n']);
    for i=1:length(inpname)
        unitTesting(oldversion,newversion,inpname{i},i,fileID);
        pause(2)
    end
    fclose(fileID);
    disp('finish')
end

function unitTesting(oldversion,newversion,inpname,looptime,fileID)
    fprintf(fileID,inpname);
    fprintf(fileID,'\n----------------\n');
    epanet1=epanet(inpname,oldversion);
    epanet1.saveInputFile('tmp1.inp');
    
    % simulate all
    HTSepanet1=epanet1.getComputedHydraulicTimeSeries; % Also are included: obj.openHydraulicAnalysis;obj.initializeHydraulicAnalysis;obj.runHydraulicAnalysis;obj.nextHydraulicAnalysisStep;obj.closeHydraulicAnalysis;
    QTSepanet1=epanet1.getComputedQualityTimeSeries;   % Also are included: obj.openQualityAnalysis;obj.initializeQualityAnalysis;obj.runQualityAnalysis;obj.stepQualityAnalysisTimeLeft;obj.closeQualityAnalysis;
    i=1;% Errors
    E=[0:6,101:106,110,120,200,202:207,223:224, 240:241, 250:251,253, 301:309];
    for e=E
        epanet1Error{i}=epanet1.getError(e);i=i+1;
    end
    libfunctionsepanet1 = libfunctions(oldversion);
    %epanet1.closeNetwork
    epanet1.unload;


    epanet2=epanet(inpname,newversion);
    % simulate all
    HTSepanet2=epanet2.getComputedHydraulicTimeSeries;
    QTSepanet2=epanet2.getComputedQualityTimeSeries; 
    
    % Save Inp File
    epanet2.saveInputFile('tmp2.inp');            
    [tlineepanet2]=regexp(fileread('tmp2.inp'), '\n', 'split');
    [tlineepanet1]=regexp(fileread('tmp1.inp'), '\n', 'split');
    if ~isempty(find(strcmp(tlineepanet1,tlineepanet2)==0))
        fprintf(fileID,'Difference in ENsaveinpfile.\n'); 
    end

    if looptime==1
        libfunctionsepanet2 = libfunctions(newversion);
        for i=1:length(libfunctionsepanet2)
            if ~strcmp(libfunctionsepanet2{i},libfunctionsepanet1)
               str = ['Function "',libfunctionsepanet2{i}, '" is a new function.'];
               fprintf(fileID,str)
            end
        end
    end
    %% Errors
    e=E;
    for i=1:length(epanet1Error) 
        if ~strcmp(epanet1Error{i},epanet2.getError(e(i)))
            fprintf(fileID,['Difference in ENgeterror.\n', num2str(e(i))]);
        end
    end
    
    %% Simulate all
    % Hydraulic analysis
    if ~isequal(HTSepanet1.Demand,HTSepanet2.Demand)
        fprintf(fileID,'Difference in demands.\n');
    end
    if ~isequal(HTSepanet1.Energy,HTSepanet2.Energy) 
        fprintf(fileID,'Difference in energy.\n');
    end
    if ~isequal(HTSepanet1.Flow,HTSepanet2.Flow)
        fprintf(fileID,'Difference in flows.\n');
    end
    if ~isequal(HTSepanet1.Head,HTSepanet2.Head)
        fprintf(fileID,'Difference in heads.\n');
    end
    if ~isequal(HTSepanet1.HeadLoss,HTSepanet2.HeadLoss)
        fprintf(fileID,'Difference in headloss.\n');
    end
    if ~isequal(HTSepanet1.Pressure,HTSepanet2.Pressure)
        fprintf(fileID,'Difference in pressures.\n');
    end
    if ~isequal(HTSepanet1.Setting,HTSepanet2.Setting)
        fprintf(fileID,'Difference in settings.\n');
    end
    if ~isequal(HTSepanet1.Status,HTSepanet2.Status) 
        fprintf(fileID,'Difference in status.\n');
    end
    if ~isequal(HTSepanet1.TankVolume,HTSepanet2.TankVolume)
        fprintf(fileID,'Difference in tank volumes.\n'); 
    end
    if ~isequal(HTSepanet1.Time,HTSepanet2.Time) 
        fprintf(fileID,'Difference in times.\n'); 
    end
    if ~isequal(HTSepanet1.Velocity,HTSepanet2.Velocity)
        fprintf(fileID,'Difference in velocity.\n');
    end    
    % Quality analysis
    if ~isequal(QTSepanet1.MassFlowRate,QTSepanet2.MassFlowRate) && ~isequal(isnan(QTSepanet1.MassFlowRate),isnan(QTSepanet2.MassFlowRate))
        fprintf(fileID,'Difference in mass flow rate.\n');
    end
    if ~isequal(QTSepanet1.Quality,QTSepanet2.Quality) 
        fprintf(fileID,'Difference in quality results.\n');
    end
    if ~isequal(QTSepanet1.Time,QTSepanet2.Time)
        fprintf(fileID,'Difference in time.\n');
    end
    
    if ~isequal(epanet1.Version,epanet2.Version)
        fprintf(fileID,'Difference in version.\n');
    end   

    % Controls
    if epanet1.ControlRulesCount
        for i=1:epanet1.ControlRulesCount
            if ~isequal(epanet1.Controls{i},epanet2.Controls{i})
                fprintf(fileID,'Difference in the controls.\n');
            end
        end
    end

    % Counts
    if epanet1.NodeCount~=epanet2.NodeCount
        fprintf(fileID,'Difference in node count.\n');
    end
    if epanet1.NodeTankReservoirCount~=epanet2.NodeTankReservoirCount
        fprintf(fileID,'Difference in tank reservoir count.\n');
    end
    if epanet1.LinkCount~=epanet2.LinkCount
        fprintf(fileID,'Difference in link count.\n');
    end
    if epanet1.PatternCount~=epanet2.PatternCount
        fprintf(fileID,'Difference in pattern count.\n');
    end
    if epanet1.CurveCount~=epanet2.CurveCount
        fprintf(fileID,'Difference in curve count.\n');
    end
    if epanet1.NodeTankCount~=epanet2.NodeTankCount
        fprintf(fileID,'Difference in tank count.\n');
    end
    if epanet1.NodeReservoirCount~=epanet2.NodeReservoirCount
        fprintf(fileID,'Difference in reservoir count.\n');
    end
    if epanet1.NodeJunctionCount~=epanet2.NodeJunctionCount
        fprintf(fileID,'Difference in junctions count.\n');
    end
    if epanet1.LinkPipeCount~=epanet2.LinkPipeCount
        fprintf(fileID,'Difference in pipe count.\n');
    end
    if epanet1.LinkPumpCount~=epanet2.LinkPumpCount
        fprintf(fileID,'Difference in pump count.\n');
    end
    if epanet1.LinkValveCount~=epanet2.LinkValveCount
        fprintf(fileID,'Difference in valve count.\n');
    end

    % Link Flow Units
    if ~strcmp(epanet2.LinkFlowUnits,epanet1.LinkFlowUnits)
        fprintf(fileID,'Difference in flow units.\n');
    end

    % Links Info
    if ~strcmp(epanet1.LinkNameID,epanet2.LinkNameID)
        fprintf(fileID,'Difference in link name ID.\n');
    end
    if ~strcmp(epanet1.LinkPipeNameID,epanet2.LinkPipeNameID)
        fprintf(fileID,'Difference in link pipe name ID.\n');
    end
    if ~strcmp(epanet1.LinkPumpNameID,epanet2.LinkPumpNameID)
        fprintf(fileID,'Difference in link pump name ID.\n');
    end
    if ~strcmp(epanet1.LinkValveNameID,epanet2.LinkValveNameID)
        fprintf(fileID,'Difference in link valve name ID.\n');
    end
    if ~strcmp(epanet1.LinkType,epanet2.LinkType)
        fprintf(fileID,'Difference in link type.\n');
    end
    if ~isequal(epanet1.LinkDiameter,epanet2.LinkDiameter)
        fprintf(fileID,'Difference in link diameter.\n');
    end
    if ~isequal(epanet1.LinkLength,epanet2.LinkLength)
        fprintf(fileID,'Difference in link length.\n');
    end
    if ~isequal(epanet1.LinkRoughnessCoeff,epanet2.LinkRoughnessCoeff)
        fprintf(fileID,'Difference in link roughness.\n');
    end
    if ~isequal(epanet1.LinkInitialStatus, epanet2.LinkInitialStatus)
        fprintf(fileID,'Difference in link status.\n');
    end
    if ~isequal(epanet1.LinkBulkReactionCoeff,epanet2.LinkBulkReactionCoeff)
        fprintf(fileID,'Difference in link bulk reaction coeff.\n');
    end
    if ~isequal(epanet1.LinkWallReactionCoeff,epanet2.LinkWallReactionCoeff)
        fprintf(fileID,'Difference in link wall reaction coeff.\n');
    end


    % Nodes Info
    if ~strcmp(epanet1.NodeNameID,epanet2.NodeNameID)
        fprintf(fileID,'Difference in node name ID.\n');
    end
    if ~strcmp(epanet1.NodeReservoirNameID,epanet2.NodeReservoirNameID)
        fprintf(fileID,'Difference in node reservoir ID.\n');
    end
    if ~strcmp(epanet1.NodeJunctionNameID,epanet2.NodeJunctionNameID)
        fprintf(fileID,'Difference in node junction ID.\n');
    end
    if ~strcmp(epanet1.NodeType,epanet2.NodeType)
        fprintf(fileID,'Difference in node type.\n');
    end
    if ~isequal(epanet1.NodeElevations,epanet2.NodeElevations)
        fprintf(fileID,'Difference in node elevation.\n');
    end
    if ~isequal(epanet1.NodeBaseDemands,epanet2.NodeBaseDemands) % epanet20013
        fprintf(fileID,'Difference in node basedemand (demand categories - epanet2).\n');
    end
    %pp = epanet2.getBinNodesInfo;
    %if ~strcmp(pp.BinNodeDemandPatternNameID,epanet2.NodeDemandPatternNameID) %epanet20013
    if ~strcmp(epanet1.NodeDemandPatternNameID,epanet2.NodeDemandPatternNameID) %epanet20013
        fprintf(fileID,'Difference in pattern name ID.\n');
    end
    if ~isequal(epanet1.NodeInitialQuality,epanet2.NodeInitialQuality)
        fprintf(fileID,'Difference in node initial quality.\n');
    end
    if ~isequal(isnan(epanet1.NodeSourceQuality),isnan(epanet2.NodeSourceQuality))
        fprintf(fileID,'Difference in node source quality.\n');
    end
    if ~isequal(isnan(epanet1.NodeSourcePatternIndex),isnan(epanet2.NodeSourcePatternIndex))
        fprintf(fileID,'Difference in node source pattern.\n');
    end
    if ~isequal(epanet1.NodeSourceTypeIndex,epanet2.NodeSourceTypeIndex)
        fprintf(fileID,'Difference in source type.\n');
    end
    if ~(strcmp(epanet1.NodesConnectingLinksID,epanet2.NodesConnectingLinksID))
        fprintf(fileID,'Difference in node connecting links.\n');
    end
        
    
    % Node Tank data
    if epanet1.NodeTankCount
        if ~isequal(epanet1.NodeTankInitialLevel,epanet2.NodeTankInitialLevel)
            fprintf(fileID,'Difference in tank initial level.\n');
        end
        if ~isequal(epanet1.NodeTankMinimumWaterVolume,epanet2.NodeTankMinimumWaterVolume)
            fprintf(fileID,'Difference in tank minimum water volume.\n');
        end
        if ~isequal(epanet1.NodeTankMixZoneVolume,epanet2.NodeTankMixZoneVolume)
            fprintf(fileID,'Difference in tank minimum water volume.\n');
        end        
        if ~strcmp(epanet1.NodeTankMixingModelType,epanet2.NodeTankMixingModelType)
            fprintf(fileID,'Difference in tank mixing model.\n');
        end
        if ~isequal(epanet1.NodeTankMinimumWaterLevel,epanet2.NodeTankMinimumWaterLevel)
            fprintf(fileID,'Difference in tank minimum water level.\n');
        end
        if ~isequal(epanet1.NodeTankMaximumWaterLevel,epanet2.NodeTankMaximumWaterLevel)
            fprintf(fileID,'Difference in tank maximum water level.\n');
        end
        if ~isequal(epanet1.NodeTankMinimumFraction,epanet2.NodeTankMinimumFraction)
            fprintf(fileID,'Difference in tank minimum fraction.\n');
        end
        if ~strcmp(epanet1.NodeTankNameID,epanet2.NodeTankNameID)
            fprintf(fileID,'Difference in node tank name Iepanet2.\n');
        end
        if ~isequal(epanet1.NodeTankDiameter,epanet2.NodeTankDiameter);% bug: Produces a different diameter % fix in epanet20013 - ok 6/3/2015
            fprintf(fileID,'Difference in tank diameter (bug epanet1).\n');
        end
    end


    % Options 
    if ~isequal(epanet1.OptionsMaxTrials,epanet2.OptionsMaxTrials)
        fprintf(fileID,'Difference in options max trials.\n');
    end
    if ~isequal(epanet1.OptionsAccuracyValue,epanet2.OptionsAccuracyValue)
        fprintf(fileID,'Difference in options accurancy value.\n');
    end
    if ~isequal(epanet1.OptionsQualityTolerance,epanet2.OptionsQualityTolerance)
        fprintf(fileID,'Difference in options quality tolerance.\n');
    end
    if ~isequal(epanet1.OptionsEmitterExponent,epanet2.OptionsEmitterExponent)
        fprintf(fileID,'Difference in options emitter exponent.\n');
    end
    if ~isequal(epanet1.OptionsPatternDemandMultiplier,epanet2.OptionsPatternDemandMultiplier)
        fprintf(fileID,'Difference in options pattern demand multiplier.\n');
    end


    % Pattern data
    if ~strcmp(epanet1.PatternNameID,epanet2.PatternNameID)
        fprintf(fileID,'Difference in pattern name Iepanet2.\n');
    end
    if ~isequal(epanet1.PatternLengths,epanet2.PatternLengths)
        fprintf(fileID,'Difference in pattern lengths.\n');
    end
    if ~isequal(epanet1.Pattern,epanet2.Pattern)
        fprintf(fileID,'Difference in pattern matrix.\n');
    end

    % Pattern average value
    if epanet2.PatternCount
        m = epanet1.Pattern;
        PatternAveragePatternValue = epanet2.PatternAveragePatternValue;
        try
            if ~isequal(PatternAveragePatternValue(1),mean(single(m(1,:))))
                fprintf(fileID,'Difference in ENgetaveragepatternvalue.\n');
            end
        catch e 
        end
    end
    
    
    % Quality 
    if ~isequal(epanet1.QualityCode,epanet2.QualityCode)
        fprintf(fileID,'Difference in quality code.\n');
    end
    if ~isequal(epanet1.QualityTraceNodeIndex,epanet2.QualityTraceNodeIndex)
        if ~isempty(epanet2.QualityTraceNodeIndex) && epanet1.QualityTraceNodeIndex
            fprintf(fileID,'Difference in quality trace node index.\n');
        end
    end
    B=epanet2.getOptionsQualityTolerance;
    epanet2.setQualityType('none'); %bug - to check DE
    A=epanet2.getOptionsQualityTolerance;
    if ~isequal(B,A)
        fprintf(fileID,'Difference in options quality tolerance.\n');
    end
    
    
    % Times
    if ~isequal(epanet1.TimeSimulationDuration,epanet2.TimeSimulationDuration)
        fprintf(fileID,'Difference in time simulation duration.\n');
    end
    if ~isequal(epanet1.TimeHydraulicStep,epanet2.TimeHydraulicStep)
        fprintf(fileID,'Difference in time hydraulic step.\n');
    end
    if ~isequal(epanet1.TimeQualityStep,epanet2.TimeQualityStep)
        fprintf(fileID,'Difference in time quality step.\n');
    end
    if ~isequal(epanet1.TimePatternStep,epanet2.TimePatternStep)
        fprintf(fileID,'Difference in time pattern step.\n');
    end
    if ~isequal(epanet1.TimePatternStart,epanet2.TimePatternStart)
        fprintf(fileID,'Difference in time pattern start.\n');
    end
    if ~isequal(epanet1.TimeReportingStep,epanet2.TimeReportingStep)
        fprintf(fileID,'Difference in time reporting step.\n');
    end
    if ~isequal(epanet1.TimeReportingStart,epanet2.TimeReportingStart)
        fprintf(fileID,'Difference in time reporting start.\n');
    end
    if ~isequal(epanet1.TimeStatisticsIndex,epanet1.TimeStatisticsIndex)
        fprintf(fileID,'Difference in time statistic index.\n');
    end
    if ~strcmp(epanet1.TimeStatisticsType,epanet2.TimeStatisticsType)
        fprintf(fileID,'Difference in time statistic type.\n');
    end
    if ~isequal(epanet1.TimeRuleControlStep,epanet2.TimeRuleControlStep)
        fprintf(fileID,'Difference in set time rule control step (bug epanet1).\n');
    end
    if ~isequal(epanet1.TimeReportingPeriods,epanet2.TimeReportingPeriods)
        fprintf(fileID,'Difference in set time reporting periods (bug epanet1).\n');
    end
    % Add Pattern
    mult=single([0.8, 1.1, 1.4, 1.1, 0.8, 0.7]);
    nPat='NewPat1';
    epanet2.addPattern(nPat,mult); 
    if ~sum(strcmp(epanet2.getPatternNameID,nPat))
        fprintf(fileID,'Difference in pattern name.\n');
    end 
    nPatInd=epanet2.getPatternIndex(nPat);
    matrixPat=epanet2.getPattern;
    if ~isequal(matrixPat(nPatInd,1:length(mult)),mult)
        fprintf(fileID,'Difference in pattern multipliers.\n');
    end  

    % Controls - set
    if epanet2.ControlRulesCount
        controlRuleIndex=1;
        controlTypeIndex=1; % Constants for control: 'LOWLEVEL','HILEVEL', 'TIMER', 'TIMEOFDAY' 0,1,2,3
        linkIndex=13;
        controlSettingValue=1;
        nodeIndex=11;
        controlLevel=150;
        epanet2.setControl(controlRuleIndex,controlTypeIndex,linkIndex,controlSettingValue,nodeIndex,controlLevel); 
        c=epanet2.getControls;
        if c{1}{2}~=controlTypeIndex && c{1}{3}~=linkIndex && c{1}{4}~=controlSettingValue && c{1}{5}~=nodeIndex && c{1}{6}~=controlLevel   
            fprintf(fileID,'Difference in set control.\n');
        end
    end

    % Get, set link parameters
    nlinkDiamBefore=2*epanet2.LinkDiameter;
    epanet2.setLinkDiameter(nlinkDiamBefore);
    nlinkDiamAfter=epanet2.getLinkDiameter;
    if ~isequal(nlinkDiamBefore,nlinkDiamAfter)
        fprintf(fileID,'Difference in set link diameter.\n');
    end

    nlinkLengthBefore=2*epanet2.LinkLength;
    epanet2.setLinkLength(nlinkLengthBefore);
    nlinkLengthAfter=epanet2.getLinkLength;
    if ~isequal(nlinkLengthBefore,nlinkLengthAfter)
        fprintf(fileID,'Difference in set link length.\n');
    end

    nlinkRoughnessCoeffBefore=2*epanet2.LinkRoughnessCoeff;
    epanet2.setLinkRoughnessCoeff(nlinkRoughnessCoeffBefore);
    nlinkRoughnessCoeffAfter=epanet2.getLinkRoughnessCoeff;
    if ~isequal(nlinkRoughnessCoeffBefore,nlinkRoughnessCoeffAfter)
        fprintf(fileID,'Difference in set link roughness coeff.\n');
    end

    nlinkMinorLossCoeffBefore=single(1.1)+epanet2.LinkMinorLossCoeff;
    epanet2.setLinkMinorLossCoeff(nlinkMinorLossCoeffBefore);
    nlinkMinorLossCoeffAfter=epanet2.getLinkMinorLossCoeff;
    if ~isequal(nlinkMinorLossCoeffBefore(1:epanet2.LinkPipeCount),nlinkMinorLossCoeffAfter(1:epanet2.LinkPipeCount)) % pipe
        fprintf(fileID,'Difference in set link minor loss coeff.\n');
    end

    q=epanet2.LinkInitialStatus;
    nlinkInitialStatusBefore=0*epanet2.LinkInitialStatus; % Cannot set status for a check valve
    epanet2.setLinkInitialStatus(nlinkInitialStatusBefore);
    nlinkInitialStatusAfter=epanet2.getLinkInitialStatus;
    if ~isequal(nlinkInitialStatusBefore,nlinkInitialStatusAfter)
        fprintf(fileID,'Difference in set link initial status.\n');
    end
    epanet2.setLinkInitialStatus(q);

    linkinitSetBefore=epanet2.LinkInitialSetting;
    linkinitSetBefore(end)=108;
    epanet2.setLinkInitialSetting(linkinitSetBefore);
    linkinitSetAfter=epanet2.getLinkInitialSetting;
    if ~isequal(linkinitSetBefore,linkinitSetAfter)
        fprintf(fileID,'Difference in set link initial setting.\n');
    end

    linkinitSetBefore=epanet2.LinkInitialSetting;
    linkinitSetBefore(end)=108;
    epanet2.setLinkInitialSetting(linkinitSetBefore);
    linkinitSetAfter=epanet2.getLinkInitialSetting;
    if ~isequal(linkinitSetBefore,linkinitSetAfter)
        fprintf(fileID,'Difference in set link initial setting.\n');
    end

    linkBulkReactionCoeffBefore=epanet2.LinkBulkReactionCoeff-single(0.055);
    epanet2.setLinkBulkReactionCoeff(linkBulkReactionCoeffBefore);
    linkBulkReactionCoeffAfter=epanet2.getLinkBulkReactionCoeff;
    if ~isequal(linkBulkReactionCoeffBefore(1:epanet2.LinkPipeCount),linkBulkReactionCoeffAfter(1:epanet2.LinkPipeCount)) % pipe
        fprintf(fileID,'Difference in set link bulk reaction coeff.\n');
    end

    linkWallReactionCoeffBefore=epanet2.LinkWallReactionCoeff*-single(1.1);
    epanet2.setLinkWallReactionCoeff(linkWallReactionCoeffBefore);
    linkWallReactionCoeffAfter=epanet2.getLinkWallReactionCoeff;
    if ~isequal(linkWallReactionCoeffBefore(1:epanet2.LinkPipeCount),linkWallReactionCoeffAfter(1:epanet2.LinkPipeCount)) % pipe
        fprintf(fileID,'Difference in set link wall reaction coeff.\n');
    end

    linkStatusBefore=ones(1,length(epanet2.LinkStatus));
    epanet2.setLinkStatus(linkStatusBefore);
    linkStatusAfter=epanet2.getLinkStatus;
    if ~isequal(linkStatusAfter,linkStatusAfter)
        fprintf(fileID,'Difference in set link status.\n');
    end

    linkSettingBefore=epanet2.getLinkSettings;
    linkSettingBefore(end)=111;
    epanet2.setLinkSettings(linkSettingBefore);
    linkSettingAfter=epanet2.getLinkSettings;
    if ~isequal(linkSettingBefore,linkSettingAfter)
        fprintf(fileID,'Difference in set link setting.\n');
    end

    % Get, set node parameters

    nodeElevationsBefore=epanet2.NodeElevations;
    nodeElevationsBefore(end)=nodeElevationsBefore(end)+10;
    epanet2.setNodeElevations(nodeElevationsBefore);
    nodeElevationsAfter=epanet2.getNodeElevations;
    if ~isequal(nodeElevationsBefore,nodeElevationsAfter)
        fprintf(fileID,'Difference in set node elevations.\n');
    end

    % NodeBaseDemandsBefore=epanet2.NodeBaseDemands;
    % NodeBaseDemandsBefore(2)=NodeBaseDemandsBefore(2)+5;
    % epanet2.setNodeBaseDemands(NodeBaseDemandsBefore);
    % NodeBaseDemandsAfter=epanet2.NodeBaseDemands;
    % if ~isequal(NodeBaseDemandsBefore,NodeBaseDemandsAfter)
    %     fprintf(fileID,'Difference in set node basedemands.\n');
    % end
    % 
    % NodeDemandPatternIndexBefore=epanet2.NodeDemandPatternIndex;
    % NodeDemandPatternIndexBefore(2)=0;
    % epanet2.setNodeDemandPatternIndex(NodeDemandPatternIndexBefore);
    % NodeDemandPatternIndexAfter=epanet2.NodeDemandPatternIndex;
    % if ~isequal(NodeDemandPatternIndexBefore,NodeDemandPatternIndexAfter)
    %     fprintf(fileID,'Difference in set node demand pattern index.\n');
    % end

    NodeEmitterCoeffBefore=epanet2.NodeEmitterCoeff;
    NodeEmitterCoeffBefore(2)=0.5;
    epanet2.setNodeEmitterCoeff(NodeEmitterCoeffBefore);
    NodeEmitterCoeffAfter=epanet2.getNodeEmitterCoeff;
    if ~isequal(NodeEmitterCoeffBefore,NodeEmitterCoeffAfter)
        fprintf(fileID,'Difference in set node emitter coeff.\n');
    end

    NodeInitialQualityBefore=epanet2.NodeInitialQuality;
    NodeInitialQualityBefore(2)=single(0.6);
    epanet2.setNodeInitialQuality(NodeInitialQualityBefore);
    NodeInitialQualityAfter=epanet2.getNodeInitialQuality;
    if ~isequal(NodeInitialQualityBefore,NodeInitialQualityAfter)
        fprintf(fileID,'Difference in set node initial quality.\n');
    end

    % Get, set node tank parameters
    if epanet2.NodeTankCount
        NodeTankInitialLevelBefore=epanet2.NodeTankInitialLevel;
        NodeTankInitialLevelBefore(end)=NodeTankInitialLevelBefore(end)+5;
        epanet2.setNodeTankInitialLevel(NodeTankInitialLevelBefore);
        NodeTankInitialLevelAfter=epanet2.getNodeTankInitialLevel;
        if ~isequal(NodeTankInitialLevelBefore,NodeTankInitialLevelAfter)
            fprintf(fileID,'Difference in set node tank initial level.\n');
        end

        NodeTankMixingModelTypeBefore=epanet2.NodeTankMixingModelType;
        NodeTankMixingModelTypeBefore{end}='MIX2';
        epanet2.setNodeTankMixingModelType(NodeTankMixingModelTypeBefore);
        NodeTankMixingModelTypeAfter=epanet2.getNodeTankMixingModelType;
        if ~isequal(NodeTankMixingModelTypeBefore,NodeTankMixingModelTypeAfter)
            fprintf(fileID,'Difference in set node tank mixing model type.\n');
        end
        NodeTankMixingModelTypeBefore{end}='FIFO';
        epanet2.setNodeTankMixingModelType(NodeTankMixingModelTypeBefore)
        NodeTankMixingModelTypeAfter=epanet2.getNodeTankMixingModelType;
        if ~isequal(NodeTankMixingModelTypeBefore,NodeTankMixingModelTypeAfter)
            fprintf(fileID,'Difference in set node tank mixing model type.\n');
        end
        NodeTankMixingModelTypeBefore{end}='LIFO';
        epanet2.setNodeTankMixingModelType(NodeTankMixingModelTypeBefore);
        NodeTankMixingModelTypeAfter=epanet2.getNodeTankMixingModelType;
        if ~isequal(NodeTankMixingModelTypeBefore,NodeTankMixingModelTypeAfter)
            fprintf(fileID,'Difference in set node tank mixing model type.\n');
        end

        NodeTankDiameterBefore=epanet2.NodeTankDiameter;
        NodeTankDiameterBefore(end)=NodeTankDiameterBefore(end)+20;
        epanet2.setNodeTankDiameter(NodeTankDiameterBefore);
        NodeTankDiameterAfter=epanet2.getNodeTankDiameter;
        if ~isequal(NodeTankDiameterBefore,NodeTankDiameterAfter)
            fprintf(fileID,'Difference in set node tank diameter.\n');
        end

        NodeTankMinimumWaterLevelBefore=epanet2.NodeTankMinimumWaterLevel;
        NodeTankMinimumWaterLevelBefore(end)=NodeTankMinimumWaterLevelBefore(end)+5;
        epanet2.setNodeTankMinimumWaterLevel(NodeTankMinimumWaterLevelBefore);
        NodeTankMinimumWaterLevelAfter=epanet2.getNodeTankMinimumWaterLevel;
        if ~isequal(NodeTankMinimumWaterLevelBefore,NodeTankMinimumWaterLevelAfter)
            fprintf(fileID,'Difference in set node tank minimum water level.\n');
        end
        
        NodeTankMinimumWaterVolumeBefore=epanet2.NodeTankMinimumWaterVolume;
        NodeTankMinimumWaterVolumeBefore(end)=NodeTankMinimumWaterVolumeBefore(end)+5;
        epanet2.setNodeTankMinimumWaterVolume(NodeTankMinimumWaterVolumeBefore);
        NodeTankMinimumWaterVolumeAfter=epanet2.getNodeTankMinimumWaterVolume;
        if ~isequal(NodeTankMinimumWaterVolumeBefore,NodeTankMinimumWaterVolumeAfter)
            fprintf(fileID,'Difference in set node tank minimum water volume.\n');
        end    

        NodeTankMaximumWaterLevelBefore=epanet2.NodeTankMaximumWaterLevel;
        NodeTankMaximumWaterLevelBefore(end)=NodeTankMaximumWaterLevelBefore(end)+20;
        epanet2.setNodeTankMaximumWaterLevel(NodeTankMaximumWaterLevelBefore);
        NodeTankMaximumWaterLevelAfter=epanet2.getNodeTankMaximumWaterLevel;
        if ~isequal(NodeTankMaximumWaterLevelBefore,NodeTankMaximumWaterLevelAfter)
            fprintf(fileID,'Difference in set node tank maximum water level.\n');
        end  

        NodeTankMinimumFractionBefore=epanet2.NodeTankMinimumFraction;
        NodeTankMinimumFractionBefore(end)=0.5; %takes values 0-1
        epanet2.setNodeTankMinimumFraction(NodeTankMinimumFractionBefore);
        NodeTankMinimumFractionAfter=epanet2.getNodeTankMinimumFraction;
        if ~isequal(NodeTankMinimumFractionBefore,NodeTankMinimumFractionAfter)
            fprintf(fileID,'Difference in set node tank minimum fraction.\n');
        end  

        NodeTankBulkReactionCoeffBefore=epanet2.NodeTankBulkReactionCoeff;
        NodeTankBulkReactionCoeffBefore(end)=1;
        epanet2.setNodeTankBulkReactionCoeff(NodeTankBulkReactionCoeffBefore);
        NodeTankBulkReactionCoeffAfter=epanet2.getNodeTankBulkReactionCoeff;
        if ~isequal(NodeTankBulkReactionCoeffBefore,NodeTankBulkReactionCoeffAfter)
            fprintf(fileID,'Difference in set node tank bulk reaction coeff.\n');
        end      
    end


    % Get, set node source parameters
    nIndex = [1 2 3 4];
    NodeSourceTypeBefore='MASS';
    epanet2.setNodeSourceType(nIndex(1),NodeSourceTypeBefore);
    NodeSourceTypeAfter=epanet2.getNodeSourceType;
    if ~isequal(NodeSourceTypeBefore,NodeSourceTypeAfter{nIndex(1)})
        fprintf(fileID,'Difference in set node source type - mass.\n');
    end     
    NodeSourceTypeBefore='CONCEN';
    epanet2.setNodeSourceType(nIndex(2),NodeSourceTypeBefore);
    NodeSourceTypeAfter=epanet2.getNodeSourceType;
    if ~isequal(NodeSourceTypeBefore,NodeSourceTypeAfter{nIndex(2)})
        fprintf(fileID,'Difference in set node source type - concen.\n');
    end       
    NodeSourceTypeBefore='SETPOINT';
    epanet2.setNodeSourceType(nIndex(3),NodeSourceTypeBefore);
    NodeSourceTypeAfter=epanet2.getNodeSourceType;
    if ~isequal(NodeSourceTypeBefore,NodeSourceTypeAfter{nIndex(3)})
        fprintf(fileID,'Difference in set node source type - setpoint.\n');
    end   
    NodeSourceTypeBefore='FLOWPACED';
    epanet2.setNodeSourceType(nIndex(4),NodeSourceTypeBefore);
    NodeSourceTypeAfter=epanet2.getNodeSourceType;
    if ~isequal(NodeSourceTypeBefore,NodeSourceTypeAfter{nIndex(4)})
        fprintf(fileID,'Difference in set node source type - flowpaceepanet2.\n');
    end   

    NodeSourceQualityBefore=epanet2.getNodeSourceQuality;
    NodeSourceQualityBefore(2)=0.5;
    epanet2.setNodeSourceQuality(NodeSourceQualityBefore);
    NodeSourceQualityAfter=epanet2.getNodeSourceQuality;
    if ~isequal(NodeSourceQualityBefore(nIndex),NodeSourceQualityAfter(nIndex))
        fprintf(fileID,'Difference in set node source quality.\n');
    end

    NodeSourcePatternIndexBefore=epanet2.getNodeSourcePatternIndex;
    NodeSourcePatternIndexBefore(4)=1;
    epanet2.setNodeSourcePatternIndex(NodeSourcePatternIndexBefore);
    NodeSourcePatternIndexAfter=epanet2.getNodeSourcePatternIndex;
    NodeSourcePatternIndexBefore(isnan(NodeSourcePatternIndexBefore))='';
    NodeSourcePatternIndexAfter(isnan(NodeSourcePatternIndexAfter))='';
    if ~isequal(NodeSourcePatternIndexBefore,NodeSourcePatternIndexAfter)
        fprintf(fileID,'Difference in set node source pattern index.\n');
    end

    
    % Get, set options parameters
    OptionsMaxTrialsBefore=45;
    epanet2.setOptionsMaxTrials(OptionsMaxTrialsBefore);
    OptionsMaxTrialsAfter=epanet2.getOptionsMaxTrials;
    if ~isequal(OptionsMaxTrialsBefore,OptionsMaxTrialsAfter)
        fprintf(fileID,'Difference in set options max trials.\n');
    end

    OptionsAccuracyValueBefore=single(0.015);
    epanet2.setOptionsAccuracyValue(OptionsAccuracyValueBefore);
    OptionsAccuracyValueAfter=epanet2.getOptionsAccuracyValue;
    if ~isequal(OptionsAccuracyValueBefore,OptionsAccuracyValueAfter)
        fprintf(fileID,'Difference in set options accurancy value.\n');
    end

    OptionsQualityToleranceBefore=single(0.02);
    epanet2.setOptionsQualityTolerance(OptionsQualityToleranceBefore);
    OptionsQualityToleranceAfter=epanet2.getOptionsQualityTolerance;
    if ~isequal(OptionsQualityToleranceBefore,OptionsQualityToleranceAfter)
        fprintf(fileID,'Difference in set options quality tolerance.\n');
    end

    OptionsEmitterExponentBefore=single(0.55);
    epanet2.setOptionsEmitterExponent(OptionsEmitterExponentBefore);
    OptionsEmitterExponentAfter=epanet2.getOptionsEmitterExponent;
    if ~isequal(OptionsEmitterExponentBefore,OptionsEmitterExponentAfter)
        fprintf(fileID,'Difference in set options emitter exponent.\n');
    end

    OptionsPatternDemandMultiplierBefore=single(1.1);
    epanet2.setOptionsPatternDemandMultiplier(OptionsPatternDemandMultiplierBefore);
    OptionsPatternDemandMultiplierAfter=epanet2.getOptionsPatternDemandMultiplier;
    if ~isequal(OptionsPatternDemandMultiplierBefore,OptionsPatternDemandMultiplierAfter)
        fprintf(fileID,'Difference in set options pattern demand multiplier.\n');
    end

    
    % Get, set time parameters
    TimeSimulationDurationBefore=single(86500);
    epanet2.setTimeSimulationDuration(TimeSimulationDurationBefore);
    TimeSimulationDurationAfter=epanet2.getTimeSimulationDuration;
    if ~isequal(TimeSimulationDurationBefore,TimeSimulationDurationAfter)
        fprintf(fileID,'Difference in set time simulation duration.\n');
    end

    TimeHydraulicStepBefore=epanet2.TimeHydraulicStep-300;   % Hstep = min(Pstep,Hstep) , Hstep = min(Rstep,Hstep) , Qstep = min(Qstep,Hstep)
    epanet2.setTimeHydraulicStep(TimeHydraulicStepBefore);      % Hstep > Rstep ==> Hstep = Rstep
    TimeHydraulicStepAfter=epanet2.getTimeHydraulicStep;
    if ~isequal(TimeHydraulicStepBefore,TimeHydraulicStepAfter)
        fprintf(fileID,'Difference in set time hydraulic step.\n');
    end

    TimeQualityStepBefore=single(250);
    epanet2.setTimeQualityStep(TimeQualityStepBefore);
    TimeQualityStepAfter=epanet2.getTimeQualityStep;
    if ~isequal(TimeQualityStepBefore,TimeQualityStepAfter)
        fprintf(fileID,'Difference in set time quality step.\n');
    end

    TimePatternStepBefore=single(3000);
    epanet2.setTimePatternStep(TimePatternStepBefore);
    TimePatternStepAfter=epanet2.getTimePatternStep;
    if ~isequal(TimePatternStepBefore,TimePatternStepAfter)
        fprintf(fileID,'Difference in set time pattern step.\n');
    end

    TimePatternStartBefore=single(100);
    epanet2.setTimePatternStart(TimePatternStartBefore);
    TimePatternStartAfter=epanet2.getTimePatternStart;
    if ~isequal(TimePatternStartBefore,TimePatternStartAfter)
        fprintf(fileID,'Difference in set time pattern start.\n');
    end

    TimeReportingStepBefore=single(3500);
    epanet2.setTimeReportingStep(TimeReportingStepBefore);
    TimeReportingStepAfter=epanet2.getTimeReportingStep;
    if ~isequal(TimeReportingStepBefore,TimeReportingStepAfter)
        fprintf(fileID,'Difference in set time reporting step.\n');
    end

    TimeReportingStartBefore=single(200);
    epanet2.setTimeReportingStart(TimeReportingStartBefore);
    TimeReportingStartAfter=epanet2.getTimeReportingStart;
    if ~isequal(TimeReportingStartBefore,TimeReportingStartAfter)
        fprintf(fileID,'Difference in set time reporting start.\n');
    end

    TimeStatisticsTypeBefore={'MINIMUM'};
    epanet2.setTimeStatisticsType(TimeStatisticsTypeBefore);
    TimeStatisticsTypeAfter=epanet2.getTimeStatisticsType;
    if ~isequal(TimeStatisticsTypeBefore,TimeStatisticsTypeAfter)
        fprintf(fileID,'Difference in set time statistic - minimum.\n');
    end   
    TimeStatisticsTypeBefore={'MAXIMUM'};
    epanet2.setTimeStatisticsType(TimeStatisticsTypeBefore);
    TimeStatisticsTypeAfter=epanet2.getTimeStatisticsType;
    if ~isequal(TimeStatisticsTypeBefore,TimeStatisticsTypeAfter)
        fprintf(fileID,'Difference in set time statistic - maximum.\n');
    end 
    TimeStatisticsTypeBefore={'RANGE'};
    epanet2.setTimeStatisticsType(TimeStatisticsTypeBefore);
    TimeStatisticsTypeAfter=epanet2.getTimeStatisticsType;
    if ~isequal(TimeStatisticsTypeBefore,TimeStatisticsTypeAfter)
        fprintf(fileID,'Difference in set time statistic - range.\n');
    end 
    TimeStatisticsTypeBefore={'AVERAGE'};
    epanet2.setTimeStatisticsType(TimeStatisticsTypeBefore);
    TimeStatisticsTypeAfter=epanet2.getTimeStatisticsType;
    if ~isequal(TimeStatisticsTypeBefore,TimeStatisticsTypeAfter)
        fprintf(fileID,'Difference in set time statistic - average.\n');
    end 
    TimeStatisticsTypeBefore={'NONE'};
    epanet2.setTimeStatisticsType(TimeStatisticsTypeBefore);
    TimeStatisticsTypeAfter=epanet2.getTimeStatisticsType;
    if ~isequal(TimeStatisticsTypeBefore,TimeStatisticsTypeAfter)
        fprintf(fileID,'Difference in set time statistic - none.\n');
    end 

    TimeRuleControlStepBefore=single(100);
    epanet2.setTimeRuleControlStep(TimeRuleControlStepBefore); 
    TimeRuleControlStepAfter=epanet2.getTimeRuleControlStep;
    if ~isequal(TimeRuleControlStepBefore,TimeRuleControlStepAfter)
        fprintf(fileID,'Difference in set time rule control step (bug epanet1.\n');
    end

    
    % Get, set patterns
    PatternBefore=single(1:0.01:2);pat=1;
    epanet2.setPattern(pat,PatternBefore); 
    PatternAfter=epanet2.getPattern;
    if ~isequal(PatternBefore,PatternAfter(1,:))
        fprintf(fileID,'Difference in set pattern.\n');
    end

    PatternMatrixBefore=epanet2.Pattern;
    PatternMatrixBefore(1,end)=3;
    epanet2.setPatternMatrix(PatternMatrixBefore); 
    PatternMatrixAfter=epanet2.getPattern;
    if ~isequal(PatternMatrixBefore(1,end),PatternMatrixAfter(1,end))
        fprintf(fileID,'Difference in set pattern matrix.\n');
    end

    PatternValueBefore=single(1.2);pat=1;mult=1;
    epanet2.setPatternValue(pat,mult,PatternValueBefore); 
    PatternValueAfter=epanet2.getPatternValue(pat,mult);
    if ~isequal(PatternValueBefore,PatternValueAfter)
        fprintf(fileID,'Difference in set pattern value.\n');
    end

    % Get, set quality
    QualityTypeBefore='none';
    epanet2.setQualityType(QualityTypeBefore);
    QualityTypeAfter=epanet2.getQualityType;
    if ~strcmp(upper(QualityTypeBefore),QualityTypeAfter)
        fprintf(fileID,'Difference in set quality type - none.\n');
    end 
    QualityTypeBefore='age';
    epanet2.setQualityType(QualityTypeBefore);
    QualityTypeAfter=epanet2.getQualityType;
    if ~strcmp(upper(QualityTypeBefore),QualityTypeAfter)
        fprintf(fileID,'Difference in set quality type - age.\n');
    end 
    QualityTypeBefore='chem';
    epanet2.setQualityType(QualityTypeBefore);
    QualityTypeAfter=epanet2.getQualityType;
    if ~strcmp(upper({QualityTypeBefore}),upper(QualityTypeAfter))
        fprintf(fileID,'Difference in set quality type - chem.\n');
    end 
    tankid=epanet2.NodeTankNameID;
    QualityTypeBefore='trace';
    epanet2.setQualityType(QualityTypeBefore,tankid{1});
    QualityTypeAfter=epanet2.getQualityType;
    if ~strcmp(upper(QualityTypeBefore),QualityTypeAfter)
        fprintf(fileID,'Difference in set quality type - trace.\n');
    end  

    
    % Report Preparation
    % Write line in report file
    % Hydraulic analysis
    epanet2.writeLineInReportFile('Line-writting testing')
    epanet2.solveCompleteHydraulics; % solves internally the hydraulics (does not return something)
    
    [SReportepanet2]=regexp(fileread([epanet2.Bintempfile(1:end-4),'.txt']), '\n', 'split');
    l=0;
    for i=1:length(SReportepanet2)
        if ~isempty(regexp(SReportepanet2{i},'Line-writting testing','match'))
            l=1;
        end
    end
    if l==0
        fprintf(fileID,'Difference in write line.\n');
    end
    %epanet2.saveHydraulicsOutputReportingFile %creates a BIN file (see EPANET documentation)

    % Case 1
    % Compute ranges (max - min) 
    %     epanet2.setTimeStatisticsType('RANGE')
    %     epanet2.setTimeStatisticsType('MINIMUM')
    % StatisticsType('AVERAGE')
    epanet2.setTimeStatisticsType('NONE')
    %     epanet2.setTimeStatisticsType('MAXIMUM')
    % Define contents of the report
    epanet2.setReportFormatReset
    epanet2.setReport('FILE TestReport.txt');
    epanet2.setReport('PAGESIZE 2')
    epanet2.setReport('NODES ALL')%/ALL/node1 node2 
    epanet2.setReport('LINKS ALL')%/ALL/link1 link2
    %epanet2.setReport('PRESSURE PRECISION 1')
    epanet2.setReport('SUMMARY NO')%YES/NO 
    epanet2.setReport('MESSAGES YES')%YES/NO 
    epanet2.setReport('ENERGY YES')%YES/NO 
    %Nodes parameters
    %YES/NO/BELOW/ABOVE/PRECISION
    epanet2.setReport('ELEVATION YES')
    epanet2.setReport('DEMAND YES')
    epanet2.setReport('HEAD YES') 
    epanet2.setReport('PRESSURE YES') 
    epanet2.setReport('QUALITY YES') 
    %Links parameters
    %BELOW/ABOVE/PRECISION
    epanet2.setReport('LENGTH YES')
    epanet2.setReport('DIAMETER YES')
    epanet2.setReport('FLOW YES')
    epanet2.setReport('LENGTH YES')
    epanet2.setReport('VELOCITY YES')
    epanet2.setReport('HEADLOSS YES')
    %epanet2.setReport('QUALITY PRECISION 1')
    epanet2.setReport('STATUS YES')%YES/NO/FULL
    epanet2.setReport('SETTING YES')
    epanet2.setReport('REACTION YES')
    epanet2.setReport('F-FACTOR YES')

    %Write the report to file 
    testrep='TestReport.txt';
    if exist(testrep)
        delete(testrep);
    end
    epanet2.writeReport
    if ~exist(testrep)
        fprintf(fileID,'Difference in write report.\n');
    else
        [SReportepanet2]=regexp(fileread(testrep), '\n', 'split');
        dataCheck = {'Page 2','Energy Usage','Elevation','Demand','Head','Pressure','Quality','Length',...
        'Diameter','Flow','Length','Velocity','Headloss','QUALITY PRECISION 1',...
        'Status','Setting','Reaction','F-Factor'};
        mm=zeros(1,length(dataCheck));
        dataSet = {'PAGESIZE 2','ENERGY YES','ELEVATION YES','DEMAND YES','HEAD YES','PRESSURE YES','QUALITY YES','LENGTH YES',...
        'DIAMETER YES','FLOW YES','LENGTH YES','VELOCITY YES','HEADLOSS YES','QUALITY PRECISION 1',...
        'STATUS YES','SETTING YES','REACTION YES','F-FACTOR YES'};
        for u=1:length(mm)
            for i=1:length(SReportepanet2)
                if ~isempty(regexp(SReportepanet2{i},dataCheck{u},'match')) 
                    mm(u)=1;break;
                end
            end
        end
        for u=1:length(mm)
            if mm(u)==0
                str = ['Difference in set report ',dataSet{u},'\n'];
                fprintf(fileID,str);
            end
        end
    end
        
    % Case 2
    epanet2.setReportFormatReset
    epanet2.setReport('FILE TestReport2.txt'); 
    epanet2.setTimeStatisticsType('AVERAGE')
    epanet2.setReport('NODES 10')
    epanet2.setReport('HEAD YES')
    epanet2.setReport('PRESSURE NO')
    epanet2.writeReport
    testrep='TestReport2.txt';
    if exist(testrep)
        delete(testrep);
    end
    epanet2.writeReport
    if ~exist(testrep)
        fprintf(fileID,'Difference in write report.\n');
    else
        [SReportepanet2]=regexp(fileread(testrep), '\n', 'split');
        dataCheck = {'AVERAGE','Selected Nodes','Head','Pressure'};
        mm=zeros(1,length(dataCheck));
        dataSet = {'AVERAGE','NODES 10','HEAD YES','PRESSURE NO'};
        for u=1:length(mm)
            for i=1:length(SReportepanet2)
                if ~isempty(regexp(SReportepanet2{i},dataCheck{u},'match')) 
                    mm(u)=1;break;
                end
            end
        end
        for u=1:length(mm)
            if mm(u)==0
                if isempty(regexp(dataSet{u},'NO','match'))
                    str = ['Difference in set report ',dataSet{u}];
                    fprintf(fileID,str);
                end
            end
        end
    end
    
    % Case 3
    epanet2.setReportFormatReset
    epanet2.setReport('FILE TestReport3.txt'); 
    epanet2.setReport('NODES ALL')
    epanet2.setReport('LINKS ALL')
    epanet2.setReport('PRESSURE ABOVE 20')
    epanet2.setReport('SUMMARY YES')
    epanet2.writeReport
    testrep='TestReport3.txt';
    if exist(testrep)
        delete(testrep);
    end
    epanet2.writeReport
    if ~exist(testrep)
        fprintf(fileID,'Difference in write report.\n');
    else
        [SReportepanet2]=regexp(fileread(testrep), '\n', 'split');
        dataCheck = {'All Nodes','All Links','Pressure above 20.00','Input Data File'};
        mm=zeros(1,length(dataCheck));
        dataSet = {'NODES ALL','LINKS ALL','PRESSURE ABOVE 20','SUMMARY YES'};
        for u=1:length(mm)
            for i=1:length(SReportepanet2)
                if ~isempty(regexp(SReportepanet2{i},dataCheck{u},'match')) 
                    mm(u)=1;break;
                end
            end
        end
        for u=1:length(mm)
            if mm(u)==0
                if isempty(regexp(dataSet{u},'NO','match'))
                    str = ['Difference in set report ',dataSet{u},'\n'];
                    fprintf(fileID,str);
                end
            end
        end
    end
    
    % Case 4
    epanet2.setReportFormatReset
    epanet2.setReport('FILE TestReport4.txt'); 
    epanet2.setTimeStatisticsType('NONE')
    epanet2.setReport('LINKS 10')
    epanet2.setReport('LINKS 11')
    epanet2.setReport('LINKS 12')
    epanet2.setReport('FLOW YES')
    epanet2.setReport('DEMAND NO')
    epanet2.setReport('QUALITY NO')
    epanet2.setReport('SUMMARY NO')
    epanet2.setReport('HEADLOSS NO') %bug: It shoes Headloss even though it shouldn't
    epanet2.setReport('VELOCITY NO')
    epanet2.writeReport
    testrep='TestReport4.txt';
    if exist(testrep)
        delete(testrep);
    end
    epanet2.writeReport
    if ~exist(testrep)
        fprintf(fileID,'Difference in write report.\n');
    else
        [SReportepanet2]=regexp(fileread(testrep), '\n', 'split');
        dataCheck = {'Demand','Qualiy','Headloss','Velocity'};
        mm=ones(1,length(dataCheck));
        dataSet = {'DEMAND NO','QUALITY NO','HEADLOSS NO','VELOCITY NO'};
        for u=1:length(mm)
            for i=1:length(SReportepanet2)
                if ~isempty(regexp(SReportepanet2{i},dataCheck{u},'match')) 
                    mm(u)=0;break;
                end
            end
        end
        for u=1:length(mm)
            if mm(u)==0
                str = ['Difference in set report ',dataSet{u}];%bug: It shoes Headloss even though it shouldn't
                fprintf(fileID,str);
            end
        end
    end
    
    % Case 5
    epanet2.setReportFormatReset
    epanet2.setReport('FILE TestReport5.txt'); 
    epanet2.setTimeStatisticsType('MINIMUM')
    epanet2.setReport('NODES ALL')
    epanet2.writeReport
    testrep='TestReport5.txt';
    if exist(testrep)
        delete(testrep);
    end
    epanet2.writeReport
    if ~exist(testrep)
        fprintf(fileID,'Difference in write report.\n');
    else
        [SReportepanet2]=regexp(fileread(testrep), '\n', 'split');
        dataCheck = {'All Nodes','MINIMUM'};
        mm=zeros(1,length(dataCheck));
        dataSet = {'NODES ALL','MINIMUM'};
        for u=1:length(mm)
            for i=1:length(SReportepanet2)
                if ~isempty(regexp(SReportepanet2{i},dataCheck{u},'match')) 
                    mm(u)=1;break;
                end
            end
        end
        for u=1:length(mm)
            if mm(u)==0
                str = ['Difference in set report ',dataSet{u}];%bug: It shoes Headloss even though it shouldn't
                fprintf(fileID,str);
            end
        end
    end
    %epanet2.closeNetwork
    epanet2.unload
end