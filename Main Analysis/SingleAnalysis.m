% load the data
% add helpers to the path
%% Analysis of order 1:
participants =51;
conditions =1;
trials = 15;
L=6;
ts_length = 1026;

for cond =1:conditions
    writematrix('Condition: ' + string(cond),'ParticipantLevelAnalysisLevel1.txt','WriteMode','append')
    writematrix('Condition: ' + string(cond),'ConditionLevelAnalysisLevel1.txt','WriteMode','append')
    
    %ORDER 2

    for alpha=1:L
        alpha
        writematrix(string(alpha),'ParticipantLevelAnalysisLevel1.txt','WriteMode','append')
        JointDistInPPT = zeros(1,1);
        JointDistOutPPT = zeros(1,1);
        currentMaxJointPPT = 1;
        for participant =1:participants
            currentMaxSingle =1;
            for trial =1:trials
                %Extract the desired time-series
                X = MEGData{participant}(:,:,trial);
                Y = X(alpha,:);
                %Build the network
                A = MultiplexVisibilityNetwork(Y);
                X = SampledData{participant}.Condition{cond}(:,:,trial);
                %Build the network
                A = MultiplexVisibilityNetwork(X);
    
                % Build a 1-d,2-d,3-d distribution for each trial
                [K_In,K_Out] = SignedDegrees(A);
                DistIn = DegreeDist(K_In);
                DistOut = DegreeDist(K_Out);
                maxIn = max(size(DistIn,1)); maxOut = max(size(DistOut,1));
                maxDegree = max(maxIn,maxOut);
                DistInTemp = DistIn;
                DistIn = zeros(maxDegree,L);
                DistOutTemp = DistOut;
                DistOut = zeros(maxDegree,L);
                DistIn(1:maxIn,:) = DistInTemp;
                DistOut(1:maxOut,:) = DistOutTemp;

                if trial ==1
                    DistInTrials = zeros(trials,maxDegree,L);
                    DistOutTrials = zeros(trials,maxDegree,L);
                    DistInTrials(trial,:,:) = DistIn;
                    DistOutTrials(trial,:,:) = DistOut;
                    currentMaxSingle = maxDegree;
                else
                    if maxDegree>currentMaxSingle
                        DistInTrialsTemp = DistInTrials;
                        DistOutTrialsTemp = DistOutTrials;
                        DistInTrials = zeros(trials,maxDegree,L);
                        DistOutTrials = zeros(trials,maxDegree,L);
                        DistInTrials(:,1:currentMaxSingle,:) = DistInTrialsTemp;
                        DistOutTrials(:,1:currentMaxSingle,:) = DistOutTrialsTemp;
                        currentMaxSingle = maxDegree;
                    else
                        DistInTrials(trial,1:maxDegree,:) = DistIn;
                        DistOutTrials(trial,1:maxDegree,:) = DistOut;
                    end
                end
                DistInMean = reshape(mean(DistInTrials),[currentMaxSingle,L]);
                DistOutMean = reshape(mean(DistOutTrials),[currentMaxSingle,L]);
                %Smooth with a Dirichlet prior
                DistInMeanSmoothed = DirichletPrior1D(DistInMean);
                DistOutMeanSmoothed = DirichletPrior1D(DistOutMean);
                JSD = CalculateJSD(DistInMeanSmoothed,DistOutMeanSmoothed);
                writematrix(JSD,'ParticipantLevelAnalysisLevel1.txt','WriteMode','append')
                
                if currentMaxSingle>currentMaxSinglePPT
                    DistInTemp = DistInPPT;
                    DistOutTemp = DistOutPPT;
                    DistInPPT = zeros(currentMaxSingle,L);
                    DistOutPPT = zeros(currentMaxSingle,L);
                    DistInPPT(1:currentMaxSinglePPT,:) = DistInTemp;
                    DistOutPPT(1:currentMaxSinglePPT,:) = DistOutTemp;
                    clear DistInTemp
                    clear DistOutTemp
                    DistInPPT = DistInPPT + DistInMeanSmoothed;
                    DistOutPPT = DistOutPPT + DistOutMeanSmoothed;
                    currentMaxSinglePPT = currentMaxSingle;
                else
                    DistInPPT(1:currentMaxSingle,:) = DistInPPT(1:currentMaxSingle,:)+DistInMeanSmoothed;
                    DistOutPPT(1:currentMaxSingle,:) = DistOutPPT(1:currentMaxSingle,:)+DistOutMeanSmoothed;
                end
            end
            DistInPPT = DistInPPT./participants;
            DistOutPPT = DistOutPPT./participants;
            %1d divergence
            clear JSD
            JSD = CalculateJSD(DistInPPT,DistOutPPT);
            writematrix(JSD,'ConditionLevelAnalysisLevel1.txt','WriteMode','append')
        end
    end
end