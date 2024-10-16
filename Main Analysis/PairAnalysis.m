% load the data
% add helpers to the path
%% Analysis of order 2:
participants =51;
conditions =1;
trials = 15;
L=6;
ts_length = 1026;

for cond =1:conditions
    writematrix('Condition: ' + string(cond),'ParticipantLevelAnalysisLevel2.txt','WriteMode','append')
    writematrix('Condition: ' + string(cond),'ConditionLevelAnalysisLevel2.txt','WriteMode','append')
    
    %ORDER 2

    for alpha=1:L-1
        alpha
        for beta=alpha+1:L
            tuple =[alpha, beta];
            writematrix(string(alpha) + ',' + string(beta),'ParticipantLevelAnalysisLevel2.txt','WriteMode','append')
            JointDistInPPT = zeros(1,1);
            JointDistOutPPT = zeros(1,1);
            currentMaxJointPPT = 1;
            for participant =1:participants
                JointDistIn = zeros(1,1);
                JointDistOut = zeros(1,1);
                currentMaxJoint =1;
                for trial =1:trials
                    %Extract the desired time-series
                    X = MEGData{participant}(:,:,trial);
                    Y = X(tuple,:);
                    %Build the network
                    A = MultiplexVisibilityNetwork(Y);
                    % Build a 2-d
                    [K_In,K_Out] = SignedDegrees(A);
                    maxDegree = max(max(max(K_In)),max(max(K_Out)));
                    if maxDegree > currentMaxJoint
                        JointDistInTemp = JointDistIn;
                        JointDistOutTemp = JointDistOut;
                        JointDistIn = zeros(maxDegree,maxDegree);
                        JointDistOut = zeros(maxDegree,maxDegree);
                        JointDistIn(1:currentMaxJoint,1:currentMaxJoint) = JointDistInTemp;
                        JointDistOut(1:currentMaxJoint,1:currentMaxJoint) = JointDistOutTemp;
                        currentMaxJoint = maxDegree;
                    end

                    for node =1:ts_length
                        if node<ts_length
                            D_vector_in = K_In(node,:);
                            JointDistIn(D_vector_in(1),D_vector_in(2)) = JointDistIn(D_vector_in(1),D_vector_in(2))+1;
                        end
                        if node>1
                            D_vector_out = K_Out(node,:);
                            JointDistOut(D_vector_out(1),D_vector_out(2)) = JointDistOut(D_vector_out(1),D_vector_out(2))+1;
                        end
                    end

                end
                %Trial average the distribution 
                JointDistIn = JointDistIn./(trials*(ts_length-1));
                JointDistOut = JointDistOut./(trials*(ts_length-1));

                %Calculate participant level divergence
                JSD= CalculateJSDJoint(JointDistOut,JointDistIn);
                writematrix(JSD,'ParticipantLevelAnalysisLevel2.txt','WriteMode','append')
                if currentMaxJoint > currentMaxJointPPT
                    JointDistInTemp = JointDistInPPT;
                    JointDistOutTemp = JointDistOutPPT;
                    JointDistInPPT = zeros(currentMaxJoint,currentMaxJoint);
                    JointDistOutPPT = zeros(currentMaxJoint,currentMaxJoint);
                    JointDistInPPT(1:currentMaxJointPPT,1:currentMaxJointPPT) = JointDistInTemp;
                    JointDistOutPPT(1:currentMaxJointPPT,1:currentMaxJointPPT) = JointDistOutTemp;
                    JointDistInPPT = JointDistInPPT + JointDistIn;
                    JointDistOutPPT = JointDistOutPPT + JointDistOut;
                    currentMaxJointPPT = currentMaxJoint;
                else
                    JointDistInPPT(1:currentMaxJoint,1:currentMaxJoint) = JointDistInPPT(1:currentMaxJoint,1:currentMaxJoint) + JointDistIn;
                    JointDistOutPPT(1:currentMaxJoint,1:currentMaxJoint) = JointDistOutPPT(1:currentMaxJoint,1:currentMaxJoint) + JointDistOut;
                end
            end
            JointDistInPPT=JointDistInPPT./participants;
            JointDistOutPPT=JointDistOutPPT./participants;
            %Calculate condition level divergence
            JSD= CalculateJSDJoint(JointDistInPPT,JointDistOutPPT);
            writematrix(JSD,'ConditionLevelAnalysisLevel2.txt','WriteMode','append')
        end
    end
end
