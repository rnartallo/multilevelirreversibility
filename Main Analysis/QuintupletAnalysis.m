% load the data
% add helpers to the path
%% Analysis of order 5
participants =51;
conditions =1;
trials = 15;
L=6;
ts_length = 1026;

for cond =1:conditions
    writematrix('Condition: ' + string(cond),'ParticipantLevelAnalysisLevel5E.txt','WriteMode','append')
    writematrix('Condition: ' + string(cond),'ConditionLevelAnalysisLevel5E.txt','WriteMode','append')

    for alpha=1:L-4
        alpha
        for beta=alpha+1:L-3
            beta
            for theta=beta+1:L-2
                for delta=theta+1:L-1
                    for epsilon=delta+1:L
                        tuple =[alpha, beta, theta, delta, epsilon];
                        writematrix(string(alpha) + ',' + string(beta) + ',' + string(theta) + ',' + string(delta) + ',' +string(epsilon),'ParticipantLevelAnalysisLevel5E.txt','WriteMode','append')
                        QuintDistInPPT = zeros(1,1,1,1,1);
                        QuintDistOutPPT = zeros(1,1,1,1,1);
                        currentMaxQuintPPT = 1;
                        for participant =1:participants
                            QuintDistIn = zeros(1,1,1,1,1);
                            QuintDistOut = zeros(1,1,1,1,1);
                            currentMaxQuint =1;
                            for trial =1:trials
                                %Extract the desired time-series
                                X = MEGData{participant}(:,:,trial);
                                Y = X(tuple,:);
                                %Build the network
                                A = MultiplexVisibilityNetwork(Y);
                                % Build a 5-d distribution for each trial
                                [K_In,K_Out] = SignedDegrees(A);
                                maxDegree = min(75,max(max(max(K_In)),max(max(K_Out))));
                                if maxDegree > currentMaxQuint
                                    QuintDistIn = padarray(QuintDistIn,[maxDegree-currentMaxQuint,maxDegree-currentMaxQuint,maxDegree-currentMaxQuint,maxDegree-currentMaxQuint,maxDegree-currentMaxQuint],0,'post');
                                    QuintDistOut = padarray(QuintDistOut,[maxDegree-currentMaxQuint,maxDegree-currentMaxQuint,maxDegree-currentMaxQuint,maxDegree-currentMaxQuint,maxDegree-currentMaxQuint],0,'post');
                                    currentMaxQuint = maxDegree;
                                end
                                for node =1:ts_length
                                    if node<ts_length
                                        D_vector_in = min(75,K_In(node,:));
                                        QuintDistIn(D_vector_in(1),D_vector_in(2),D_vector_in(3),D_vector_in(4),D_vector_in(5)) = QuintDistIn(D_vector_in(1),D_vector_in(2),D_vector_in(3),D_vector_in(4),D_vector_in(5))+1;
                                    end
                                    if node >1
                                        D_vector_out = min(75,K_Out(node,:));
                                        QuintDistOut(D_vector_out(1),D_vector_out(2),D_vector_out(3),D_vector_out(4),D_vector_out(5)) = QuintDistOut(D_vector_out(1),D_vector_out(2),D_vector_out(3),D_vector_out(4),D_vector_out(5))+1;
                                    end
                                end
                            end
                            %Trial average the distribution 
                            QuintDistIn = QuintDistIn./(trials*ts_length);
                            QuintDistOut = QuintDistOut./(trials*ts_length);
    
                            %Calculate participant level divergence
                            JSD= CalculateJSD5(QuintDistOut,QuintDistIn);
                            writematrix(JSD,'ParticipantLevelAnalysisLevel5E.txt','WriteMode','append')
                            if currentMaxQuint > currentMaxQuintPPT
                                QuintDistInPPT = padarray(QuintDistInPPT,[currentMaxQuint-currentMaxQuintPPT,currentMaxQuint-currentMaxQuintPPT,currentMaxQuint-currentMaxQuintPPT,currentMaxQuint-currentMaxQuintPPT,currentMaxQuint-currentMaxQuintPPT],0,'post');
                                QuintDistOutPPT = padarray(QuintDistOutPPT,[currentMaxQuint-currentMaxQuintPPT,currentMaxQuint-currentMaxQuintPPT,currentMaxQuint-currentMaxQuintPPT,currentMaxQuint-currentMaxQuintPPT,currentMaxQuint-currentMaxQuintPPT],0,'post');
                                QuintDistInPPT = QuintDistInPPT + QuintDistIn;
                                QuintDistOutPPT = QuintDistOutPPT + QuintDistOut;
                                currentMaxQuintPPT = currentMaxQuint;
                            else
                                QuintDistInPPT(1:currentMaxQuint,1:currentMaxQuint,1:currentMaxQuint,1:currentMaxQuint,1:currentMaxQuint) = QuintDistInPPT(1:currentMaxQuint,1:currentMaxQuint,1:currentMaxQuint,1:currentMaxQuint,1:currentMaxQuint) + QuintDistIn;
                                QuintDistOutPPT(1:currentMaxQuint,1:currentMaxQuint,1:currentMaxQuint,1:currentMaxQuint,1:currentMaxQuint) = QuintDistOutPPT(1:currentMaxQuint,1:currentMaxQuint,1:currentMaxQuint,1:currentMaxQuint,1:currentMaxQuint) + QuintDistOut;
                            end
                        end
                        QuintDistInPPT=QuintDistInPPT./participants;
                        QuintDistOutPPT=QuintDistOutPPT./participants;
                        %Calculate condition level divergence
                        JSD= CalculateJSD5(QuintDistInPPT,QuintDistOutPPT);
                        writematrix(JSD,'ConditionLevelAnalysisLevel5E.txt','WriteMode','append')
                    end
                end
            end
        end
    end
end
