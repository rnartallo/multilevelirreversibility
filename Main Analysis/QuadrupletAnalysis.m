% load the data
% add helpers to the path
%% Analysis of order 4:
participants =51;
conditions =1;
trials = 15;
L=6;
ts_length = 1026;

for cond =1:conditions
    writematrix('Condition: ' + string(cond),'ParticipantLevelAnalysisLevel4E.txt','WriteMode','append')
    writematrix('Condition: ' + string(cond),'ConditionLevelAnalysisLevel4E.txt','WriteMode','append')
    
    %ORDER 4

    for alpha=1:L-3
        alpha
        for beta=alpha+1:L-2
            beta
            for theta=beta+1:L-1
                for delta=theta+1:L
                    tuple =[alpha, beta, theta, delta];
                    writematrix(string(alpha) + ',' + string(beta) + ',' + string(theta) + ',' + string(delta),'ParticipantLevelAnalysisLevel4E.txt','WriteMode','append')
                    QuadDistInPPT = zeros(1,1,1,1);
                    QuadDistOutPPT = zeros(1,1,1,1);
                    currentMaxQuadPPT = 1;
                    for participant =1:participants
                        QuadDistIn = zeros(1,1,1,1);
                        QuadDistOut = zeros(1,1,1,1);
                        currentMaxQuad =1;
                        for trial =1:trials
                            %Extract the desired time-series
                            X = MEGData{participant}(:,:,trial);
                            Y = X(tuple,:);
                            %Build the network
                            A = MultiplexVisibilityNetwork(Y);
                            % Build a 4-d
                            [K_In,K_Out] = SignedDegrees(A);
                            maxDegree = max(max(max(K_In)),max(max(K_Out)));
                            if maxDegree > currentMaxQuad
                                QuadDistInTemp = QuadDistIn;
                                QuadDistOutTemp = QuadDistOut;
                                QuadDistIn = zeros(maxDegree,maxDegree,maxDegree,maxDegree);
                                QuadDistOut = zeros(maxDegree,maxDegree,maxDegree,maxDegree);
                                QuadDistIn(1:currentMaxQuad,1:currentMaxQuad,1:currentMaxQuad,1:currentMaxQuad) = QuadDistInTemp;
                                QuadDistOut(1:currentMaxQuad,1:currentMaxQuad,1:currentMaxQuad,1:currentMaxQuad) = QuadDistOutTemp;
                                currentMaxQuad = maxDegree;
                            end

                            for node =1:ts_length
                                if node<ts_length
                                    D_vector_in = K_In(node,:);
                                    QuadDistIn(D_vector_in(1),D_vector_in(2),D_vector_in(3),D_vector_in(4)) = QuadDistIn(D_vector_in(1),D_vector_in(2),D_vector_in(3),D_vector_in(4))+1;
                                end
                                if node>1
                                    D_vector_out = K_Out(node,:);
                                    QuadDistOut(D_vector_out(1),D_vector_out(2),D_vector_out(3),D_vector_out(4)) = QuadDistOut(D_vector_out(1),D_vector_out(2),D_vector_out(3),D_vector_out(4))+1;
                                end
                            end

                        end
                        %Trial average the distribution 
                        QuadDistIn = QuadDistIn./(trials*(ts_length-1));
                        QuadDistOut = QuadDistOut./(trials*(ts_length-1));

                        %Calculate participant level divergence
                        JSD= CalculateJSD4(QuadDistOut,QuadDistIn);
                        writematrix(JSD,'ParticipantLevelAnalysisLevel4E.txt','WriteMode','append')
                        if currentMaxQuad > currentMaxQuadPPT
                            QuadDistInTemp = QuadDistInPPT;
                            QuadDistOutTemp = QuadDistOutPPT;
                            QuadDistInPPT = zeros(currentMaxQuad,currentMaxQuad,currentMaxQuad,currentMaxQuad);
                            QuadDistOutPPT = zeros(currentMaxQuad,currentMaxQuad,currentMaxQuad,currentMaxQuad);
                            QuadDistInPPT(1:currentMaxQuadPPT,1:currentMaxQuadPPT,1:currentMaxQuadPPT,1:currentMaxQuadPPT) = QuadDistInTemp;
                            QuadDistOutPPT(1:currentMaxQuadPPT,1:currentMaxQuadPPT,1:currentMaxQuadPPT,1:currentMaxQuadPPT) = QuadDistOutTemp;
                            QuadDistInPPT = QuadDistInPPT + QuadDistIn;
                            QuadDistOutPPT = QuadDistOutPPT + QuadDistOut;
                            currentMaxQuadPPT = currentMaxQuad;
                        else
                            QuadDistInPPT(1:currentMaxQuad,1:currentMaxQuad,1:currentMaxQuad,1:currentMaxQuad) = QuadDistInPPT(1:currentMaxQuad,1:currentMaxQuad,1:currentMaxQuad,1:currentMaxQuad) + QuadDistIn;
                            QuadDistOutPPT(1:currentMaxQuad,1:currentMaxQuad,1:currentMaxQuad,1:currentMaxQuad) = QuadDistOutPPT(1:currentMaxQuad,1:currentMaxQuad,1:currentMaxQuad,1:currentMaxQuad) + QuadDistOut;
                        end
                    end
                    QuadDistInPPT=QuadDistInPPT./participants;
                    QuadDistOutPPT=QuadDistOutPPT./participants;
                    %Calculate condition level divergence
                    JSD= CalculateJSD4(QuadDistInPPT,QuadDistOutPPT);
                    writematrix(JSD,'ConditionLevelAnalysisLevel4E.txt','WriteMode','append')
                end
            end
        end
    end
end