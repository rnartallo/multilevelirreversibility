% load the data
% add helpers to the path
%% Analysis of order 3:
participants =51;
conditions =1;
trials = 15;
L=6;
ts_length = 1026;

for cond =1:conditions
    writematrix('Condition: ' + string(cond),'ParticipantLevelAnalysisLevel3.txt','WriteMode','append')
    writematrix('Condition: ' + string(cond),'ConditionLevelAnalysisLevel3.txt','WriteMode','append')
    
    %ORDER 3

    for alpha=1:L-2
        alpha
        for beta=alpha+1:L-1
            beta
            for theta=beta+1:L
                tuple =[alpha, beta, theta];
                writematrix(string(alpha) + ',' + string(beta) + ',' + string(theta),'ParticipantLevelAnalysisLevel3.txt','WriteMode','append')
                TripDistInPPT = zeros(1,1,1);
                TripDistOutPPT = zeros(1,1,1);
                currentMaxTripPPT = 1;
                for participant =1:participants
                    TripDistIn = zeros(1,1,1);
                    TripDistOut = zeros(1,1,1);
                    currentMaxTrip =1;
                    for trial =1:trials
                        %Extract the desired time-series
                        X = MEGData{participant}(:,:,trial);
                        Y = X(tuple,:);
                        %Build the network
                        A = MultiplexVisibilityNetwork(Y);
                        % Build a 3-d distribution
                        [K_In,K_Out] = SignedDegrees(A);
                        maxDegree = max(max(max(K_In)),max(max(K_Out)));
                        if maxDegree > currentMaxTrip
                            TripDistInTemp = TripDistIn;
                            TripDistOutTemp = TripDistOut;
                            TripDistIn = zeros(maxDegree,maxDegree,maxDegree);
                            TripDistOut = zeros(maxDegree,maxDegree,maxDegree);
                            TripDistIn(1:currentMaxTrip,1:currentMaxTrip,1:currentMaxTrip) = TripDistInTemp;
                            TripDistOut(1:currentMaxTrip,1:currentMaxTrip,1:currentMaxTrip) = TripDistOutTemp;
                            currentMaxTrip = maxDegree;
                        end

                        for node =1:ts_length
                            if node<ts_length
                                D_vector_in = K_In(node,:);
                                TripDistIn(D_vector_in(1),D_vector_in(2),D_vector_in(3)) = TripDistIn(D_vector_in(1),D_vector_in(2),D_vector_in(3))+1;
                            end
                            if node>1
                                D_vector_out = K_Out(node,:);
                                TripDistOut(D_vector_out(1),D_vector_out(2),D_vector_out(3)) = TripDistOut(D_vector_out(1),D_vector_out(2),D_vector_out(3))+1;
                            end
                        end

                    end
                    %Trial average the distribution 
                    TripDistIn = TripDistIn./(trials*(ts_length-1));
                    TripDistOut = TripDistOut./(trials*(ts_length-1));

                    %Calculate participant level divergence
                    JSD= CalculateJSD3(TripDistOut,TripDistIn);
                    writematrix(JSD,'ParticipantLevelAnalysisLevel3.txt','WriteMode','append')
                    if currentMaxTrip > currentMaxTripPPT
                        TripDistInTemp = TripDistInPPT;
                        TripDistOutTemp = TripDistOutPPT;
                        TripDistInPPT = zeros(currentMaxTrip,currentMaxTrip,currentMaxTrip);
                        TripDistOutPPT = zeros(currentMaxTrip,currentMaxTrip,currentMaxTrip);
                        TripDistInPPT(1:currentMaxTripPPT,1:currentMaxTripPPT,1:currentMaxTripPPT) = TripDistInTemp;
                        TripDistOutPPT(1:currentMaxTripPPT,1:currentMaxTripPPT,1:currentMaxTripPPT) = TripDistOutTemp;
                        TripDistInPPT = TripDistInPPT + TripDistIn;
                        TripDistOutPPT = TripDistOutPPT + TripDistOut;
                        currentMaxTripPPT = currentMaxTrip;
                    else
                        TripDistInPPT(1:currentMaxTrip,1:currentMaxTrip,1:currentMaxTrip) = TripDistInPPT(1:currentMaxTrip,1:currentMaxTrip,1:currentMaxTrip) + TripDistIn;
                        TripDistOutPPT(1:currentMaxTrip,1:currentMaxTrip,1:currentMaxTrip) = TripDistOutPPT(1:currentMaxTrip,1:currentMaxTrip,1:currentMaxTrip) + TripDistOut;
                    end
                end
                TripDistInPPT=TripDistInPPT./participants;
                TripDistOutPPT=TripDistOutPPT./participants;
                %Calculate condition level divergence
                JSD= CalculateJSD3(TripDistInPPT,TripDistOutPPT);
                writematrix(JSD,'ConditionLevelAnalysisLevel3.txt','WriteMode','append')
            end
        end
    end
end
