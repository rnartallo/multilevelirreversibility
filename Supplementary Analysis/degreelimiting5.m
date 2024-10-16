% load the data add the helpers
%% Analysis of order 4:6
participants =51;
conditions =5;
trials = 15;
L=6;
ts_length = 1026;
cond=1;
for degprop=0.1:0.1:0.9
    writematrix('Degree proportion: ' + string(degprop),'DegreeLimited5.txt','WriteMode','append')

    for alpha=1:L-4
        alpha
        for beta=alpha+1:L-3
            beta
            for theta=beta+1:L-2
                for delta=theta+1:L-1
                    for epsilon=delta+1:L
                        tuple =[alpha, beta, theta, delta, epsilon];
                        QuintDistInPPT = zeros(1,1,1,1,1);
                        QuintDistOutPPT = zeros(1,1,1,1,1);
                        currentMaxQuintPPT = 1;
                        for participant =1:participants
                            QuintDistIn = zeros(1,1,1,1,1);
                            QuintDistOut = zeros(1,1,1,1,1);
                            currentMaxQuint =1;
                            for trial =1:trials
                                %Extract the desired time-series
                                X = SampledData{participant}.Condition{cond}(:,:,trial);
                                Y = X(tuple,:);
                                %Build the network
                                A = MultiplexVisibilityNetwork(Y);
                                % Build a 4-d,5-d,6-d distribution for each trial
                                [K_In,K_Out] = SignedDegrees(A);
                                
                                maxDegree = min(75,max(max(max(K_In)),max(max(K_Out))));
                                maxDegree = max(floor(degprop*maxDegree),1);

                                if maxDegree > currentMaxQuint
                                    QuintDistIn = padarray(QuintDistIn,[maxDegree-currentMaxQuint,maxDegree-currentMaxQuint,maxDegree-currentMaxQuint,maxDegree-currentMaxQuint,maxDegree-currentMaxQuint],0,'post');
                                    QuintDistOut = padarray(QuintDistOut,[maxDegree-currentMaxQuint,maxDegree-currentMaxQuint,maxDegree-currentMaxQuint,maxDegree-currentMaxQuint,maxDegree-currentMaxQuint],0,'post');
                                    currentMaxQuint = maxDegree;
                                end
                                for node =1:ts_length
                                    if node<ts_length
                                        D_vector_in = min(maxDegree,K_In(node,:));
                                        QuintDistIn(D_vector_in(1),D_vector_in(2),D_vector_in(3),D_vector_in(4),D_vector_in(5)) = QuintDistIn(D_vector_in(1),D_vector_in(2),D_vector_in(3),D_vector_in(4),D_vector_in(5))+1;
                                    end
                                    if node >1
                                        D_vector_out = min(maxDegree,K_Out(node,:));
                                        QuintDistOut(D_vector_out(1),D_vector_out(2),D_vector_out(3),D_vector_out(4),D_vector_out(5)) = QuintDistOut(D_vector_out(1),D_vector_out(2),D_vector_out(3),D_vector_out(4),D_vector_out(5))+1;
                                    end
                                end
                            end
                            %Trial average the distribution 
                            QuintDistIn = QuintDistIn./(trials*ts_length);
                            QuintDistOut = QuintDistOut./(trials*ts_length);
                               
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
                        writematrix(JSD,'DegreeLimited5.txt','WriteMode','append')
                    end
                end
            end
        end
    end
end
