%% load data and add helpers
repeats = 20;
sample_sizes =[15,12,9];
L=6;
ts_length = 1026;
for r=1:repeats
    r
    sample{1}=1:15;
    for i=2:4
        sample{i} = randsample(sample{i-1},sample_sizes(i));
    end
    parfor alpha=1:L
        alpha
        for i=1:3
            DistIn = zeros(1);
            DistOut = zeros(1);
            currentMax =1;
            no_trials = sample_sizes(i)*51;
            for p=1:51
                mvts_all_trials = MEGData{p};
                for trial =1:sample_sizes(i)
                    X = mvts_all_trials(:,:,sample{i}(trial));
                    Y = X(alpha,:);
                    A = MultiplexVisibilityNetwork(Y);
                    [K_In,K_Out] = SignedDegrees(A);
                    maxDegree = max(max(max(K_In)),max(max(K_Out)));
                    if maxDegree > currentMax
                        DistInTemp = DistIn;
                        DistOutTemp = DistOut;
                        DistIn = zeros(maxDegree);
                        DistOut = zeros(maxDegree);
                        DistIn(1:currentMax) = DistInTemp;
                        DistOut(1:currentMax) = DistOutTemp;
                        currentMax = maxDegree;
                    end
                    for node =1:ts_length
                        if node<ts_length
                            D_vector_in = K_In(node);
                            DistIn(D_vector_in) = DistIn(D_vector_in)+1;
                        end
                        if node>1
                            D_vector_out = K_Out(node);
                            DistOut(D_vector_out(1)) = DistOut(D_vector_out(1))+1;
                        end
                    end
                end
            end
            %Trial average the distribution 
            DistIn = DistIn./(no_trials*(ts_length-1));
            DistOut = DistOut./(no_trials*(ts_length-1));

            %Calculate divergence
            DistInMeanSmoothed = DirichletPrior1D(DistIn);
            DistOutMeanSmoothed = DirichletPrior1D(DistOut);
            JSD = CalculateJSD(DistInMeanSmoothed,DistOutMeanSmoothed);
            writematrix(string(JSD) + ',' + string(alpha) + ',' + string(r)+ ',' + string(i),'Finite_Correction1.txt','WriteMode','append')
        end
    end
end
