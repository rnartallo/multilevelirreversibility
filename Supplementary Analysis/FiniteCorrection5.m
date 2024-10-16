%% load the data and add the helpers

repeats = 20;
sample_sizes =[15,12,9];
L=6;
ts_length = 1026;
for r=15:repeats
    r
    sample{1}=1:15;
    for i=2:4
        sample{i} = randsample(sample{i-1},sample_sizes(i));
    end
    for alpha=1:L-4
        alpha
        for beta=(alpha+1):L-3
            for theta=(beta+1):L-2
                for delta=(theta+1):L-1
                    for eps=(delta+1):L
                        tuple =[alpha,beta,theta,delta,eps];
                        for i=1:3
                            DistIn = zeros(1,1,1,1,1);
                            DistOut = zeros(1,1,1,1,1);
                            currentMax =1;
                            no_trials = sample_sizes(i)*51;
                            for p=1:51
                                mvts_all_trials = MEGData{p};
                                for trial =1:sample_sizes(i)
                                    X = mvts_all_trials(:,:,sample{i}(trial));
                                    Y = X(tuple,:);
                                    A = MultiplexVisibilityNetwork(Y);
                                    [K_In,K_Out] = SignedDegrees(A);
                                    maxDegree = max(max(max(K_In)),max(max(K_Out)));
                                    if maxDegree>45
                                        maxDegree=45;
                                    end
                                    if maxDegree > currentMax
                                        DistIn = padarray(DistIn,[maxDegree-currentMax,maxDegree-currentMax,maxDegree-currentMax,maxDegree-currentMax,maxDegree-currentMax],0,'post');
                                        DistOut = padarray(DistOut,[maxDegree-currentMax,maxDegree-currentMax,maxDegree-currentMax,maxDegree-currentMax,maxDegree-currentMax],0,'post');
                                        currentMax = maxDegree;
                                    end
                                    for node =1:ts_length
                                        if node<ts_length
                                            D_vector_in = min(45,K_In(node,:));
                                            DistIn(D_vector_in(1),D_vector_in(2),D_vector_in(3),D_vector_in(4),D_vector_in(5)) = DistIn(D_vector_in(1),D_vector_in(2),D_vector_in(3),D_vector_in(4),D_vector_in(5))+1;
                                            clear D_vector_in
                                        end
                                        if node>1
                                            D_vector_out = min(45,K_Out(node,:));
                                            DistOut(D_vector_out(1),D_vector_out(2),D_vector_out(3),D_vector_out(4),D_vector_out(5)) = DistOut(D_vector_out(1),D_vector_out(2),D_vector_out(3),D_vector_out(4),D_vector_out(5))+1;
                                            clear D_vector_out
                                        end
                                    end
                                end
                            end
           
                        %Trial average the distribution 
                        DistIn = DistIn./(no_trials*(ts_length-1));
                        DistOut = DistOut./(no_trials*(ts_length-1));
        
                        %Calculate divergence
                        JSD= CalculateJSD5(DistOut,DistIn);
                        writematrix(string(JSD) + ',' + string(alpha) + ',' + string(beta) + ',' + string(theta) + ',' + string(delta) + ',' + string(eps) + ',' + string(r)+ ',' + string(i),'Finite_Correction5.txt','WriteMode','append')
                        end
                    end
                end
            end
        end
    end
end