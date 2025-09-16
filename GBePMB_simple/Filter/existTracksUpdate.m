function tracks_upd = existTracksUpdate(mbm,W,clusters,gating_matrix_d,measmodel,paras)
%
n_track = length(mbm.track);
n_clusters = size(clusters,2);
n_local_hypo = cellfun('length',mbm.track);
tracks_upd = cell(n_track,1);
for i = 1:n_track
    tracks_upd{i} = cell(n_local_hypo(i),1);
    for j = 1:n_local_hypo(i)
        %misdetection for the jth single target hypothesis under the
        %ith local hypothesis tree
        [l_missed,bern_missed] = BGGbern_miss(mbm.track{i}(j),measmodel,paras);
        tracks_upd{i}{j}.c = 0;
        tracks_upd{i}{j}.lik = l_missed;
        tracks_upd{i}{j}.bern = bern_missed;
        %measurement update for the jth single target hypothesis under
        %the ith local hypothesis tree
        for c = 1:n_clusters
            %check if the cth cluster is in the gate of the
            %corresponding single target hypothesis
            if sum(gating_matrix_d{i}(:,j)-clusters(:,c)<0) == 0
                [l_upd,bern_updated] = ...
                    BGGbern_upd(mbm.track{i}(j),W(:,clusters(:,c)),measmodel,paras);
                tracks_upd{i}{j}.c(end+1) = c;
                tracks_upd{i}{j}.lik(end+1) = l_upd;
                tracks_upd{i}{j}.bern(end+1) = bern_updated;
            end
        end
    end
end

end