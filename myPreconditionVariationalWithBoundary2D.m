function v1 = myPreconditionVariationalWithBoundary2D(b, v, sigma, v0)

v1 = v0; % zeros(size(v));

idx_tr = (b ~= 0);
nb_label = max(b(:));

w = (b>0);

for ii = 1:3
     
    
     g1 = -2 * w .* (v - v1); %  / maxw;
     g = myGaussianLPF2D(g1, sigma); % g is soblev gradient 
     
%        v1 = v1 - 1./(sigma.^2) * g;
%        continue;
     
     
     idx_gnonzero = (g~=0);
     if isempty(find(idx_gnonzero, 1))
         disp 'converge';
         break; % converge
     end;
     
     idx_tr1 = idx_tr & idx_gnonzero; % points inside trajectories
     
     
     l = (v1 - v) ./ (g + 5);
     
     
     
     if find(isnan(l(idx_tr)))
         disp 'NaN found inside velocities of trajectories  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n';
         % l(isnan(l)) = 1;
         
%          l = medfilt2(l, [5, 5]); 
%          l(l > 1e3) = 1e3;
%          l(isnan(l)) = 1e3;
         
          
          la = l;
          la(l > 1e2) = NaN;
          l = inpaint_nans(la);
     end;
% 
%           l(l > 10) = 10;
%           l(l < 0) = 0;
%           l = medfilt2(l, [5, 5]); 
    
     
     idx_middle = (~idx_tr) & idx_gnonzero;
     nb_middle = sum(idx_middle(:));
     logwlist = zeros(nb_middle, nb_label);
     

     
     for jj = 1:nb_label;
         [D, L] = bwdist(b==jj);
         basis_value = l(L(idx_middle));
         dist_to_basis = D(idx_middle);
         % w = 1 / (2*pi*sigma) * exp(-1/(sigma^2) .* dist_to_basis.^2 );
         logwlist(:, jj) = -1 * log( (2*pi*sigma) ) + (-1/(sigma^2) .* dist_to_basis.^2 );
     end;
     
     lmiddle = zeros(nb_middle, 1);
     for jj = 1:nb_label
         [D, L] = bwdist(b==jj);
         basis_value = l(L(idx_middle));

         sw = zeros(nb_middle, 1);
         for kk = 1:nb_label
             sw = sw + exp( logwlist(:, kk) - logwlist(:, jj) );
         end;
         lmiddle = lmiddle + basis_value ./ sw;
     end;
     
     
     
     lrbf = zeros(size(l));
     lrbf(idx_middle) = lmiddle;
     lrbf(idx_tr1) = l(idx_tr1);
 
     v1 = v1 - lrbf .* g;
     
     
     
     
%     figure; clf;imagesc(v1);
%     figure; clf; imagesc(l);
end;

return;
