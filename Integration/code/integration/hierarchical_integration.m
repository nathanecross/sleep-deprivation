function HI = hierarchical_integration(fmri, network_ids, mask_rois, ...
    opt_inference, verbose)
% HIERARCHICAL_INTEGRATION - Estimates total, inter and intra integration
%
%
% Examples:
% % Creates artificial data:
% n_samples = 1000;
% C = chol([2, 1.5; 1.5, 3]);
% n_rois = 2;
% network_ids = 1 : 3; % 3 networks
% mask_rois = randi([min(network_ids), max(network_ids)], [n_rois, 1]);
% fmri = {...
%     randn(n_samples, n_rois) * C; ...
%     randn(n_samples, n_rois) * C; ...
%     randn(n_samples, n_rois) * C; ...
%     randn(n_samples, n_rois) * C}; % 4 runs
%
% % Estimate integration measures #1
% HI = hierarchical_integration(fmri, network_ids, mask_rois, ...
%     struct('method', 'standard-inference', 'n_samplings', 1000), true)
%
% % Estimate integration measures #2
% HI = hierarchical_integration(fmri, network_ids, mask_rois, ...
%     struct('method', 'hierarchical-inference', 'n_samplings', 1000, ...
%     'nu_max', 500), true)
%
%
% Inputs:
% * fmri: a cell of matrices with dimensions [(number observations), (number
%   regions)]; each matrix represents one fMRI run
% * network_ids: a vector of non-null integers of dimension [(number networks)]
% * mask_rois: a vector of scalars of dimension [(number regions)] where each
%   ROI has been assigned to one network
% * opt_inference: a structure specifying the options for inference on the
%   sample covariance with the following fields:
%   * method: a string {'standard-inference'; 'hierarchical-inference'}
%   * n_samplings: a scalar (e.g. 1000)
%   * nu_max: a non-null positive integrer (e.g. 500) representing the greatest
%     degrees of freedom for the group only when 'method' is
%     'hierarchical-inference'            NC - !! Must be greater than (number observations) !!
%   * force_hi: a boolean (false by default) to force 'hierarchical-inference'
%     even if only a single fMRI run can be analyzed
% * verbose: a boolean indicating if to print additional informations
%
% Output:
% * HI: a cell of structures where each structure gives total, inter and intra
%   integration; each structure corresponds to one fMRI run in the same order;
%   if 'hierarchical-inference' was chosen, then the last structure gives the
%   group measures
%
% 
% Maintainer: @frazavi, Fatemeh Razavipour
% Last revision: 2020, June.


% Copyright (c) 2020 Fatemeh Razavipour
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.


% - Data covariance should be symmetric positive definite
% mask: filters ROIs from the specified networks
mask = ismember(mask_rois, network_ids);
[S, w, discard] = sample_cov(fmri, mask);

% - Sample covariance matrix
if ~isfield(opt_inference, 'method')
    error(['Missing the field ''method'' in ''opt_inference'', ', ...
        'type ''help %s'' for more info.'], mfilename)
end
if ~isfield(opt_inference, 'force_hi')
    opt_inference.force_hi = false;
end
if (size(S, 1) == 1) && ~opt_inference.force_hi && ...
        ~strcmp(opt_inference.method, 'standard-inference')
    warning(['Performing ''standard-inference'', to change this behavior ', ...
        'set to true ''opt_inference.force_hi'', ', ...
        'type ''help %s'' for more info.'], mfilename)
    opt_inference.method = 'standard-inference';
end
C = estimate_cov(S, w, opt_inference, verbose);

% - Integration measures: total, inter, intra
mask = mask_rois(mask);
HI = estimate_int(C, discard, network_ids, mask, verbose);
end



%% ------------ LOCAL FUNCTIONS ------------------------------------------------

function [S, w, discard] = sample_cov(fmri, mask_rois)
if ~iscell(fmri)
    error('Input ''fmri'' in not a cell, type ''help %s'' for more info.', ...
        mfilename)
elseif ~isvector(mask_rois)
    error(['Input ''mask_rois'' in not a vector, type ''help %s''', ...
        'for more info.'], mfilename)
elseif all(~mask_rois)
    error(['ROIs in ''mask_rois'' are not assigned to any network, ', ...
        'type ''help %s'' for more info.'], mfilename)
end

S = cell(numel(fmri), 1);
w = zeros(numel(fmri), 1);
discard = false(numel(fmri), 1);
for d = 1 : size(S, 1)
    if size(fmri{d}, 2) ~= numel(mask_rois)
        error(['fMRI data #%d does not have %d ROIs ', ...
            '(size of ''mask_rois''), type ''help %s'' for more info.'], d, ...
            numel(mask_rois), mfilename)
    else
        w(d) = size(fmri{d}, 1) - 1;
        S{d} = w(d) * cov(fmri{d}(:, mask_rois));
        try
            chol(S{d});
        catch
            discard(d) = true;
            warning(['Data covariance #%d is not symmetric positive ', ...
                'definite, ignoring it...'])
        end
    end
end
S(discard) = [];
w(discard) = [];
end


function C = estimate_cov(S, w, opt, verbose)
if strcmp(opt.method, 'standard-inference')
    if verbose
        fprintf('- Sampling from a Wishart distribution...\n')
    end
    C = standard_inference(S, w, opt);
    
elseif strcmp(opt.method, 'hierarchical-inference')
    if verbose
        fprintf('- Gibbs sampling on a hierarchical model...\n')
    end
    C = hierarchical_inference(S, w, opt);
    
else
    error('Cannot estimate the covariance, unknown method:\n%s', opt.method)
    
end
end


function C = standard_inference(S, w, opt)
if ~isfield(opt, 'n_samplings')
    error(['Missing the field ''n_samplings'' in ''opt'', type ''help %s''', ...
        'for more info.'], mfilename)
end
n_data = size(S, 1);
C.Sigma = cell(opt.n_samplings, n_data);
for d = 1 : n_data
    D = inv(chol(S{d})).';
    for k = 1 : opt.n_samplings
        C.Sigma{k, d} = iwishrnd(S{d}, w(d), D);
    end
end
end


function C = hierarchical_inference(S, w, opt)
if ~isfield(opt, 'n_samplings')
    error(['Missing the field ''n_samplings'' in ''opt'', type ''help %s''', ...
        'for more info.'], mfilename)
elseif ~isfield(opt, 'nu_max')
    error(['Missing the field ''nu_max'' in ''opt'', type ''help %s''', ...
        'for more info.'], mfilename)
end

% Initialisation
n_data = size(S, 1);
C.Sigma = cell(opt.n_samplings, n_data);
C.Sigma_0 = cell(opt.n_samplings, 1);
C.nu0 = zeros(opt.n_samplings, 1);

C.Sigma_0{1} = 0;
for d = 1 : n_data
    C.Sigma{1, d} = S{d};
    C.Sigma_0{1} = C.Sigma_0{1} + S{d};
end
C.Sigma_0{1} = C.Sigma_0{1} / n_data;

C.nu0(1) = size(S{1}, 1) + 1;
nu = (size(S{1}, 1) : opt.nu_max).';

% Iterations
for k = 2 : (2 * opt.n_samplings)
    % group covariance matrix
    s = 0;
    for d = 1 : n_data
        s = s + inv(C.Sigma{k - 1, d});
    end
    is = abs(round(inv(s),6));                                        %% NC- ensure inv matrix is symmetric
    isj = nearestSPD(is);
    C.Sigma_0{k} = wishrnd(isj, n_data * C.nu0(k - 1));
    % individual covariance matrices
    for d = 1 : n_data
        C.Sigma{k, d} = iwishrnd(S{d} + C.Sigma_0{k}, w(d) + C.nu0(k - 1));
    end
    
    % degrees of freedom for the group
    p = 0;
    for d = 1 : n_data
        try chol(C.Sigma_0{k});
            %disp('Matrix is symmetric positive definite.')
        catch ME
            %disp('Matrix is not symmetric positive definite')
            C.Sigma_0{k} = nearestSPD(C.Sigma_0{k}); 
        end
        try chol(C.Sigma{k, d});
            %disp('Matrix is symmetric positive definite.')
        catch ME
            %disp('Matrix is not symmetric positive definite')
            C.Sigma{k, d} = nearestSPD(C.Sigma{k, d}); 
        end
        size(cholcov(C.Sigma{k, d}));
        size(cholcov(C.Sigma_0{k}));
        p = p + ln_inv_wishart(C.Sigma{k, d}, C.Sigma_0{k}, nu);
    end
    p = exp(p - max(p));
    C.nu0(k) = [nu(find(cumsum(p / sum(p)) >= rand, 1))];
end
k_burn_in = opt.n_samplings + 1;
C.Sigma = C.Sigma(k_burn_in : end, :);
C.Sigma_0 = C.Sigma_0(k_burn_in : end);
C.nu0 = C.nu0(k_burn_in : end);
end


function p = ln_inv_wishart(x, Sigma, nu)
n_rois = size(Sigma, 1);
p = sum(gammaln((repmat(nu, 1, n_rois) - ...
    repmat(0 : (n_rois - 1), size(nu, 1), 1)) / 2), 2);
p = (nu * log_det_cov(Sigma) - (nu + n_rois + 1) * log_det_cov(x) - ...
    trace(Sigma / x) - nu * n_rois * log(2) - n_rois * (n_rois - 1) * ...
    log(pi) / 2) / 2 - p;
end


function d = log_det_cov(C)
    %d = 2 * log(det(cholcov(C)));
    d = 2 * real(logdet(cholcov(C)));             %% NC - logdet ensures doesn't return Inf for large matrices
end


function HI = estimate_int(C, discard, network_ids, mask_rois, verbose)
if verbose
    fprintf('- Integration measures...\n')
end

HI = cell(size(discard, 1), 1);
for d = 1 : size(discard, 1)
    if discard(d)
        continue
    end
    HI{d} = integragtion_measures(C.Sigma(:, d), network_ids, mask_rois);
end
if isfield(C, 'Sigma_0')
    HI{end + 1} = integragtion_measures(C.Sigma_0, network_ids, mask_rois);
end
end


function HI = integragtion_measures(Sigma, network_ids, mask_rois)
n_samples = size(Sigma, 1);
n_networks = numel(network_ids);
int_total = zeros(n_samples, 1);
int_inter = zeros(n_samples, 1);
int_intra = zeros(n_networks, n_networks, n_samples);
int_matrix = zeros(n_networks);
for k = 1 : n_samples
    for k_1 = 1 : n_networks
        mask_1 = (mask_rois == network_ids(k_1));
        S1 = Sigma{k}(mask_1, mask_1);
        d1 = log_det_cov(S1);
        
        int_matrix(k_1, k_1) = sum(log(diag(S1))) - d1;
        for k_2 = (k_1+1) : n_networks
            mask_2 = (mask_rois == network_ids(k_2));
            mask_1_2 = mask_1 | mask_2;
            
            int_matrix(k_1, k_2) = d1 + ...
                log_det_cov(Sigma{k}(mask_2, mask_2)) - ...
                log_det_cov(Sigma{k}(mask_1_2, mask_1_2));
            int_matrix(k_2, k_1) = int_matrix(k_1, k_2);
        end
    end
    int_intra(:, :, k) = int_matrix / 2;
    int_total(k) = (sum(log(diag(Sigma{k}))) - log_det_cov(Sigma{k})) / 2;
    int_inter(k) = int_total(k) - sum(diag(int_intra(:, :, k)));
end
if (n_samples > 1)
    HI = struct(...
        ... % total
        'int_total', struct(...
        'samples', int_total, ...
        'mean', mean(int_total), ...
        'var', std(int_total)), ...
        ... % inter
        'int_inter', struct(...
        'samples', int_inter, ...
        'mean', mean(int_inter), ...
        'var', std(int_inter)), ...
        ... % intra
        'int_intra', struct(...
        'samples', int_intra, ...
        'mean', mean(int_intra, 3), ...
        'var', var(int_intra, 1, 3)));
else
    HI = struct(...
        'int_total', int_total, ...
        'int_inter', int_inter, ...
        'int_intra', int_intra);
end
end
