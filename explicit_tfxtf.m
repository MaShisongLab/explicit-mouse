classdef explicit_tfxtf
	properties
		beta =  [];   %  intercept + tf (row)  x tf (col)
		beta_pvalue = [];  % intercept + tf (row)  x tf (col)
		TF_name = [];
		NRMSE = []; % normalized root mean squre error for the training dataset
		Correlation_by_gene = []; % between predicted values vs. actual values of every TF gene across all samples
		SigEdges = []; % edges with pValue smaller than cutoff;
	end

	methods (Static)
		function obj=explicit_tfxtf (x, tf_name) % x, tf matrix. sample in rows, tf gene in columns.
			
			cut_off = 0.00001;

			TF = [ones(size(x,1),1) x];
			p = size(TF, 2);
			tr = ones(p,1);
			tr = tr == 1;
			
			I = eye(size(TF, 1));
			b = zeros(size(TF,2), size(TF,2) - 1);
			bp = ones(size(TF,2), size(TF,2) - 1);

			time_trend = zeros(10,1);
			fprintf('Calculating predictor models for %d TFs.\n', p - 1 );

			for i = 2 : p
				
				loop_start_t = clock;

				idx = tr;
				idx(i) = false;
				A = TF(:,idx);
				B = TF(:,i);

				AtA = A' * A;

				beta = AtA \ (A' * B);

				invAtA = inv( AtA );

				H = A * invAtA * A';

				SSE = diag( B' * ( I - H ) * B );

				dof_sse = size(A, 1) - size(A, 2);
				se_beta = sqrt( diag( invAtA ) * (SSE/dof_sse)');
				tStat = beta ./ se_beta;
				beta_pvalue = tcdf( -abs(tStat), dof_sse) * 2;

				b(idx,i - 1) = beta;
				bp(idx, i - 1) = beta_pvalue;

				loop_time = etime(clock, loop_start_t);
				idx_time  = mod(i - 2, 10) + 1;
				time_trend(idx_time) = loop_time;
				average_loop_time = mean(time_trend);
				time_left = (p - i) * average_loop_time / 3600 ;


				if i == 11
					fprintf('Estimated to complete in %.2f hours.\n', time_left);
				end

				if mod(i-1,50) == 0 & i < p
					fprintf('%d TFs done. %.2f hours to go.\n', i - 1, time_left);
				end

				if i == p
					fprintf('%d TFs done.\n', i - 1);
				end

			end

			beta = b;
			beta_pvalue = bp;
			obj.beta = b;
			obj.beta_pvalue =bp;

			TFp = TF * beta;
			TFo = TF(:,2:p);
			cc = arrayfun(@(k) corr(TFp(:,k),TFo(:,k)),1:size(TFp,2),'Uni',1);
			obj.Correlation_by_gene = cc';

			r = TFp - TFo;
			obj.NRMSE = sqrt(sum(sum(r.^2)) / sum(sum(TFo.^2)));

			idx = find(beta_pvalue <= cut_off );
			e3 = beta(idx);
			e3 = round(e3,4);
			e4 = beta_pvalue(idx);
			e4 = round(e4,4,'significant');
			[e1, e2] = ind2sub(size(beta), idx);

			if nargin > 1
				tf_name = string(tf_name);
				gene_name = tf_name;
				tf_name = ["intercept" tf_name'];
				tf_name = regexp( tf_name,'([a-zA-Z0-9_\.]+)','once','match');
				gene_name = string(gene_name);
				gene_name = gene_name';
				gene_name = regexp( gene_name,'([a-zA-Z0-9_\.]+)','once','match');
				obj.TF_name = gene_name';
				e1 = tf_name(e1);
				e2 = gene_name(e2);
				e1 = e1';
				e2 = e2';
			else
				tf_name = 1:size(TF,2);
			end

			colName = {'Gene';'TF';'beta';'beta_pvalue'};
			SigEdges = table(e2,e1,e3,e4,'VariableNames',colName);
			obj.SigEdges = SigEdges;
			fprintf('\n\n');

		end
	end
end

			
