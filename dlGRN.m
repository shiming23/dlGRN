function recgn = dlGRN(expr_matrix_tf, expr_matrix_tg, latents_num, ...
    samp_percentage, samp_times, alpha, T_max)
% dlGRN(expr_matrix_tf, expr_matrix_tg, latents_num,samp_percentage, samp_times,alpha,T_max)
% expr_matrix_tf:ת¼���ӵı��ˮƽ����
% expr_matrix_tg:Ŀ�����ı��ˮƽ����
% latents_num: ��l��Ҳ����AR����Ŀ��
% samp_percentage:�ز�������ռ������������
% samp_times:�ز���������
if size(expr_matrix_tf,1)~=size(expr_matrix_tg,1)
    error('ת¼������Ŀ��������������ͬ��');
end
samp_num=round(size(expr_matrix_tf,1)*samp_percentage);
recgn=zeros(size(expr_matrix_tf,2), size(expr_matrix_tg,2));%��ʼ����������

for ii=1:samp_times
    iitsample=randperm(size(expr_matrix_tf,1));
    iisample=iitsample(1:samp_num);
    iiexpr_tf=zscore(expr_matrix_tf(iisample,:));
    iiexpr_tg=zscore(expr_matrix_tg(iisample,:));
    
    params.data=iiexpr_tg;
    params.Tdata=T_max;
    params.dictsize=latents_num;
    params.iternum=200;
    params.memusage='high';
    
    [Dksvd, g, ~]=ksvd(params,'');
    logg= (g~=0);
    tf2atom=abs(iiexpr_tf'*Dksvd);
    
    for jj=size(recgn,1)
        for kk=size(recgn,2)
            recgn(jj,kk)=recgn(jj,kk)+quantile(tf2atom(jj,logg(:,kk)), alpha);
        end
    end
end

recgn=recgn./( sqrt(samp_num-1)*samp_times);
end