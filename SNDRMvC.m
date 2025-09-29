function [LABEL] = SNDRMvC(X, K , options, gnd,eta,etaH,beta)
rng('default')
rng(5489)

% ��ʼ������
C = options.nClass; % �����Ŀ
% ���õ�������
Round = 0;
maxRound=30;
% ��ȡ��ͼ����
nView = numel(X);
lumdav = 1/nView;

obj_NMF1 = cell(1, nView);
obj_NMF2 = cell(1, nView);
obj_NMF3 = cell(1, nView);
obj_NMF4 = cell(1, nView);
obj_NMF = cell(1, nView);
Ux = cell(1, nView);
dX=cell(1, nView);


% ��ʼ�� H Ϊ��Ԫ�����飬ÿ����Ԫ��洢һ����ͼ�ĵ�ά��ʾ
H = cell(1, nView);
% ��ʼ�� m Ϊ��Ԫ�����飬ÿ����Ԫ��洢һ����ͼ��ƫ����
m = cell(1, nView);
% ��ʼ�� h �� m
for viewIdx = 1:nView
    [mFea, nSmp] = size(X{viewIdx});
    H{viewIdx} = logsig(rand(mFea, K)); % ��ʼ�� h
    m{viewIdx} = logsig(rand(mFea, 1)); % ��ʼ�� m 
end
% centroidH = logsig(rand(mFea, C));
centroidH=H{1};


% ��ʼ�� ACC �� NMI �������
acc = zeros(1, 30);
NMI = zeros(1, 30);
LABEL = cell(1, 30); % �洢ÿ�ֵı�ǩ

% �������� h �� m
while Round < maxRound
    Round = Round + 1;
   
    % ѭ������ÿ����ͼ
    for viewIdx = 1:nView

        %ȡ����ǰ��ͼ
        current_X=X{viewIdx};
        current_H=H{viewIdx};
        current_m=m{viewIdx};

        % ��ȡ��ǰ��ͼ�ĳߴ�
        [mFea, nSmp] = size(current_X);

        % ������Լ�� SNMF (NNMF)
        for i = 1:mFea
            for j = 1:nSmp
                sum1 = 0;

                for kk = 1:K
                    sum1 = sum1 + 1/(1 + exp(-current_H(i, kk))) * 1/(1 + exp(-current_H(j, kk)));
                end
                temp = current_X(i, j) - sum1 - 1/(1 + exp(-current_m(i)));

                % ����ƫ���� m
                C1 = 1/(1 + exp(-current_m(i)));
                current_m(i) = current_m(i) - eta * (temp * (-1) * C1 * (1 - C1));

                % ���µ�ά���� h ��Ԫ��
                for kk = 1:K
                    A = 1/(1 + exp(-current_H(i, kk)));
                    CH = 1/(1 + exp(-centroidH(i, kk)));
                    current_H(i, kk) = current_H(i, kk) - eta * (temp * (-1)*1/(1+exp(-current_H(j, kk))) + lumdav * (A - CH) + beta * A )* (1 - A) * A;
                end
            end
        end

       
         % �����ά���� h ��ֵ
        H{viewIdx} = current_H;
        % ����ƫ�� m ��ֵ
        m{viewIdx} = current_m;
        
    end
 % ���㹲ʶ�����ֵ
 for viewIdx = 1:nView
        %ȡ����ǰ��ͼ21
        current_H=H{viewIdx}; 
        for i = 1:mFea
            for kk = 1:K
                A = 1/(1 + exp(-current_H(i, kk)));
                CA = 1/(1 + exp(-centroidH(i, kk)));
                centroidH(i, kk) = centroidH(i, kk) - etaH * lumdav * ((-1) *(A-CA))*CA*(1-CA);
            end
        end
 end

 for viewIdx = 1:nView
    H{viewIdx} = logsig(H{viewIdx});
    m{viewIdx} = logsig(m{viewIdx});
 end
 centroidH = logsig(centroidH);

    obj_NMFall(Round) = 0;
    O1=0;
    O2=0;
    O3=0;
 for viewIdx = 1:nView
        current_H=H{viewIdx};
        current_m=m{viewIdx};
        Ux{viewIdx} = [current_H, current_m];
        dX{viewIdx} = Ux{viewIdx}*(Ux{viewIdx})'-X{viewIdx};
        obj_NMF1{viewIdx}(Round) = sum(sum(dX{viewIdx}.^2));
        obj_NMF2{viewIdx}(Round) = beta*sum(sum(current_H.^2));
        obj_NMF3{viewIdx}(Round) = lumdav*sum(sum(current_H-centroidH).^2);
        obj_NMF4{viewIdx}(Round)=obj_NMF1{viewIdx}(Round)+obj_NMF2{viewIdx}(Round)+obj_NMF3{viewIdx}(Round);
        obj_NMF{viewIdx}(Round)=sqrt(obj_NMF4{viewIdx}(Round));
        obj_NMFall(Round) = obj_NMFall(Round)+obj_NMF{viewIdx}(Round);
        O1=O1+obj_NMF1{viewIdx}(Round);
        O2=O2+obj_NMF2{viewIdx}(Round);
        O3=O3+obj_NMF3{viewIdx}(Round);  
 end
  disp(['obj_NMF1��',  num2str(O1)]);
    disp(['obj_NMF2��',   num2str(O2)]);
    disp(['obj_NMF3: ',   num2str(O3)]);
 
    % ��þ����ǩ
    rng(5489)
    label = litekmeans(centroidH, C, 'Replicates', 20);

 % ȷ����ǩ������������
    gnd   = gnd(:);
    label = label(:);

    % �������ӳ�䣨���� ACC��
    idx22 = bestMap(gnd, label);

    % ���� ACC
    LABEL{Round} = idx22;
    acc(Round) = sum(gnd == idx22) / length(gnd);

    % ---- ���� NMI ǰ����ǩӳ�䵽 1..K ----
    [~, ~, gnd_tmp]   = unique(gnd,   'stable');   % �� gnd �ȳ�
    [~, ~, label_tmp] = unique(label, 'stable');   % �� label �ȳ�

    % ����쳣
    if length(gnd_tmp) ~= length(label_tmp)
        error('��ǩ���Ȳ�һ�£�gnd=%d, label=%d', length(gnd_tmp), length(label_tmp));
    end
    if any(isnan(gnd_tmp)) || any(isnan(label_tmp))
        error('��ǩ�к� NaN���������ݡ�');
    end

 % ���� NMI
    NMI(Round) = MutualInfo(gnd_tmp, label_tmp);
    disp(['Number of Rounds:', num2str(Round)]);
    disp(['Object Function Value: ',num2str(obj_NMFall(Round))]);
    disp(['Acc:', num2str(acc(Round))]);
    disp(['NMI:', num2str(NMI(Round))]);
    disp('-------------------------------------------------------');
end
end


