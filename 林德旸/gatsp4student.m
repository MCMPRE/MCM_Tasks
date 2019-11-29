function gatsp()
	hold off
	nc = 20;		% ���еĸ���
	N = 100;		% Ⱦɫ��Ⱥ���ģ
	MAXITER = 10^3;	% ���ѭ������
	SelfAdj = 0; 	%SelfAdj = 1ʱΪ����Ӧ
    
	pc = 0.85;		% �������
	pw = 0.15;		% �������
	
	locationx=floor(rand(nc,1)*100);
	locationy=floor(rand(nc,1)*100);
	locations=[locationx,locationy];
    %load locations10.mat;
    load locations20.mat;
    %load locations100.mat;
	R = distCal(locations);
    
	
	% ����1������N��Ⱦɫ��ĳ�ʼȺ��,������pop����
	pop = initPop(N,nc);
	
	iter=0;
    tstart=clock;
	while iter<MAXITER
        iter = iter+1;
		trajLength=calLen(pop,R);				%����2������ÿ��Ⱦɫ���·������		
		fitness=calFitness(trajLength);			%      ����ÿ��Ⱦɫ�����Ӧֵ	%������Ҫ��д
        [b,f,ind]=findBest(pop,fitness);
		if satified(fitness), break; end		%����3���������ĳЩ��׼���㷨ֹͣ		
		pop = chooseNewP(pop,fitness);			%      ����ѡȡ��һ���µ�Ⱥ��	%������Ҫ��д
		pop = crossPop(pop,pc,fitness,SelfAdj);	%����4�����������Ⱦɫ�壬�õ���Ⱥ�� %������Ҫ��д		
		pop = mutPop(pop,pw,fitness,SelfAdj); 	%����5���������	%������Ҫ��д
		pop(1,:) = b(1,:);						% ������һ������Ӧֵ��ߵ�Ⱦɫ��             
    end	    
	
	%�������Ⱦɫ��/·��
	disp(strcat('���ŵ�·������Ϊ��',num2str(calLen(b(1,:),R)),'����ӦֵΪ��',num2str(f),'����Ӧ��·��Ϊ��'));
	disp(b(1,:));			
    
	%��������ʱ��
    tend=clock;
	disp(strcat('The program need : ',num2str(etime(tend,tstart)),' seconds.'));
	
	%��ʾ����·��ͼ��
	drawTSP(locations, b(1,:));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ���ݸ����е�location�������������
% ��ľ������R
function R = distCal(location)
	nc=size(location,1);
	R = zeros(nc);
	for i=1:nc
		for j=i+1:nc
			R(i,j) = sqrt(sum((location(i,:)-location(j,:)).^2));
		end
	end
	R = R+R';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ����N��Ⱦɫ��ĳ�ʼȺ��,������pop��
% ÿ��Ⱦɫ�����TSP�����ĳһ���⣨�����г��ж�����һ�εĹ켣��
function pop=initPop(N,nc)
	pop = zeros(N,nc);
	for i=1:N,pop(i,2:end) = randperm(nc-1)+1;end %ʹ�������������N��Ⱦɫ��
	pop(:,1)=ones(N,1);  %����Ⱦɫ�嶼�ӳ���1��ʼ�����ص�����1.
end	



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%���ݾ������R ������pop��ÿ��Ⱦɫ��������Ĺ켣�ĳ���
function trajLength = calLen(pop,R)
	[N,nc]=size(pop);
    trajLength = zeros(1,N);
	for i=1:N
		f=0;
		traj = pop(i,:);
		traj_2 = [traj(2:end),traj(1)];
		for j=1:nc , f = f+R(traj(j),traj_2(j)); end
		trajLength(i) = f;
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ��������Ⱦɫ��ֱ�Ĺ켣���ȣ��������Ⱦɫ�����Ӧֵ
% trajLength	: ����Ⱦɫ��Ĺ켣����
% fitness		������Ⱦɫ�����Ӧֵ
function fitness=calFitness(trajLength)  %������Ҫ��д
	%fitness = ones(size(trajLength)); 
    fitness = 1./trajLength; %trajLength��i��Ϊ������
    Cmult=1.6;
    favg = mean(fitness(:));
    fmin = min(fitness);
    fmax = max(fitness);
   % a = (Cmult-1)*favg/(fmin-favg);
   % b = (fmax-Cmult*favg)*favg/(fmax-favg);
    a = favg / (favg - fmin);
    b = -fmin*favg / (favg - fmin);
    fitness = a*fitness + b; 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ����Ⱦɫ�����Ӧֵ������һ���ĸ��ʣ�������һ��Ⱦɫ��Ⱥ��
% pop		������Ⱦɫ��Ⱥ�壬 N x nc�ľ���
% fitness	����Ⱦɫ�����Ӧֵ������ΪN������
% newpop	�������µ�Ⱦɫ��Ⱥ�壬N x nc�ľ���
function newpop = chooseNewP(pop,fitness)	%������Ҫ��д
    Ps = fitness ./ sum(fitness); %�����ѡ�����
    newPs=sort(Ps,'descend'); %��Ps�������򣬱������ۼƸ���
    [Ps,index] = sort(Ps); %�ҵ�����ǰ������
    newpop = zeros(100,10);
    %�ۼƸ���sumPs
    sumPs = newPs;  %sumPsΪ�ۼƸ���
    for i=2:100
        sumPs(i) = sumPs(i)+sumPs(i-1);
    end
    b=1;
    for i = 1:100
        p = rand(1,1); 
        for j = 1:100
            if(sumPs(j) >= p) %��ѡ����
               out = index(j);
               newpop(b,:) = pop(out,:);  %�ҵ�����ѡ�ĸ��壬��ֵ���¸���Ⱥ��
               b=b+1;
               break
            end   
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ���ݽ������pc���Լ���Ⱦɫ�����Ӧֵfitness��ͨ������ķ�ʽ������Ⱥ��
% SelfAdj = 1ʱΪ����Ӧ������ȡ�̶��Ľ������pc
% pop		������Ⱦɫ��Ⱥ�壬 N x nc�ľ���
% fitness	����Ⱦɫ�����Ӧֵ������ΪN������
% cpop		�������µ�Ⱦɫ��Ⱥ�壬N x nc�ľ���
% pc		: �������
function cpop = crossPop(pop,pc,fitness,SelfAdj) %������Ҫ��д	
	cpop = pop;
    for i=1:100
        k = randi(100);
        j = randi(100);
        p = rand(1,1);
        if( p <=pc )
            [cpop(k,:),cpop(j,:)] = genecross(pop(k,:),pop(j,:));
        end
    end
	%[child1,child2] = genecross(partner1,partner2); %���genecross�������ܻ�����
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ��Ⱦɫ��partner1,partner2��ͨ�����淽ʽ
% ����������Ⱦɫ��child1,child2
% partner1/2	: ����ǰ������Ⱦɫ��
% child1/2		������������Ⱦɫ��
function [child1,child2] = genecross(partner1,partner2)
	len = length(partner1);
	idx1 = randi(len-1)+1;
	idx2 = randi(len-1)+1;
	ind1 = min(idx1,idx2);
	ind2 = max(idx1,idx2);
	
	child1 = partner1;
	child2 = partner2;
	
	tem1 = child1(ind1:ind2);
	tem2 = child2(ind1:ind2);	
	
	temdff1 = setdiff(tem1,tem2);
	temdff2 = setdiff(tem2,tem1);	

	for i=1:length(temdff1)
		child1(find(child1==temdff2(i)))=temdff1(i);
		child2(find(child2==temdff1(i)))=temdff2(i);		
	end	
	
	child1(ind1:ind2) = tem2;
	child2(ind1:ind2) = tem1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ���ݱ������pw���Լ���Ⱦɫ�����Ӧֵfitness��ͨ������ķ�ʽ������Ⱥ��
% %SelfAdj = 1ʱΪ����Ӧ������ȡ�̶��ı������pw
% pop		������Ⱦɫ��Ⱥ�壬 N x nc�ľ���
% fitness	����Ⱦɫ�����Ӧֵ������ΪN������
% mpop		�������µ�Ⱦɫ��Ⱥ�壬N x nc�ľ���
% pw		: �������
function mpop = mutPop(pop,pw,fitness,SelfAdj) %������Ҫ��д
	mpop = pop;
    for i=1:100
        k = randi(10);
        j = randi(10);
        p = rand(1,1);
        if( p <= pw )
          mpop(:,[k,j]) = mpop(:,[j,k]);
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ����Ⱦɫ��Ⱥ��pop�Ѿ���Ӧ����Ӧֵfitness��
% �ҵ���ߵ���Ӧֵf���Լ���Ӧ��Ⱦɫ��bst������pop�еı��/�±�ind
function [bst,f,ind] = findBest(pop,fitness) 
    f = max(fitness);
	ind = find(fitness == f);
	bst = pop(ind,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% �������е�ĳ��״̬������ĳЩ��׼����ֹͣ���С�
% �ڵ�ǰ�����£������þ���ֹͣ��׼�����Է���sat=0
function [sat]=satified(fitness)
	sat=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ����Ⱦɫ��traj���涨�Ĺ켣���Լ������е�λ��locations�������ù켣
function drawTSP(locations, traj)
	nc = size(locations, 1);
	plot(0,0,'.'); hold on; plot(100,100,'.');
	for i = 1:nc
		indP = traj(i);
		strPx = locations(traj(i), 1);
		strPy = locations(traj(i), 2);
		endPx = locations(traj(mod(i,nc)+1), 1);
		endPy = locations(traj(mod(i,nc)+1), 2);		
		plot([strPx,endPx],[strPy,endPy],'ms-','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g');
		text(strPx,strPy,['  ', int2str(indP)]);		
	end
end
 