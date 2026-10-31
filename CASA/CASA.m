classdef CASA < ALGORITHM
% <multi/many> <real> <expensive> <constrained>
% Kriging-assisted evolutionary algorithm with two search modes
% lambda --- 0.7 --- Parameter determining the evolving stages 
    
    methods
        function main(Algorithm,Problem)
            %% Parameter setting on KTA2
            phi1  = 0.1;
            wmax1 = 10;
            mu1   = 5;
            %% Parameter setting on KTS
            gap2=10;
            [lambda] = Algorithm.ParameterSet(0.7);
            %% Initialization
            p      = 1/Problem.M;
            CAsize = Problem.N;
            N      = Problem.N;
            P      = UniformPoint(N, Problem.D, 'Latin');
            Population = Problem.Evaluation(repmat(Problem.upper-Problem.lower,N,1).*P+repmat(Problem.lower,N,1));
            A1     = Population;
            UA1    = A1;

            Iter = 1;
            IdealPoints = [];
            search_mode = 0;

            count_case0 = 0;
            count_case1 = 0;
            count_case2 = 0;
            % CA DA in KTA2
            CA = UpdateCA([],Population,CAsize);
            DA1 = UpdateDA(Population,[],Problem.N,p);
            
            % P1 P2 in KCCMO
            P1 = Population;
            P2 = Population;
            DA = DA1;

           
            Dec = [DA.decs;CA.decs;P2.decs];
            Obj = [DA.cons;CA.cons;P2.cons];
            [~,index] = unique( round(Dec*1e4)/1e4,'rows');
            Dec   =  Dec(index,:);
            Obj   = Obj(index,:);
            THETA = 5.*ones((Problem.M + size(DA.cons,2)+1),Problem.D);
            Model = cell(1,(Problem.M + size(DA.cons,2)+1));
           
            for i = Problem.M+1: (Problem.M + size(DA.cons,2))
                j = i - Problem.M;
                dmodel     = dacefit(Dec,Obj(:,j),'regpoly0','corrgauss',THETA(i,:),1e-5.*ones(1,Problem.D),100.*ones(1,Problem.D));
                Model{i}   = dmodel;
                THETA(i,:) = dmodel.theta;
            end
           
            while Algorithm.NotTerminated(UA1)
                IdealPoints(Iter,:) = min(A1.objs,[],1);
                if Iter > gap2 && search_mode == 0
                    max_change = calc_maxchange(IdealPoints,Iter,gap2);
                    MaxCon = max(max(0,A1.cons));
                    Ninf   = length(find(MaxCon>0));
                    if max_change <= lambda && Ninf == size(A1.cons,2)
                        search_mode = 1;
                    elseif max_change <= lambda && Ninf < size(A1.cons,2)
                        search_mode = 2;
                    elseif max_change <= lambda && Ninf  == 0
                        search_mode = 0;
                    end
                end
                if search_mode == 0 
                    DA = DA1;
                    count_case0 = count_case0 + 1;
                elseif search_mode == 1
                    DA = P1;
                    count_case1 = count_case1 + 1;
                elseif search_mode == 2
                    DA = P1;
                    count_case2 = count_case2 + 1;                    
                end
 
                Dec = UA1.decs;
                Obj = UA1.objs;
                for i = 1:(Problem.M )%+ size(DA.cons,2))
                    dmodel     = dacefit(Dec,Obj(:,i),'regpoly0','corrgauss',THETA(i,:),1e-5.*ones(1,Problem.D),100.*ones(1,Problem.D));
                    Model{i}   = dmodel;
                    THETA(i,:) = dmodel.theta;
                end

                if search_mode == 1
                   
                    Dec = [CA.decs;P2.decs];
                    Obj = [CA.cons;P2.cons];
                    [~,index] = unique( roundn(Dec,-4) ,'rows');
                    Dec =  Dec(index,:);
                    Obj =  Obj(index,:);
                    for i = Problem.M+1: (Problem.M + size(DA.cons,2))
                        j = i - Problem.M;
                        dmodel     = dacefit(Dec,Obj(:,j),'regpoly0','corrgauss',THETA(i,:),1e-5.*ones(1,Problem.D),100.*ones(1,Problem.D));
                        Model{i}   = dmodel;
                        THETA(i,:) = dmodel.theta;
                    end
                end

              if search_mode == 2 
              Dec = [DA.decs; P2.decs]; 
              Obj = [DA.cons; P2.cons];  
              [~,index] = unique( roundn(Dec,-4) ,'rows');
               Dec =  Dec(index,:);
               Obj =  Obj(index,:);
                    for i = Problem.M+1: (Problem.M + size(DA.cons,2))
                        j = i - Problem.M;
                        dmodel     = dacefit(Dec,Obj(:,j),'regpoly0','corrgauss',THETA(i,:),1e-5.*ones(1,Problem.D),100.*ones(1,Problem.D));
                        Model{i}   = dmodel;
                        THETA(i,:) = dmodel.theta;
                    end
            end
            
                % Set the CCA and CDA as the current CA and DA
                if search_mode == 0
                    CCA.obj = CA.objs; CCA.dec = CA.decs; CCA.con = CA.cons;  CCA.MSE = zeros(size(CCA.con,1),Problem.M+size(CCA.con,2));
                    CP2.obj = P2.objs; CP2.dec = P2.decs; CP2.con = P2.cons; CP2.MSE = zeros(size(CP2.con,1),Problem.M+size(CP2.con,2));
                    CDA.obj = DA.objs; CDA.dec = DA.decs; CDA.con = DA.cons; CDA.MSE = zeros(size(CDA.con,1),Problem.M+size(CDA.con,2));
                elseif search_mode == 1
                    CCA.obj = CA.objs; CCA.dec = CA.decs; CCA.con = CA.cons;  CCA.MSE = zeros(size(CCA.con,1),Problem.M+size(CCA.con,2));
                    CP2.obj = P2.objs; CP2.dec = P2.decs; CP2.con = P2.cons; CP2.MSE = zeros(size(CP2.con,1),Problem.M+size(CP2.con,2));
                    CDA.obj = DA.objs; CDA.dec = DA.decs; CDA.con = DA.cons; CDA.MSE = zeros(size(CDA.con,1),Problem.M+size(CDA.con,2));
                elseif search_mode == 2
    
    CCA.obj = CA.objs;
    CCA.dec = CA.decs;
    CCA.con = CA.cons;
    CCA.MSE = zeros(size(CCA.con,1), Problem.M + size(CCA.con,2));
  
    DAConMax = max(DA.cons, [], 2);   
    P2ConMax = max(P2.cons, [], 2);

    idxDA_inf = find(DAConMax > 0);    
    idxP2_inf = find(P2ConMax > 0);  

    if isempty(idxDA_inf)
        idxDA_inf = 1:length(DA);
    end
    if isempty(idxP2_inf)
        idxP2_inf = 1:length(P2);
    end

    DA_inf = DA(idxDA_inf);
    CDA.obj = DA_inf.objs;
    CDA.dec = DA_inf.decs;
    CDA.con = DA_inf.cons;
    CDA.MSE = zeros(size(CDA.con,1), Problem.M + size(CDA.con,2));

    P2_inf = P2(idxP2_inf);
    CP2.obj = P2_inf.objs;
    CP2.dec = P2_inf.decs;
    CP2.con = P2_inf.cons;
    CP2.MSE = zeros(size(CP2.con,1), Problem.M + size(CP2.con,2));
end

                w = 1;
                while w <= wmax1
                    if search_mode == 0
                        [~,ParentCdec,~,ParentMdec] = MatingSelection_KTA2(CCA.obj,CCA.dec,CDA.obj,CDA.dec,Problem.N);
                        OffspringDec = [OperatorGA(Problem,ParentCdec,{1,20,0,0});OperatorGA(Problem,ParentMdec,{0,0,1,20})];
                    else
                        Fitness1 = CalFitness(CP2.obj,CP2.con);
                        Fitness2 = CalFitness(CDA.obj);
                        
                        MatingPool1 = TournamentSelection(2,Problem.N,Fitness1);
                        MatingPool2 = TournamentSelection(2,Problem.N,Fitness2);
                        
                        Offspring1  = OperatorGA(Problem,CP2.dec(MatingPool1,:));
                        Offspring2  = OperatorGA(Problem,CDA.dec(MatingPool2,:));
                        OffspringDec = [Offspring1;Offspring2];
                    end
                    
                   
                        Pop.dec = OffspringDec;
                        N       = size(Pop.dec,1);
                        Pop.obj = zeros(N,Problem.M);
                        Pop.con = zeros(N,size(DA.cons,2));
                        Pop.MSE = zeros(N,Problem.M+ size(DA.cons,2));
                        PopObj  = zeros(N,Problem.M+ size(DA.cons,2));
                    
                        for i = 1 : N
                            for j = 1 : (Problem.M+size(DA.cons,2))
                                [PopObj(i,j),~,Pop.MSE(i,j)] = predictor(Pop.dec(i,:),Model{j});
                            end
                        end
                    
                        Pop.obj = PopObj(:,1:Problem.M );
                        Pop.con = PopObj(:,Problem.M+1 :end );
                    
                    Pop = orderfields(Pop, {'obj', 'dec', 'con', 'MSE'});

                    
                        Pop.MSE = Pop.MSE(:, 1:Problem.M+size(DA.cons,2));
                        PopC = cat_struct(CCA,Pop);
                        CCA  = K_UpdateCA(PopC,CAsize);
                        PopD = cat_struct(CDA,Pop);
                   

                    if search_mode == 0
                        CDA = K_UpdateDA(PopD,Problem.N,p);
                        PopC1   = cat_struct(CP2,Pop);
                        [CP2,~] = K_UpdateP(PopC1,Problem.N,true);
                    elseif search_mode == 1
                        if size(CP2.obj, 2) ~= size(Pop.obj, 2)
                            Pop.obj = Pop.obj(:, 1:Problem.M);
                        end
                        if size(CP2.MSE, 2) ~= size(Pop.MSE, 2)
                           CP2.MSE = CP2.MSE(:, 1:Problem.M+size(DA.cons,2));
                           Pop.MSE = Pop.MSE(:, 1:Problem.M+size(DA.cons,2));
                        end

                     
                        CDA = K_UpdateP(PopD, Problem.N, false);
                        PopC1 = cat_struct(CP2, Pop);
                        [CP2, ~] = K_UpdateP(PopC1, Problem.N, true);
                    elseif search_mode == 2
                        CDA = K_UpdateP(PopD,Problem.N,false);
                        PopC1   = cat_struct(CP2,Pop);
                        [CP2,~] = K_UpdateP(PopC1,Problem.N,true);
                    end

                    w = w + 1;
                end
                
                if search_mode == 0
                    [~,ia,~] = setxor(CCA.dec,A1.decs,'rows');
                    CCA = givevalue(CCA,ia);
                    [~,ia,~] = setxor(CDA.dec,A1.decs,'rows');
                    CDA = givevalue(CDA,ia);
                    
                    Offspring01 = Adaptive_sampling(CCA.obj,CDA.obj,CCA.dec,CDA.dec,CDA.MSE,DA,P2,mu1,p,phi1);
                else
                    % KCCMO
                    [CCA2,~]    = KCCMO_sampling(CP2,P2,mu1);
                    Offspring01 = CCA2.dec;

                end
                [Offspring01] = unique(Offspring01, 'rows', 'stable');
                [~,index] = unique(Offspring01 ,'rows');
                PopNew = Offspring01(index,:);
                Offspring02 = PopNew;
                
                if ~isempty(Offspring02)  
                    Offspring02_rounded = round(Offspring02, 4);
                    Offspring02 = unique(Offspring02_rounded, 'rows', 'stable');

                    Offspring = Problem.Evaluation(Offspring02);
                    temp =  A1.decs;
                    for i = 1 : size(Offspring,2)
                        dist2 = pdist2(Offspring(i).decs,temp);
                        if min(dist2) > 1e-5
                            A1 = [A1,Offspring(i)];
                            UA1 = Update_P(A1,Problem.N,true);
                        end
                        temp = A1.decs;
                    end
                     
                    CA     = UpdateCA(CA,Offspring,CAsize);
                    DA1    = UpdateDA(DA1,Offspring,Problem.N,p);
                    [P1,~] = Update_P([P1,Offspring],Problem.N,false);
                    [P2,~] = Update_P([P2,Offspring],Problem.N,true);
                    if Problem.FE >= Problem.maxFE
                       fprintf('Case 0 was used %d times.\n', count_case0);
                       fprintf('Case 1 was used %d times.\n', count_case1);
                       fprintf('Case 2 was used %d times.\n', count_case2);
                    end
                end
                Iter = Iter + 1;
            end
        end
    end
end
%% 


function CV = NormalizeCV(PopCon)
    % Calculate the normalized overall constraints violation
    PopCon = max(0,PopCon);
    PopCon = (PopCon-min(PopCon,[],1))./(max(PopCon,[],1)-min(PopCon,[],1));
    PopCon(:,isnan(PopCon(1,:))) = 0;
    CV = sum(PopCon,2);
end

function max_change = calc_maxchange(ideal_points,Iter,gap)
    % Calculate the maximum change rate of ideal points
    delta = 1e-6 * ones(1,size(ideal_points,2));
    rz    = abs((ideal_points(Iter,:) - ideal_points(Iter - gap,:)) ./ max(ideal_points(Iter - gap,:),delta));  
    max_change = max(rz);
end