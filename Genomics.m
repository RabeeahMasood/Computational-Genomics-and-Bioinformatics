clear all
close all
fprintf('..........Question 1..........')

% Getting the sequence from Genbank and displaying it

walrus_seq = getgenbank('NC_004029', 'SEQUENCEONLY', true);
walrus = getgenbank('NC_004029')
len = length(walrus_seq)

ntdensity(walrus_seq)

bases = basecount(walrus_seq)

figure
basecount(walrus_seq, 'chart', 'pie')
title('Distribution of Nurcleotide bases for Walrus mitochondrion')

%Finding all ORFs (Protien Coding regions)
ORFs = seqshoworfs(walrus,'GeneticCode','Vertebrate Mitochondrial',...
    'frames','all','minimumlength',75)

%Finding and translating CYTB and CYTC(COX2) from all available genes 
Features = featureparse(walrus,'Sequence',true);
Codingsequence =  Features.CDS;
codingsequence_ID = sprintf('%s',Codingsequence.gene)

%Translating CYTB
CYTB_ID = Codingsequence(13);
gCYTB = walrus.Sequence(14173:15312)
CYTB = nt2aa(gCYTB,'GeneticCode',2)

%Translating (CYTC)COX2
CYTC_ID = Codingsequence(4)
gCYTC = walrus.Sequence(7034:7718)
CYTC = nt2aa(gCYTC,'GeneticCode',2)

fprintf('..............Question 3..........')

        % Species Description      Accession Number 
data = {'Atlantic Walrus'          'NP_659349';       
        'Grey Seal'                'ACZ28998';
        'Wolverine'                'YP_001382271';
        'Harbour Seal'             'BAI60013';
        'Bearded Seal'             'YP_778837';
        'Long-tailed weasel'       'AFP32883';
        };
    
    for ind = 1:length(data)
        primates1(ind).Header = data{ind,1};
        primates1(ind).Sequence = getgenpept(data{ind,2},'sequenceonly','true');
    end
    
    
    
    distances = seqpdist(primates1,'Method','Jukes-Cantor','Alpha','DNA')
    UPGMAtree = seqlinkage(distances,'UPGMA',primates1)
    
    
    h = plot(UPGMAtree,'orient','top');
    title('Rooted Phylogentic tree using Jukes-Cantor model')
    ylabel('Evolutionary Distance')

 fprintf('...........Question 4..........')
%ACCESSION NUMBRS USED FOR INDIVIDUAL PROTIEN SINCE COULD NOT FIND FULL
%GENOME

 %Distance between AA of Atlantic Walrus 
 AW1={'CYTB'   'NP_659349';
      'CYTC'   'NP_659340.3';};
  
    for ind = 1:length(AW1)
        AWal(ind).Header = AW1{ind,1};
        AWal(ind).Sequence = getgenpept(AW1{ind,2},'sequenceonly','true');
    end
   
  Aw1Dis = seqpdist(AWal,'method','jukes-cantor','alpha','DNA')
  
  %ditance between NT Atlantic Walrus
  
 AW2={'CYTB'   'X82299';
      'CYTC'   'AY377171';};
  
    for ind = 1:length(AW2)
        AWal1(ind).Header = AW2{ind,1};
        AWal1(ind).Sequence = getgenbank(AW2{ind,2},'sequenceonly','true');
    end
   
  Aw2Dis = seqpdist(AWal1,'method','jukes-cantor','alpha','DNA')
  
  
 %distance between AA of Grey seal
   
 GS1={'CYTB'   'ACZ28998';
      'CYTC'   'NP_007072';};
  
    for ind = 1:length(GS1)
        GSeal(ind).Header = GS1{ind,1};
        GSeal(ind).Sequence = getgenpept(GS1{ind,2},'sequenceonly','true');
    end
   
  GS1Dis = seqpdist(GSeal,'method','jukes-cantor','alpha','DNA')

 %distance between NT of grey seal 
   
 GS2={'CYTB'   'GU733679';
      'CYTC'   'GU733706';};
  
    for ind = 1:length(GS2)
        Gseal1(ind).Header = GS2{ind,1};
        Gseal1(ind).Sequence = getgenbank(GS2{ind,2},'sequenceonly','true');
    end
   
  GS2Dis = seqpdist(Gseal1,'method','jukes-cantor','alpha','DNA')

 %distance between AA of wolverine 
   
 wr={'CYTB'   'YP_001382271';
      'CYTC'  'YP_001382262';};
  
    for ind = 1:length(wr)
        Wol(ind).Header = wr{ind,1};
        Wol(ind).Sequence = getgenpept(wr{ind,2},'sequenceonly','true');
    end
   
  wrDis = seqpdist(Wol,'method','jukes-cantor','alpha','DNA')

 %distance between NT of wolverine 
   
 wr2={'CYTB'   'DQ206375';
      'CYTC'   'AY377174';};
  
    for ind = 1:length(wr2)
        Wol1(ind).Header = wr2{ind,1};
        Wol1(ind).Sequence = getgenbank(wr2{ind,2},'sequenceonly','true');
    end
   
  wr2Dis = seqpdist(Wol1,'method','jukes-cantor','alpha','DNA')

 %distance between AA of Dog 
   
dog1={'CYTB'   'ABY61049';
      'CYTC'   'NP_001183974';};
  
    for ind = 1:length(dog1)
        Dog1(ind).Header = dog1{ind,1};
        Dog1(ind).Sequence = getgenpept(dog1{ind,2},'sequenceonly','true');
    end
   
  Dog1Dis = seqpdist(Dog1,'method','jukes-cantor','alpha','DNA')

 %distance between NT of Dog
   
 dog2={'CYTB'   'EU352854';
       'CYTC'   'NM_001197045.1';};
  
    for ind = 1:length(dog2)
        Dog2(ind).Header = dog2{ind,1};
        Dog2(ind).Sequence = getgenbank(dog2{ind,2},'sequenceonly','true');
    end
   
  dogDis = seqpdist(Dog2,'method','jukes-cantor','alpha','DNA')

 %distance between AA of Spectacled Bear
   
  SB={'CYTB'   'AAB50570';
      'CYTC'   'YP_001542732';};
  
    for ind = 1:length(SB)
        Sbear(ind).Header = SB{ind,1};
        Sbear(ind).Sequence = getgenpept(SB{ind,2},'sequenceonly','true');
    end
   
  SBDis = seqpdist(Sbear,'method','jukes-cantor','alpha','DNA')

 %cannot find distance between NT of Spectacled Bear due to no CYTC blast

 %distance between AA of SA sea Lion
   
 sl={'CYTB'   'AAQ95107';
      'CYTC'  'AAR00312';};
  
    for ind = 1:length(sl)
        sl1(ind).Header = sl{ind,1};
        sl1(ind).Sequence = getgenpept(sl{ind,2},'sequenceonly','true');
    end
   
  sl1Dis = seqpdist(sl1,'method','jukes-cantor','alpha','DNA')

 %distance between NT of SA sea Lion 
   
 SLION={'CYTB'   'AY377328.1';
        'CYTC'   'AY377172.1';};
  
    for ind = 1:length(SLION)
        SL1(ind).Header = SLION{ind,1};
        SL1(ind).Sequence = getgenbank(SLION{ind,2},'sequenceonly','true');
    end
   
  sl2Dis = seqpdist(SL1,'method','jukes-cantor','alpha','DNA')

 
 
%Using the top two-three results along with the orignal walrus we will compute
%the distance between species for each protien
%Not able to find complete genomes for most of the selected species, using
%direct accession numbers instead

%CYTB AA results


    % Species Description      Accession Number 
Data = {'Atlantic Walrus'          'NP_659349';       
        'Grey Seal'                'ACZ28998';
        'Wolverine'                'YP_001382271';
        'Dog'                      'ABY61049';
       % 'Spectacled Bear'          'AAB50570';
        'South American Sea Lion'  'AAQ95107';};
    
    for ind = 1:length(Data)
        primates2(ind).Header = Data{ind,1};
        primates2(ind).Sequence = getgenpept(Data{ind,2},'sequenceonly','true');
    end
    
    data2= {'Atlantic Walrus';
            'Grey Seal';
            'Wolverine';
            'Dog';
          %  'Spectacled Bear';
            'South American Sea Lion';}
          
    
    Distances = seqpdist(primates2,'Method','Jukes-Cantor','Alpha','DNA') 
    Y =  squareform(Distances)
    Coord = cmdscale(Y)
    plot(Coord(:,1),Coord(:,2),'X')
    text(Coord(:,1)+0.005,Coord(:,2),data2)
    title('Multi-demsional Scaling CYTB')
    
    
    %CYTC AA results
    
         %Species Description    Accession Number
 DATA ={'Atlantic Walrus'          'NP_659340.3';       
        'Grey Seal'                'NP_007072';
        'Wolverine'                'YP_001382262';
        'Dog'                      'NP_001183974';
       % 'Spectacled Bear'          'YP_001542732';
        'South American Sea Lion'  'AAR00312';};
    
    for ind = 1:length(DATA)
        primates3(ind).Header = DATA{ind,1};
        primates3(ind).Sequence = getgenpept(DATA{ind,2},'sequenceonly','true');
    end
    
    DATA2= {'Atlantic Walrus';
            'Grey Seal';
            'Wolverine';
            'Dog';
           % 'Spectacled Bear';
            'South American Sea Lion';}
          
    
    Distance = seqpdist(primates3,'Method','Jukes-Cantor','Alpha','DNA')  
    W =  squareform(Distance)
    Coord2 = cmdscale(W)
    plot(Coord2(:,1),Coord2(:,2),'X')
    text(Coord2(:,1)+0.001,Coord2(:,2),DATA2)
    title('Multi-demsional Scaling CYTC')
    
    fprintf('...........Question 5.........')
    
    %CYTB AA Results
    
      AAtree = seqlinkage(Distances,'UPGMA',primates2)   
    t = plot(AAtree,'orient','top');
    title('Rooted Phylogentic tree for CYTB AA')
        ylabel('Evolutionary Distance')
        
   %CYTC AA Results
   
     actree = seqlinkage(Distance,'UPGMA',primates3)  
    u = plot(actree,'orient','top');
    title('Rooted Phylogentic tree for CYTC AA')
    ylabel('Evolutionary Distance')


    
    % CYTB NT Results 
     
        %Species Description    Accession Number
  B  = {'Atlantic Walrus'          'X82299';       
        'Grey Seal'                'GU733679';
        'Wolverine'                'DQ206375';
        'Dog'                      'EU352854';
      %  'Spectacled Bear'          'U23554.1';
        'South American Sea Lion'  'AY377328.1';};
    
    for ind = 1:length(B)
        primates4(ind).Header = B{ind,1};
        primates4(ind).Sequence = getgenbank(B{ind,2},'sequenceonly','true');
    end
    
        
    D = seqpdist(primates4,'Method','alignment-score','Indels','pairwise-delete',...
        'scoringmatrix','pam250','pairwisealignment',true); 
    UPGMA1 = seqlinkage(D,'UPGMA',primates4)
       
    W = plot(UPGMA1,'orient','top');
    ylabel('Evolutionary Distance')
    title('Rooted Phylogenetic tree for CYTB NT')
    
    %Results for CYTC NT
    %Could not find spectacled bear CYTC NT on NCBI 
    
      
       %Species Description    Accession Number
  C  = {'Atlantic Walrus'          'AY377171';       
        'Grey Seal'                'GU733706';
        'Wolverine'                'AY377174';
        'Dog'                      'NM_001197045.1';
        'South American Sea Lion'  'AY377172.1';};
    
    for ind = 1:length(C)
        primates5(ind).Header = C{ind,1};
        primates5(ind).Sequence = getgenbank(C{ind,2},'sequenceonly','true');
    end
          
    
    D1 = seqpdist(primates5,'Method','alignment-score','Indels','pairwise-delete',...
        'scoringmatrix','pam250','pairwisealignment',true); 
    UPGMA2 = seqlinkage(D1,'UPGMA',primates5)
       
    V = plot(UPGMA2,'orient','top');
    ylabel('Evolutionary Distance')
    title('Rooted Phylogenetic tree for CYTC NT')
    
    fprintf('..........Question 7.......')
    
    w = fastaread('Document.fasta.txt')
    dist = seqpdist(w,'scoringmatrix','GONNET')
    tree = seqlinkage(dist,'average',w)
    m = multialign(w,tree,'scoringmatrix',{'pam150','pam200','pam250'})
    showalignment(m)