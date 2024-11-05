clear all;clc;
Data_Mask='mask.nii';
[Data_Mask,header]=y_Read(Data_Mask);
table=unique(Data_Mask);
table(1)=[];
Z=length(table);
Data_Path='PET_data';
File=dir(Data_Path);
File(1:2)=[];
L=length(File);n=2^9;
Maxtric_PDF=zeros(L,n,Z);
for i=1:L
    disp(i)
    Img_data = y_Read([Data_Path , File(i).name ]);
    for e=1:Z
        Data_Mask='mask.nii';
        [Data_Mask,header]=y_Read(Data_Mask);
        j=table(e);
        Region=find(Data_Mask~=j);  Data_Mask(Region)=0;
        region=find(Data_Mask==j);  Data_Mask(region)=1;
        number=find(Data_Mask~=0);
        Result=Img_data(find(Data_Mask));
        Range=max(Result)-min(Result);
        MIN= min (Result)-Range/10; MAX=max(Result)+Range/10;
        [bandwidth,density,xmesh,cdf]=kde(Result,n,MIN,MAX);
        Maxtric_PDF(i,:,e)=density;
    end
end
%%
%%Calculate the KL divergence
Maxtric=zeros(L,Z,Z);
Maxtric_KLS=zeros(L,Z,Z);
JS=zeros(L,Z,Z);
for l=1:L
    for m=1:Z
        for k=1:Z
            vec1=Maxtric_PDF(l,:,m); vec1 = vec1/sum(vec1);
            vec2=Maxtric_PDF(l,:,k); vec2 = vec2/sum(vec2);% Make sure vec1 and vec2 sum to 1
            [score_KL, score_JS] = KL_JS_div(vec1, vec2);
            KLS=exp(-score_KL);JS(l,m,k)=score_JS;
            Maxtric(l,m,k)=score_KL;
            Maxtric_KLS(l,m,k)=KLS;
        end
    end
end
