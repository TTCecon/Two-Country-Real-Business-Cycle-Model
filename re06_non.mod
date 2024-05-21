// Institute of Economics, NSYSU 
// Kuan Lun Chen (2024.4)
// This is a two country RBC model (non-linear adjust b) 
// Replicate Iacoviello & Minetti (2006)"National Business Cycles With Domestic and Foreign Lenders"(JME)
// This replicate of the model is two country RBC model with endogeneous alpha but without investment adjustment cost 
// i代表Household面對的變數，如消費ci...f代表國外，如國外家戶房產數量hf
//bH,bF為國內企業家面臨之借貸限制，其s.s.為bHs,bFs
//bHf,bFf為國外企業家面臨之借貸限制，其s.s.為bHsf,bFsf
//endogeneous alpha case and no adjustment costs

var c ci h hi hf hif y yf cf cif bb bH bF bHf bFf q qf R alp alpf A Af; 
varexo e_Af;

parameters b g mH mF j nu eta tau psi alps alpsf Rs bHsy bFsy bHsfy bFsfy bbs cisy cisfy ys ysf qs qsf qshs qsfhsf bHs bFs bHsf bFsf cis cs csf cisf hs hsf his hisf As Asf rho_A rho_Af tau_A tau_Af;
b           = 0.99; //*beta*
g           = 0.98; //*gamma*
mH          = 0.9;  
mF          = 0.8;
j           = 0.1;
nu          = 0.1;
eta         = 21;
tau         = 0.156;//不確定
rho_A       = 0;
rho_Af      = 0;
psi         = 0;
As           =1;
Asf          =1;
tau_A     =0;
tau_Af    =0;

//steady-state

    alps = 1-((1-mH)/(2*(1-mF))); // alpha
    alpsf = 1-((1-mH)/(2*(1-mF))); //domestic and foreign are symmetric 
    Rs = 1/b; //interest rate between two countries are equal

    bHsy = (mH*alps*g*nu)/(Rs*(1-(b-g)*mH-g));//bHsy means bH over y
    bFsy = (1-alps)*(1-(1-mF)*(1-alps))*(bHsy)/(alps*mH);//bFsy means bF over y
    bHsfy = bHsy;//domestic and foreign are symmetric 
    bFsfy = bFsy;//domestic and foreign are symmetric 
    bbs = 0;//Maybe 0.00001

    cisy = (1-b)*Rs*(bHsy+bFsy)+(1-nu);
    cisfy = cisy;//domestic and foreign are symmetric 
    hs = (1-b)*g*nu/((1-b)*g*nu+j*(1-(b-g)*mH-g)*cisy);
    hsf = hs; //ok
    ys = (hs^nu)*((1-nu)/(tau*cisy))^((1-nu)/(eta)); //ok
    ysf = ys; //ok

    bHs = bHsy*ys;//ok
    bFs = bFsy*ys;//ok
    bHsf = bHs;
    bFsf = bFs;
    qs = g*nu*ys/(hs*(1-g-(b-g)*mH)); // ok bHs/alp*b*mH*hs
    qsf = qs; //ok
    qshs = qs*hs; //ok
    qsfhsf = qshs; //ok
    
    
    cs = nu*ys+(1-Rs)*(bHs+bFs);//ok
    cis = cisy*ys;//ok
    csf = cs;
    cisf = cis;

    his = 1-hs; //ok
    hisf = 1-hsf; //ok

    //his = (j*cis)/((1-b)*qs); //ok
    //hisf = his; //ok
   //hs備選：1-(j*cis)/((1-b)*qs)
   //hsf備選：1-(j*cisf)/((1-b)*qsf)
   //hs備選2:(g*nu*(1-b)*ys)/(j*(1-g-(b-g)*mH*cis))/(1+(g*nu*(1-b)*ys)/(j*(1-g-(b-g)*mH*cis)))

model;

///Home

nu*y+bH+bF=
bH(-1)*R(-1)+bF(-1)*R(-1)+q*(h-h(-1))+c;//廠商的預算限制式

y=A*h(-1)^nu*((1-nu)*y/(tau*ci))^((1-nu)/(eta));//廠商生產函數

bH=(alp*mH*q(+1)*h)/R;//廠商和國內借貸面臨之借貸限制式

bF=(q(+1)*(1-alp)*h*(1-(1-mF)*(1-alp)*q(+1)*h/qshs))/R;//廠商和國外借貸面臨之借貸限制式

q=c*(1/(R*c)-g/c(+1))*mH*q(+1)+c*(g/c(+1))*((nu*y(+1)/h)+q(+1));//廠商最適房地產需求

ci+bH+bFf+bb-q*(h-h(-1))+psi*(bb-bbs)^2/2=
bH(-1)*R(-1)+bFf(-1)*R(-1)+R(-1)*bb(-1)+(1-nu)*y;//家戶的預算限制式，wl已被一階最適汰換ci+bH+bFf+bb+q*(hi-hi(-1))+psi*(bb-bbs)^2/2

(1/ci)*(1+psi*(bb-bbs))=b*(R/ci(+1));//去除投資調整成本(1/ci)*(1+psi*(bb-bbs))，將psi設為0可去除

q/ci=(j/hi)+b*(q(+1)/ci(+1));//採用（hi=1-h）設定

h=1-hi;

mH=1-(2*(1-mF)*(1-alp)*q(+1)*h)/qshs;//alpha is endogeneous variable

///Foreign

nu*yf+bHf+bFf=
bHf(-1)*R(-1)+bFf(-1)*R(-1)+qf*(hf-hf(-1))+cf;//廠商的預算限制式

yf=Af*(hf(-1)^nu)*((1-nu)*yf/(tau*cif))^((1-nu)/(eta));//廠商生產函數
bHf=alpf*mH*qf(+1)*hf/R;//廠商和國內借貸面臨之借貸限制式

bFf=(qf(+1)*(1-alpf)*hf*(1-(1-mF)*(1-alpf)*qf(+1)*hf/qsfhsf))/R;//廠商和國外借貸面臨之借貸限制式

qf=cf*(1/(R*cf)-g/cf(+1))*mH*qf(+1)+cf*(g/cf(+1))*((nu*yf(+1)/hf)+qf(+1));//廠商最適房地產需求

cif+bHf+bF-bb-qf*(hf-hf(-1))-psi*(bb-bbs)^2/2=
bHf(-1)*R(-1)+bF(-1)*R(-1)-R(-1)*bb(-1)+(1-nu)*yf;//家戶的預算限制式，wl已被一階最適汰換

(1/cif)*(1+psi*(bb-bbs))=b*(R/cif(+1));//加入投資調整成本(1/cif)*(1+psi*(bb-bbs))，將psi設為0可去除

qf/cif=(j/hif)+b*(qf(+1)/cif(+1));//採用（hi=1-h）設定

hf=1-hif; 

mH=1-(2*(1-mF)*(1-alpf)*qf(+1)*hf)/qsfhsf;//foreign alpha is endogeneous variable



log(A)=rho_A*log(A(-1))+tau_Af*log(Af(-1));    //log(A)=rho_A*log(A(-1))+tau_Af*Af(-1);

log(Af)=rho_Af*log(Af(-1))+tau_A*log(A(-1))+e_Af;


end;

initval;
c=cs;
ci=cis;
cf=csf;
cif=cisf;
y=ys;
yf=ysf;
q=qs;
qf=qsf;
h=hs;
hi=his;
hf=hsf;
hif=hisf;
bH=bHs;
bF=bFs;
bHf=bHsf;
bFf=bFsf;
R=Rs;
bb=bbs;
alp=alps;
alpf=alpsf;
A=As;
Af=Asf;


end;

resid;
check;steady;

shocks;

var e_Af; stderr 1;


end;



stoch_simul(order=1,hp_filter=1600,irf=10,nograph);
//stoch_simul(order=1,periods=200000,irf=10);





