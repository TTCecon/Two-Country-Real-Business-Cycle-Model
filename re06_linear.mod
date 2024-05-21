// Institute of Economics, NSYSU 
// Kuan Lun Chen (2024.4)
// This is a two country RBC model (linear endogeneous alpha) 
// Replicate Iacoviello & Minetti (2006)"National Business Cycles With Domestic and Foreign Lenders"(JME)
// This replicate of the model is two country RBC model with endogeneous alpha but without investment adjustment cost 
// i代表Household面對的變數，如消費ci..
// f代表國外，如國外家戶房產數量hf
// bH,bF為國內企業家面臨之借貸限制，其s.s.為bHs,bFs
// bHf,bFf為國外企業家面臨之借貸限制，其s.s.為bHsf,bFsf
// endogeneous alpha case and no adjustment costs

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

ys*nu*y+bHs*bH+bFs*bF=
Rs*bHs*R(-1)+Rs*bFs*R(-1)+Rs*bHs*bH(-1)+Rs*bFs*bF(-1)+cs*c+qshs*(h-h(-1));//廠商的預算限制式

(1-(1-nu)/(eta))*y=A+nu*h(-1)-((1-nu)/(eta))*ci;//廠商生產函數

R+bH=alp+q(+1)+h;//廠商和國內借貸面臨之借貸限制式1/alps*(q(+1)+h)

Rs*bFs*(R+bF)=((1-alps)*(q(+1)+h)-alps*alp)*(qshs-2*qshs*(1-mF)*(1-alps));//廠商和國外借貸面臨之借貸限制式Rs*bFsf*(R+bFf)=(((-alps)*qshs))*(1-alp+q(+1)*h)*(1-2*(1-mF)*(-alps))

q-(g+(b-g)*mH)*q(+1)=
(1-g-(b-g)*mH)*(y(+1)-h)+(1-b*mH)*(c-c(+1))-(b*mH)*R;//廠商最適房地產需求，原[g-(b-g)*mH](double check q-c=y(+1)-c(+1)-h+q(+1)-c(+1)+q(+1)+h)

cis*ci+bHs*bH+bFsf*bFf+bb=
Rs*bHs*R(-1)+Rs*bFsf*R(-1)+Rs*bHs*bH(-1)+Rs*bFsf*bFf(-1)+Rs*(bb(-1))+(1-nu)*ys*y+qshs*(h-h(-1));//家戶的預算限制式，wl已被一階最適汰換

ci(+1)=ci+R;//去除投資調整成本(1/cif)*(1+psi*(bb-bbs)),將psi設為0可去除

q=b*q(+1)+(hs/(1-hs))*(1-b)*h+ci-b*ci(+1);//採用（hi=1-h）設定

hs*h=-his*hi;

alp=((1-alps)/alps)*(q(+1)+h);

///Foreign

ysf*nu*yf+bHsf*bHf+bFsf*bFf=
Rs*bHsf*R(-1)+Rs*bFsf*R(-1)+Rs*bHsf*bHf(-1)+Rs*bFsf*bFf(-1)+csf*cf+qsfhsf*(hf-hf(-1));//廠商的預算限制式

(1-((1-nu)/(eta)))*yf=Af+nu*hf(-1)-((1-nu)/(eta))*cif;//廠商生產函數

R+bHf=alpf+qf(+1)+hf;//廠商和國內借貸面臨之借貸限制式1/alpsf*(qf(+1)+hf)

Rs*bFsf*(R+bFf)=((1-alpsf)*(qf(+1)+hf)-alpsf*alpf)*(qsfhsf-2*qsfhsf*(1-mF)*(1-alpsf));//廠商和國外借貸面臨之借貸限制式(((-alpsf)*qsfhsf))*(1-alpf+qf(+1)*hf)*(1-2*(1-mF)*(-alpsf))

qf-(g+(b-g)*mH)*qf(+1)=
(1-g-(b-g)*mH)*(yf(+1)-hf)+(1-b*mH)*(cf-cf(+1))-(b*mH)*R;//廠商最適房地產需求(qf-cf=yf(+1)-cf(+1)-hf+qf(+1)-cf(+1)+qf(+1)+hf)

cisf*cif+bHsf*bHf+bFs*bF-bb=
Rs*bHsf*R(-1)+Rs*bFs*R(-1)-Rs*(bb(-1))+Rs*bHsf*bHf(-1)+Rs*bFs*bF(-1)+(1-nu)*ysf*yf+qsfhsf*(hf-hf(-1));//家戶的預算限制式，wl已被一階最適汰換

cif(+1)=cif+R;//去除投資調整成本cif(+1)=cif+R+psi*bbs*b,將psi設為0可去除

qf=b*qf(+1)+(hsf/(1-hsf))*(1-b)*hf+cif-b*cif(+1);//採用（hi=1-h）設定

hsf*hf=-hisf*hif;

alpf=((1-alpsf)/alpsf)*(qf(+1)+hf);

A=rho_A*A(-1)+tau_Af*Af(-1);

Af=rho_Af*Af(-1)+tau_A*A(-1)+e_Af;

end;

initval;
c=0;
ci=0;
cf=0;
cif=0;
y=0;
yf=0;
q=0;
qf=0;
h=0;
hf=0;
bH=0;
bF=0;
bHf=0;
bFf=0;
R=0;
bb=0;
alp=0;
alpf=0;
A=0;
Af=0;


end;

resid;
check;steady;

shocks;

var e_Af; stderr 1;


end;


stoch_simul(order=1,hp_filter=1600,irf=10,nograph);
//stoch_simul(order=1,periods=200000,irf=10);