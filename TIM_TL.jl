#
#  ************************************************************
#  *         Code by Prof. Dr. Xue-Feng Zhang at CQU          *
#  *      To calculate the 1-d Transeverse Ising model        *
#  *                    Loop Algorithm                        *
#  ************************************************************
#---------Begining of part1 RNG--------
#Initialization of the 64-bits LXM LCG
seed::Int64=2022;  #The seed of RNG LCG64 range from [-9223372036854775808 9223372036854775807]
#The parameter of the LCG64
const lcg_a=Int64(2862933555777941757) 
const lcg_c=Int64(3037000493)
#The precision of the LCG64
const lcg_eps=-0.5/(Int64(2)^63)
#The function of the LCG64
function lcg(x0)
    x0=lcg_a*x0+lcg_c;    #The equation of the iteration
    return x0;            #X_n -> X_{n+1}
end
function rndm()           #The rng for float64 in [0. 1.]
    global seed=lcg(seed);
    return seed*lcg_eps+0.5
end
function rndm(n::Int64)   #The rng for Int64 in [1 n]
    global seed=lcg(seed);
    return min(floor(Int64,(seed*lcg_eps+0.5)*n),n-1)+1
end
#-----------End of part1 RNG-------------
#  ************************************************************
#----------Begining of part2 lattice-----
# The parameter of the lattice
ns::Int64=20; #Number of the sites or the system length
nb=ns;        #Number of bond is same as the number of sites
#The information of bond with periodical boundary condition
bond=zeros(Int,2,nb);bond[1,:]=1:ns;bond[2,:]=[2:ns;1];
#Initial the spin configuration with random states
conf=[rndm(2)-1 for i=1:ns];
#storing the value of the spin 0 for spin down and 1 for spin up
#-----------End of part2 lattice----------
#  ************************************************************
#-----------Begining of part3 weight--------
# Parameter of the physical system
Γ::Float64=0.2;      #Transverse field
B::Float64=0.5;      #Longitudinal field
V::Float64=1.;       #Ising interaction
beta::Float64=50.0;  #The inverse temperature 1/T
#--------Initialization of vertex------
# Configuration of vertexes
# Bond operators are 1: [0 0 0 0] 2: [1 0 1 0] 3: [0 1 0 1] 4:[1 1 1 1]
# Site operators I operators 5: [0 0] 8:[1 1] off-diagonal operator 6: [1 0] 7: [0 1]
# Bond operators
weight_b=zeros(4); energy_shift=V/4+B/2+.5;
# Diagonal operator
weight_b[1]=-V/4-B/2;weight_b[2]=V/4;
weight_b[3]=V/4;weight_b[4]=-V/4+B/2;
weight_b[1:4]=weight_b[1:4].+energy_shift;
#--------End of part3 vertex and weight----------------
#  ************************************************************
#--------Begining of part4 Initial operator list-------------
# Initialization of the operator list
#The length for storing the operator list
ll::Int=floor(Int,beta*nb*4); 
#The truncation number lm
lm::Int=10;
#The nubmer of non-zero operators, there is no non-zero operators at beginning
nh::Int=0;
#The operator list, index 1 store the type of vertex [1:8]
#index 2 store the position  of operator
#(0 means zero operator)
opl=zeros(Int,2,ll); 
#--------End of part4 Initial operator list----------------
#  ************************************************************
#--------Begining of part5 dupdate-------------
function dupdate()
opn=nb+ns; #number of the possible operators
for i=1:lm
    #The type of the vertex    
    vtp=opl[1,i];      
    if vtp==0          #zero operator
        r=rndm(opn);   #select one random position
        if r≤nb    # Bond operater
            tp=conf[bond[1,r]]+conf[bond[2,r]]*2+1;  #binary representation
            ap=(weight_b[tp]*beta*opn)/(lm-nh); #The probability
            if ap>rndm()   #zero -> non-zero
                #storing the vertex information into the operator list
                opl[1,i]=tp;opl[2,i]=r;    
                #non-zero operator plus one
                global nh=nh+1;                   
            end
        else             # Site operaters 
            r=r-nb;
            ap=(Γ*beta*opn)/(lm-nh); #The probability
            if ap>rndm()   #zero -> non-zero
                #storing the vertex information into the operator list
                opl[1,i]=conf[r]*3+5;opl[2,i]=r;    
                #non-zero operator plus one
                global nh=nh+1;                   
            end
        end
    elseif vtp≠6 && vtp≠7       #diagonal operator
        if vtp<5 # Bond operator
            ap=(lm-nh+1)/(weight_b[vtp]*beta*opn);
        else
            ap=(lm-nh+1)/(Γ*beta*opn);
        end
        if ap>rndm()   #nonzero -> zero
            #ereasing the vertex information
            opl[1,i]=0;opl[2,i]=0;
            #non-zero operator minus one
            global nh=nh-1;                   
        end
    else           #off-diagonal operator  
        r=opl[2,i];
        conf[r]=1-conf[r];
    end
end
end
#--------End of part5 dupdate----------------
#  ************************************************************
#--------Begining of part6 lupdate-------------
function lupdate()
    #--------Firstly construct the link table and merging the operator
    opl3=zeros(Int,2,nh);  #The operator list without zero operators
    link=zeros(Int,4*nh);  #Link table
    ft=-ones(Int,ns);      #storing the free legs at the beginning
    lt=-ones(Int,ns);      #storing the free legs at the end (-1 mean no leg connected)
    is=0;  #The index of non-zero operator excluding the unitary operator
    ln=0;  #The number of legs
    nm=0; #The number of merging operator
    pos_nm=zeros(Int,nh);  #storing the τ-position of merging operator
    opcfg=zeros(Int,nh*4);  #storing the configuration of the link
    for i=1:lm
        tp=opl[1,i]; #The type of the vertex
    	#!!!!Note that the type in opl3 1:4 bond 0: constant -1 off-diagonal operator
        if tp≠0   #Non-zero operator
            r=opl[2,i]; #position of the operator
            if tp <5 # Bond operator
                is=is+1;
                opl3[:,is]=opl[:,i];
                opcfg[ln+1]=conf[bond[1,r]];opcfg[ln+2]=conf[bond[2,r]]
                opcfg[ln+3]=conf[bond[1,r]];opcfg[ln+4]=conf[bond[2,r]]
            else #    Site operator
                is=is+1;
                if rndm()>0.5   #The appending site is in the left
                    if r==1
                        r=ns;
                    else
                        r=r-1;
                    end
                end
                opcfg[ln+1]=conf[bond[1,r]];opcfg[ln+2]=conf[bond[2,r]]
                if tp==6 || tp==7   #off-diagonal
                    opl3[1,is]=-1;
                    conf[opl[2,i]]=1-conf[opl[2,i]];  #update the configuration
                else
                    opl3[1,is]=0;
                end
                opcfg[ln+3]=conf[bond[1,r]];opcfg[ln+4]=conf[bond[2,r]]
                opl3[2,is]=r;
                nm=nm+1;
                pos_nm[nm]=is;
            end
            for leg=1:2
                bl=bond[leg,r];
                if  ft[bl]==-1; ft[bl]=ln+leg;end
                if  lt[bl]==-1
                    lt[bl]=ln+2+leg; #storing the last slide
                else
                     #Link the operator before and after
                    link[lt[bl]]=ln+leg; 
                    link[ln+leg]=lt[bl];  
                    #move to next slide
                    lt[bl]=ln+2+leg;      
                end
           end
           ln=ln+4; #Increase the number of legs
        end
    end
    #After that, link the legs at the boundary
    for i=1:ns
        if ft[i]≠-1
            link[ft[i]]=lt[i];
            link[lt[i]]=ft[i];
        end
    end
    #starting process
    exitP=0.5;
    pass_leg=[3 4 1 2];
    pass_off=[2 3 4;1 3 4; 1 2 4; 1 2 3];
    for n_loop=1:nm
        st=pos_nm[rndm(nm)];# randomly choose a single site operator to start the update
        vtx0=opl3[1,st];
        if vtx0==0 #i-operator
            # random choose a leg to start
            # because constant operator has same probability to appeart on the two side, so it has same probability to start at 4 leg
            outleg=rndm(4);
        elseif opcfg[(st-1)*4+1]≠opcfg[(st-1)*4+3] #left-off-operator
            outleg=rndm(2)*2-1
        else    #right-operator
            outleg=rndm(2)*2
        end
        opl3[1,st]=-1-vtx0;j0=(st-1)*4+outleg;
        opcfg[j0]=1-opcfg[j0];
        j2::Int64=j0;
        fflag=true;
        while fflag       
            j1=link[j2];
            inleg=mod((j1-1),4)+1;
            st=floor(Int,(j1-1)/4)+1;  #position of the vertex
            vtx0=opl3[1,st];
            if vtx0>0 #bond-operator
                if isodd(inleg)
                    vtx2=(1-opcfg[(st-1)*4+1])+opcfg[(st-1)*4+2]*2+1;
                else
                    vtx2=opcfg[(st-1)*4+1]+(1-opcfg[(st-1)*4+2])*2+1;
                end    
                ap=rndm();
                if ap<weight_b[vtx2]/weight_b[vtx0]
                    opcfg[j1]=1-opcfg[j1];
                    outleg=pass_leg[inleg];j2=j1-inleg+outleg;
                    opcfg[j2]=1-opcfg[j2];
                    opl3[1,st]=vtx2
                else
                    j2=j1;
                end
            elseif vtx0==0 #i-operator
                ap=rndm();opcfg[j1]=1-opcfg[j1];
                if ap≤0.5*exitP # it has a 1/2 probability to meet the current position of the const operator
                    opl3[1,st]=-1-vtx0;
                    fflag=false;
                else   #passing through
                    outleg=pass_leg[inleg];j2=j1-inleg+outleg;
                    opcfg[j2]=1-opcfg[j2];                
                end
            elseif opcfg[(st-1)*4+1]≠opcfg[(st-1)*4+3] #left-off-operator
                opcfg[j1]=1-opcfg[j1];
                if isodd(inleg)&(rndm()≤exitP)  #possible stop
                    opl3[1,st]=-1-vtx0;
                    fflag=false;
                else     #passing through
                    an=rndm(3);outleg=pass_off[inleg,an];j2=j1-inleg+outleg;
                    opcfg[j2]=1-opcfg[j2];
                end
            else     #right-operator
                opcfg[j1]=1-opcfg[j1];
                if iseven(inleg)&(rndm()≤exitP)  #possible stop
                    opl3[1,st]=-1-vtx0;
                    fflag=false;
                else    #passing through
                    an=rndm(3);outleg=pass_off[inleg,an];j2=j1-inleg+outleg;
                    opcfg[j2]=1-opcfg[j2];
                end
            end 
        end
    end
    #---------------recovering the configuration---------
    for i=1:ns
        if ft[i]≠-1 
            conf[i]=opcfg[ft[i]];
        end
    end
    #--------------updating the original operator list-------
    #               disjoin the operator
    is=0;
    for i=1:lm
        #The type of the vertex    
        vtp=opl[1,i];   
        if vtp≠0
            if vtp<5 #Bond operator
                is=is+1;
                opl[:,i]=opl3[:,is];
            else               #site operator
                is=is+1;
                r=opl3[2,is];
                vtx0=opl3[1,is];
                if vtx0==0 #i-operator
                    if rndm()<0.5 #left side left
                        opl[1,i]=opcfg[(is-1)*4+1]*3+5;
                        opl[2,i]=bond[1,r];
                    else #right side left
                        opl[1,i]=opcfg[(is-1)*4+2]*3+5;
                        opl[2,i]=bond[2,r];
                    end
                elseif opcfg[(is-1)*4+1]≠opcfg[(is-1)*4+3] #left-off-operator
                    opl[1,i]=7-opcfg[(is-1)*4+1]
                    opl[2,i]=bond[1,r];
                else #right operator
                    opl[1,i]=7-opcfg[(is-1)*4+2]
                    opl[2,i]=bond[2,r];
                end
            end
        end
        if is==nh; break;end
    end
end
#--------End of part6 lupdate----------------
#  ************************************************************
#--------Begining of part7 Main-------------
istp=500;         #The number of thermalization step
mstp=1000;        #The number of Measuring steps
    en=zeros(mstp);
    mag=zeros(mstp);
    mag_s=zeros(mstp);
for i=1:istp
    dupdate();
    lupdate();
    lt=floor(Int,1.25*nh);
    if lt>lm 
        global lm=lt;
    end
end
for i=1:mstp
    dupdate();
    lupdate();
    en[i]=nh;
    rho1=sum(conf[1:2:ns]);
    rho2=sum(conf[2:2:ns]);
    mag[i]=rho1+rho2;
    mag_s[i]=(rho1-rho2)^2;
end
display("The Energy is:"*string(sum(en.*1.)/mstp/beta/nb-energy_shift-Γ))
display("The magnetization is:"*string(sum(mag.*1.)/mstp/nb))
display("The stagger magnetization is:"*string(sum(mag_s.*1.)/mstp/nb))
