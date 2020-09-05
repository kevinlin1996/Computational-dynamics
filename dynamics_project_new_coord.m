% Kevin Lin
% Determine dynamics for 3 jointed RRP manipulator for project
% Inputs: degree of freedom, DH parameters,
% gravity,type of joint, local centre of mass vectors, and inertial
% matricies
% Outputs: Dynamic model of 3 DoF RRP manipulator 

clear
clc

%% Input DH param, create DH tables, and transformation matricies
joints = input('Enter the type of serial robot to determine Jacobian(RRR, PRP, etc.): ','s');
n = strlength(joints);

          
DHt = cell(n,4); % initialize DH table
for i=1:n % create DH table row by row
            
    fprintf('\nFor row %d of the DH table (can be numeric or symbolic):\n\n',i)
    d_t = input('Enter the value for d: ','s');
    th_t = input('Enter the value for theta in degrees: ','s');
    a_t = input('Enter the value for a: ','s');
    alf_t = input('Enter the value for alpha in degrees: ','s');
    DHt(i,:) = {d_t th_t a_t alf_t}; 
end

DHt = string(DHt);

%Initialization
sz = size(DHt);
d = sym('d%d',[1 n]); % lists syms d1 to dn
th = sym('th%d',[1 n]); % lists syms th1 to thn
a = sym('a%d',[1 n]); % lists syms a1 to an
alf = sym('alf%d',[1 n]); % lists syms alf1 to alfn
q=sym('q',[n,1],'real'); % joint variables

DHtd = str2double(DHt); % converts string to double, strings become NaN
d_d = DHtd(1:n,1);
th_d = deg2rad(DHtd(1:n,2)); % convert deg to rad for the DH matrix function
a_d = DHtd(1:n,3);
alf_d = deg2rad(DHtd(1:n,4));
T=sym('T',[4,4],'real');
Ti = sym(zeros(4,4)); % initialize T (contains T matrices between frames)
T0n = sym(eye(4,4)); % initialize T0n (contains T matrices of frames WRT 0 frame)


for i=1:n % runs loop for DH matrix and jacobian components
    % checks if values in DH are numbers to know what values to sub    
    if (joints(i)=='P' || joints(i)=='p')
        d(i) = q(i);
        fprintf('hey');
    elseif (isnan(d_d(i)) == 0) % if input is a number then store in corresponding variable
        d(i) = d_d(i);
    end
    
    if (joints(i)=='R'|| joints(i)=='r')
        th(i)=q(i);
    elseif (isnan(th_d(i)) == 0)
        th(i) = th_d(i);
    end
    
    if isnan(a_d(i)) == 0
        a(i) = a_d(i);
    end
    
    if isnan(alf_d(i)) == 0
        alf(i) = alf_d(i);
    end
    
    % creates transformation matrices
    Ti(:,:,i) = DHmatrix(d(i),th(i),a(i),alf(i)); % T matrix
    T0n(:,:,i) = simplify(T0n(:,:,i)*Ti(:,:,i)); % multiplies T0 to Tn saving between matrices
    T0n(:,:,i+1) = T0n(:,:,i);  %To be used in next iteration
    
end

T(:,:,1) = eye(4,4); % First element is T00
for i=1:n
    T(:,:,i+1)=T0n(:,:,i);
end

%% Variables and constants 

I=sym('I',[4,4],'real');
U=sym('U',[4,4],'real');
UU = sym('UU',[4,4],'real');
M=sym('M',[1,1],'real'); 
Gam=sym('Gam',[1,1],'real');
C = sym('C',[1,1],'real'); 
G = sym('G',[1,1],'real'); 
F = sym('F',[n,1],'real');

Qr=[0 -1 0 0; 1 0 0 0; 0 0 0 0; 0 0 0 0];
Qp=[0 0 0 0; 0 0 0 0; 0 0 0 1; 0 0 0 0];


mass=sym('m',[n,1],'real');
mass(1) = 0.15;
mass(2) = 0.122;
mass(3) = 0.0295;

dq=sym('dq',[n,1],'real');
ddq=sym('ddq',[n,1],'real');

% gravity stuff
gx = 0;
gy = 0; 
gz = -9.81; % m/s^2 wrt global frame F0
g = [gx,gy,gz,0];

Q = zeros (4,4,n);

% Creates Q arrays depending on the joint
for i=1:n
    % joints char array
    if joints(i)=='R' || joints(i)=='r'
      Q(:,:,i)=Qr;
    else 
      Q(:,:,i)=Qp;
    end
end


% local centre of mass vectors (*** change/add depending on robot)
rbar(:,:,1) = [0;-0.06238;0.00062;1];  %vector from local sys to mass of the 1st link
rbar(:,:,2) = [0;-0.000767;0.059;1]; 
rbar(:,:,3) = [0;0;-0.0297;1];

% Inertial matricies (*** Change/add depending on robot)
I(:,:,1)=[0.0000059,	0,	0,	0;
0,	0.00017733,	0.00000533,	-0.009357;
0,	0.00000533,	0.0000058,	0.000093;
0,	-0.009357,	0.000093,	0.15];
  
I(:,:,2)=[0.00000559,	0,	0,	0;
0,	0.0000049,	0.00000502,	-0.000093574;
0,	0.00000502,	0.00015913,	0.007198;
0,	-0.000093574,	0.007198,	0.122];
  
I(:,:,3) = [0.000000805,	0,	0,	0;
0,	0.000000615,	0,	0;
0,	0,	0.000023065,	-0.00087615;
0,	0,	-0.00087615,	0.0295];
 


%    --- Begin calculation of required matricies for dynamics ---

% U matrix, i = rows, j = columns 4 x 4 matricies placed within an i x j
% matrix
for i=1:n
    for j=1:n
          
        % when n = 1, only executes j = 1, not j = 2, see condition for U matrix in sheet 
        % j <= i
        if j<=i
            % starts at (j-1), aka 0. j here = j-1. for i, it is i+1
            % (same with j) for transformation matrix only! 
            % start at T01 (T 0 to j-1)
            % Q1 = Q(:,:,1), T01 = T(:,:,2)
            U(:,:,i,j)=simplify(T(:,:,j)*Q(:,:,j)*inv(T(:,:,j))*T(:,:,i+1));
        % j > i 
        else 
            U(:,:,i,j)=zeros(4,4);
        end
    end
end

% M matrix i = row, k = column note that this is just a ixk matrix with
% scalar terms M is n x n 
for i=1:n
    for k=1:n
        M(:,:,i,k)=0;       
        
        for j=max(i,k):n
            % prev plus new trace for summation 
            M(:,:,i,k)=M(:,:,i,k)+trace(U(:,:,j,k)*I(:,:,j)*U(:,:,j,i)');
        end
    end
end

 M=simplify(M);
 
% UU matrix
% 3 dimensions: 4 x 4 matricies within each cell of ixjxk array 
  for i=1:n
      for j=1:n
          for k=1:n
               
              % i>=k>=j
              if (i>=k && k>=j)
                  % starts at (j-1), aka 0. j here = j-1. for i, it is i+1
                  % (same with j) for transformation matrix only! 
                  % Q1 = Q(:,:,1), T01 = T(:,:,2)
                  UU(:,:,i,j,k)=simplify(T(:,:,j)*Q(:,:,j)*inv(T(:,:,j))*T(:,:,k)*Q(:,:,k)*inv(T(:,:,k))*T(:,:,i+1));
                  
              % i>=j>=k
              elseif (i>=j && j>=k)
                  UU(:,:,i,j,k)=simplify(T(:,:,k)*Q(:,:,k)*inv(T(:,:,j))*T(:,:,k)*Q(:,:,j)*inv(T(:,:,j))*T(:,:,i+1));
                  
              % i<j or i<k
              else
                  UU(:,:,i,j,k)=zeros(4,4);
              end
          end
      end
  end
  
  % Gamma matrix
  
  for i=1:n
      for k=1:n
           for m=1:n
               
               Gam(:,:,i,k,m)=0;       
        
               for j=max(i,max(k,m)):n
                   Gam(:,:,i,k,m)=Gam(:,:,i,k,m)+trace(UU(:,:,j,k,m)*I(:,:,j)*U(:,:,j,i)');
                
               end
           end
      end
  end
  

 % C matrix, C is n x n 
 for i=1:n
      for j=1:n
          C(:,:,i,j)=0;
           for k=1:n
               C(:,:,i,j)=C(:,:,i,j)+Gam(:,:,i,j,k)*dq(k);
           end
      end
 end
 
 C=simplify(C);
 
 % G matrix, should be n x 1 
 for i=1:n
     G(:,:,i)=0;
     for j=1:n
         G(:,:,i)=G(:,:,i)+mass(j)*g*U(:,:,j,i)*rbar(:,:,j);
     end
     G(:,:,i)=-1*G(:,:,i);
 end
 
 G=simplify(G);
 
 
 % Final dynamics model: Forces and torques (this vec is n x 1)
 
 for i=1:n
     
    F(i) = 0;
     
    for k=1:n
        F(i)=F(i)+ddq(k)*M(:,:,i,k)+dq(k)*C(:,:,i,k);
    end
    
    F(i) = F(i)+G(:,:,i);
 end
 
 
 F= vpa(simplify(F,4))
 
 %% See the M matrix and test 
 Mmatrix=sym('Mmatrix',[3,3],'real'); 
 
 for i=1:3
     for j=1:3
         Mmatrix(i,j)= M(:,:,i,j);
     end
 end
 
 vpa(simplify(Mmatrix),4)
 
 %%
 m_dot=0;
 for i=1:n
     
    m_dot = m_dot+diff(M(:,:,1,2),q(i))*dq(i)
    
 end

 x = [1 2];
 check = simplify(x*(m_dot-2*C(:,:,1,2))*x')
 
 %% Symbolic DH matrix function

function Mdh = DHmatrix(d,theta,a,alpha)
% Note: uses trig(sym(angle)) to produce 0 when required as cases involving
% cos(pi/2) produce a non-zero answer when it should be 0
Mdh = [cos(sym(theta)) -sin(sym(theta))*cos(sym(alpha)) sin(sym(theta))*sin(sym(alpha)) a*cos(sym(theta));
    sin(sym(theta)) cos(sym(theta))*cos(sym(alpha)) -cos(sym(theta))*sin(sym(alpha)) a*sin(sym(theta));
    0 sin(sym(alpha)) cos(sym(alpha)) d;
    0 0 0 1];
end

 