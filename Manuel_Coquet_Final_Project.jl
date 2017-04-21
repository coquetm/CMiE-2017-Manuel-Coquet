
using JuMP, Ipopt, Plots, DataFrames
gr()

# Since we are not changing resistance or reactance, the admittance matrix will always be constant

Zij = zeros(5,5)
Zij = complex(Zij)
Zij[1,2] = 1+1im; Zij[2,1] = 1+1im;
Zij[2,3] = 1+1im; Zij[3,2] = 1+1im;
Zij[3,4] = 1+1im; Zij[4,3] = 1+1im;
Zij[4,5] = 1+1im; Zij[5,4] = 1+1im;
Zij[5,1] = 1+1im; Zij[1,5] = 1+1im;
Y = zeros(5,5);
Y = complex(Y)
for i=1:5, j=1:5
    if i != j
        if Zij[i,j] != 0
            Y[i,j] = -Zij[i,j]^-1;
        else
            Y[i,j] = 0;
        end
    end
end
for i = 1:5
    Y[i,i] = -sum(Y[i,:]);
end
Y

function feval(Y, V; f0 = zeros(7))

function f(x::Vector)
fx = zeros(7)
θ₁ = 0; #slack reference bus
fx[1] = V[2]*V[1]*(real(Y[2,1])*cos(x[1]-θ₁)+imag(Y[2,1])*sin(x[1]-θ₁))+ V[2]*V[2]*real(Y[2,2])+
        V[2]*x[5]*(real(Y[2,3])*cos(x[1]-x[2])+imag(Y[2,3])*sin(x[1]-x[2])) - f0[1]
fx[2] = x[5]*V[2]*(real(Y[3,2])*cos(x[2]-x[1])+imag(Y[3,2])*sin(x[2]-x[1]))+ x[5]*x[5]*real(Y[3,3])+
        x[5]*x[6]*(real(Y[3,4])*cos(x[2]-x[3])+imag(Y[3,4])*sin(x[2]-x[3])) - f0[2]
fx[3] = x[6]*x[5]*(real(Y[4,3])*cos(x[3]-x[2])+imag(Y[4,3])*sin(x[3]-x[2]))+ x[6]*x[6]*real(Y[4,4])+
        x[6]*x[7]*(real(Y[4,5])*cos(x[3]-x[4])+imag(Y[4,5])*sin(x[3]-x[4])) - f0[3]
fx[4] = x[7]*x[6]*(real(Y[5,4])*cos(x[4]-x[3])+imag(Y[5,4])*sin(x[4]-x[3]))+ x[7]*x[7]*real(Y[5,5])+
        x[7]*V[1]*(real(Y[5,1])*cos(x[4]-θ₁)+imag(Y[5,4])*sin(x[4]-θ₁)) - f0[4]
fx[5] = x[5]*V[2]*(real(Y[3,2])*sin(x[2]-x[1])-imag(Y[3,2])*cos(x[2]-x[1])) -x[5]*x[5]*imag(Y[3,3])+
        x[5]*x[6]*(real(Y[3,4])*sin(x[2]-x[3])-imag(Y[3,4])*cos(x[2]-x[3])) - f0[5]
fx[6] = x[6]*x[5]*(real(Y[4,3])*sin(x[3]-x[2])-imag(Y[4,3])*cos(x[3]-x[2]))-x[6]*x[6]*imag(Y[4,4])+
        x[6]*x[7]*(real(Y[4,5])*sin(x[3]-x[4])-imag(Y[4,5])*cos(x[3]-x[4])) - f0[6]
fx[7] = x[7]*x[6]*(real(Y[5,4])*sin(x[4]-x[3])-imag(Y[5,4])*cos(x[4]-x[3]))-x[7]*x[7]*imag(Y[5,5])+
        x[7]*V[1]*(real(Y[5,1])*sin(x[4]-θ₁)-imag(Y[5,1])*cos(x[4]-θ₁)) - f0[7]
return fx
end
    
end

function Jeval(Y, V)

function J(x::Vector)
θ₁ = 0 #slack reference bus
df1x1 = V[2]*V[1]*(-real(Y[2,1])*sin(x[1]-θ₁)+imag(Y[2,1])*cos(x[1]-θ₁))+
            V[2]*x[5]*(-real(Y[2,3])*sin(x[1]-x[2])+imag(Y[2,3])*cos(x[1]-x[2]))
df1x2 = V[2]*x[5]*(real(Y[2,3])*sin(x[1]-x[2])-imag(Y[2,3])*cos(x[1]-x[2]))
df1x3 = 0
df1x4 = 0
df1x5 = V[2]*(real(Y[2,3])*cos(x[1]-x[2])+imag(Y[2,3])*sin(x[1]-x[2]))
df1x6 = 0
df1x7 = 0
df2x1 = x[5]*V[2]*(real(Y[3,2])*sin(x[2]-x[1])-imag(Y[3,2])*cos(x[2]-x[1]))
df2x2 = x[5]*V[2]*(-real(Y[3,2])*sin(x[2]-x[1])+imag(Y[3,2])*cos(x[2]-x[1]))+
            x[5]*x[6]*(-real(Y[3,4])*sin(x[2]-x[3])+imag(Y[3,4])*cos(x[2]-x[3]))
df2x3 = x[5]*x[6]*(real(Y[3,4])*sin(x[2]-x[3])-imag(Y[3,4])*cos(x[2]-x[3]))
df2x4 = 0;
df2x5 = V[2]*(real(Y[3,2])*cos(x[2]-x[1])+imag(Y[3,2])*sin(x[2]-x[1]))+
            2*x[5]*real(Y[3,3])+x[6]*(real(Y[3,4])*cos(x[2]-x[3])+imag(Y[3,4])*sin(x[2]-x[3]))
df2x6 = x[5]*(real(Y[3,4])*cos(x[2]-x[3])+imag(Y[3,4])*sin(x[2]-x[3]))
df2x7 = 0
df3x1 = 0
df3x2 = x[6]*x[5]*(real(Y[4,3])*sin(x[3]-x[2])-imag(Y[4,3])*cos(x[3]-x[2]))
df3x3 = x[6]*x[5]*(-real(Y[4,3])*sin(x[3]-x[2])+imag(Y[4,3])*cos(x[3]-x[2]))+
            x[6]*x[7]*(-real(Y[4,5])*sin(x[3]-x[4])+imag(Y[4,5])*cos(x[3]-x[4]))
df3x4 = x[6]*x[7]*(real(Y[4,5])*sin(x[3]-x[4])-imag(Y[4,5])*cos(x[3]-x[4]))
df3x5 = x[6]*(real(Y[4,3])*cos(x[3]-x[2])+imag(Y[4,3])*sin(x[3]-x[2]))
df3x6 = x[5]*(real(Y[4,3])*cos(x[3]-x[2])+imag(Y[4,3])*sin(x[3]-x[2]))+
            2*x[6]*real(Y[4,4])+x[7]*(real(Y[4,5])*cos(x[3]-x[4])+imag(Y[4,5])*sin(x[3]-x[4]))
df3x7 = x[6]*(real(Y[4,5])*cos(x[3]-x[4])+imag(Y[4,5])*sin(x[3]-x[4]))
df4x1 = 0
df4x2 = 0
df4x3 = x[7]*x[6]*(real(Y[5,4])*sin(x[4]-x[3])-imag(Y[5,4])*cos(x[4]-x[3]))
df4x4 = x[7]*x[6]*(-real(Y[5,4])*sin(x[4]-x[3])+imag(Y[5,4])*cos(x[4]-x[3]))+
            x[7]*V[1]*(-real(Y[5,1])*sin(x[4]-θ₁)+imag(Y[5,1])*cos(x[4]-θ₁))
df4x5 = 0
df4x6 = x[7]*(real(Y[5,4])*cos(x[4]-x[3])+imag(Y[5,4])*sin(x[4]-x[3]))
df4x7 = x[6]*(real(Y[5,4])*cos(x[4]-x[3])+imag(Y[5,4])*sin(x[4]-x[3]))+
            2*x[7]*real(Y[5,5])+V[1]*(real(Y[5,1])*cos(x[4]-θ₁)+imag(Y[5,1])*sin(x[4]-θ₁))
df5x1 = x[5]*V[2]*(-real(Y[3,2])*cos(x[2]-x[1])-imag(Y[3,2])*sin(x[2]-x[1]))
df5x2 = x[5]*V[2]*(real(Y[3,2])*cos(x[2]-x[1])+imag(Y[3,2])*sin(x[2]-x[1]))+
                        x[5]*x[6]*(real(Y[3,4])*cos(x[2]-x[3])+imag(Y[3,4])*sin(x[2]-x[3]))
df5x3 = x[5]*x[6]*(-real(Y[3,4])*cos(x[2]-x[3])-imag(Y[3,4])*sin(x[2]-x[3]))
df5x4 = 0;
df5x5 = V[2]*(real(Y[3,2])*sin(x[2]-x[1])-imag(Y[3,2])*cos(x[2]-x[1]))-
            2*x[5]*imag(Y[3,3])+x[6]*(real(Y[3,4])*sin(x[2]-x[3])-imag(Y[3,4])*cos(x[2]-x[3]))
df5x6 = x[5]*(real(Y[3,4])*sin(x[2]-x[3])-imag(Y[3,4])*cos(x[2]-x[3]))
df5x7 = 0
df6x1 = 0
df6x2 = x[6]*x[5]*(-real(Y[4,3])*cos(x[3]-x[2])-imag(Y[4,3])*sin(x[3]-x[2]))
df6x3 = x[6]*x[5]*(real(Y[4,3])*cos(x[3]-x[2])+imag(Y[4,3])*sin(x[3]-x[2]))+
            x[6]*x[7]*(real(Y[4,5])*cos(x[3]-x[4])+imag(Y[4,5])*sin(x[3]-x[4]))
df6x4 = x[6]*x[7]*(-real(Y[4,5])*cos(x[3]-x[4])-imag(Y[4,5])*sin(x[3]-x[4]))
df6x5 = x[6]*(real(Y[4,3])*sin(x[3]-x[2])-imag(Y[4,3])*cos(x[3]-x[2]))
df6x6 = x[5]*(real(Y[4,3])*sin(x[3]-x[2])-imag(Y[4,3])*cos(x[3]-x[2]))-
            2*x[6]*imag(Y[4,4])+x[7]*(real(Y[4,5])*sin(x[3]-x[4])-imag(Y[4,5])*cos(x[3]-x[4]))
df6x7 = x[6]*(real(Y[4,5])*sin(x[3]-x[4])-imag(Y[4,5])*cos(x[3]-x[4]))
df7x1 = 0
df7x2 = 0
df7x3 = x[7]*x[6]*(-real(Y[5,4])*cos(x[4]-x[3])-imag(Y[5,4])*sin(x[4]-x[3]))
df7x4 = x[7]*x[6]*(real(Y[5,4])*cos(x[4]-x[3])+imag(Y[5,4])*sin(x[4]-x[3]))+
            x[7]*V[1]*(real(Y[5,1])*cos(x[4]-θ₁)+imag(Y[5,1])*sin(x[4]-θ₁))
df7x5 = 0
df7x6 = x[7]*(real(Y[5,4])*sin(x[4]-x[3])-imag(Y[5,4])*cos(x[4]-x[3]))
df7x7 = x[6]*(real(Y[5,4])*sin(x[4]-x[3])-imag(Y[5,4])*cos(x[4]-x[3]))-
            2*x[7]*imag(Y[5,5])+V[1]*(real(Y[5,1])*sin(x[4]-θ₁)-imag(Y[5,1])*cos(x[4]-θ₁))
Jx = [
    df1x1 df1x2 df1x3 df1x4 df1x5 df1x6 df1x7;
    df2x1 df2x2 df2x3 df2x4 df2x5 df2x6 df2x7;
    df3x1 df3x2 df3x3 df3x4 df3x5 df3x6 df3x7;
    df4x1 df4x2 df4x3 df4x4 df4x5 df4x6 df4x7;
    df5x1 df5x2 df5x3 df5x4 df5x5 df5x6 df5x7;
    df6x1 df6x2 df6x3 df6x4 df6x5 df6x6 df6x7;
    df7x1 df7x2 df7x3 df7x4 df7x5 df7x6 df7x7;
]
return Jx
                                
end

end

V = 0.9+0.2*rand(5,1) # random number between 0.9-1.1
θ = (-2+2*2*randn(5,1))*(π/180) # random angle between -2º to 2º
V[1] = 1 # slack bus voltage magnitude (reference)
θ[1] = 0 # slack bus voltage angle (reference)

f0 = feval(Y,V) # generates a function to evaluate power flows in terms of x
x_true = [θ[2]; θ[3]; θ[4]; θ[5]; V[3]; V[4]; V[5]]
constraints = f0(x_true) # computes constraints based on initial vector x0

function PF_newton(feval, Jeval, V, θ, Y; x0 = [0,0,0,0,1,1,1],tol = 1e-9, xtol = 1e-4, max_iter = 15)

x_true = [θ[2]; θ[3]; θ[4]; θ[5]; V[3]; V[4]; V[5]] # original x vector (angles in radians)
    
f0 = feval(Y,V) # generates a function to evaluate power flows in terms of x
constraints = f0(x_true) # computes constraints based on initial vector x0

# compute new fx function incorporating constraints and Jacobian
fx = feval(Y,V, f0 = constraints)
Jx = Jeval(Y,V)

#initialize algorithm
x = x0

for i = 1:max_iter
    xn = x - Jx(x)\fx(x)
    x = xn

# Was a solution found?
if sum(fx(x).^2) < tol
    break
end
    
end

# Is it the correct soultion? It is possible to have more than 1 local solution
if sum((x-x_true).^2) < xtol
    success_newt = true
else
    success_newt = false
end

    return success_newt, x
end

(success_newt, x_newt) = PF_newton(feval, Jeval, V, θ, Y)

function PF_JuMP(feval, V, θ, Y; xtol = 1e-4)

x_true = [θ[2]; θ[3]; θ[4]; θ[5]; V[3]; V[4]; V[5]] # original x vector (angles in radians)    

# JuMP does not recognize real and imaginary functions so I have to define the Real and Imaginary Parts of Y
YR = real(Y); YI = imag(Y); θ₁ = 0;
    
f0 = feval(Y,V) # generates a function to evaluate power flows in terms of x
constraints = f0(x_true) # computes constraints based on initial vector x0

# Initialize Ipopt Solver
m = Model(solver=IpoptSolver(print_level=0))
    
# Define our x variable
@variable(m, x[1:7])

# Set initial guess
for i in 1:4
    setvalue(x[i], 0.0)
end

for i in 5:7
    setvalue(x[i], 1.0)
end


@NLobjective(m, Max, 10) #doesn't matter we can write constant

# Active Power & Reactive Power non-linear constraints

@NLconstraints(m, begin
        
    V[2]*V[1]*((YR[2,1])*cos(x[1]-θ₁)+(YI[2,1])*sin(x[1]-θ₁))+ V[2]*V[2]*(YR[2,2])+ 
        V[2]*x[5]*((YR[2,3])*cos(x[1]-x[2])+(YI[2,3])*sin(x[1]-x[2])) == constraints[1]
        
    x[5]*V[2]*((YR[3,2])*cos(x[2]-x[1])+(YI[3,2])*sin(x[2]-x[1]))+ x[5]*x[5]*(YR[3,3])+ 
        x[5]*x[6]*((YR[3,4])*cos(x[2]-x[3])+(YI[3,4])*sin(x[2]-x[3])) == constraints[2]
        
    x[6]*x[5]*((YR[4,3])*cos(x[3]-x[2])+(YI[4,3])*sin(x[3]-x[2]))+ x[6]*x[6]*(YR[4,4])+ 
        x[6]*x[7]*((YR[4,5])*cos(x[3]-x[4])+(YI[4,5])*sin(x[3]-x[4])) == constraints[3]
    
    x[7]*x[6]*((YR[5,4])*cos(x[4]-x[3])+(YI[5,4])*sin(x[4]-x[3]))+ x[7]*x[7]*(YR[5,5])+ 
        x[7]*V[1]*((YR[5,1])*cos(x[4]-θ₁)+(YI[5,4])*sin(x[4]-θ₁)) == constraints[4]
        
    x[5]*V[2]*((YR[3,2])*sin(x[2]-x[1])-(YI[3,2])*cos(x[2]-x[1]))- x[5]*x[5]*(YI[3,3])+
        x[5]*x[6]*((YR[3,4])*sin(x[2]-x[3])-(YI[3,4])*cos(x[2]-x[3])) == constraints[5]
    
    x[6]*x[5]*((YR[4,3])*sin(x[3]-x[2])-(YI[4,3])*cos(x[3]-x[2]))- x[6]*x[6]*(YI[4,4])+
        x[6]*x[7]*((YR[4,5])*sin(x[3]-x[4])-(YI[4,5])*cos(x[3]-x[4])) == constraints[6]
        
    x[7]*x[6]*((YR[5,4])*sin(x[4]-x[3])-(YI[5,4])*cos(x[4]-x[3]))- x[7]*x[7]*(YI[5,5])+
        x[7]*V[1]*((YR[5,1])*sin(x[4]-θ₁)-(YI[5,1])*cos(x[4]-θ₁)) == constraints[7]
        
end)

solve(m; suppress_warnings=true)
x_JuMP = getvalue(x);

# Is it the correct soultion? It is possible to have more than 1 local solution
if sum((x_JuMP-x_true).^2) < xtol
    success_JuMP = true
else
    success_JuMP = false
end
    return success_JuMP, x_JuMP
end

(success_JuMP, x_JuMP) = PF_JuMP(feval, V, θ, Y)

Results = DataFrame(x_true = x_true, x_newton = x_newt, x_JuMP = x_JuMP)

@time PF_JuMP(feval, V, θ, Y)

@time PF_newton(feval,Jeval, V, θ, Y) # although it took me ages to do the Jacobian manually

for i = 5:10:45
    
θ_val = i;
s_newt = 0;
s_JuMP = 0;

# perform 100 iterations for each angle range
for k=1:100
    V = 0.9+0.2*rand(5,1) # random number between 0.9-1.1
    θ = (-θ_val+2*θ_val*randn(5,1))*(π/180) # random angle between -5º to 5º
    V[1] = 1 # slack bus voltage magnitude (reference)
    θ[1] = 0 # slack bus voltage angle (reference)
    
    (success_JuMP, x_JuMP) = PF_JuMP(feval, V, θ, Y)
    (success_newt, x_newt) = PF_newton(feval,Jeval, V, θ, Y)
    
    if success_JuMP == true
        s_JuMP = s_JuMP +1;
    end
    if success_newt == true
        s_newt = s_newt +1;
    end
end
    
println("With θ = -$θ_val degrees to $θ_val degrees")
println("Newton success rate is $s_newt")
println("JuMP success rate is $s_JuMP")
println("__________________________________________")    
end
    

# We can also plot the results

angle_amplitude = Int64[]
ss_JuMP = Int64[]
ss_newt = Int64[]

for i = 2:2:90
    
θ_val = i;
s_newt = 0;
s_JuMP = 0;

# perform 100 iterations for each angle range
for k=1:100
    V = 0.9+0.2*rand(5,1) # random number between 0.9-1.1
    θ = (-θ_val+2*θ_val*randn(5,1))*(π/180) # random angle between -5º to 5º
    V[1] = 1 # slack bus voltage magnitude (reference)
    θ[1] = 0 # slack bus voltage angle (reference)
    
    (success_JuMP, x_JuMP) = PF_JuMP(feval, V, θ, Y)
    (success_newt, x_newt) = PF_newton(feval,Jeval, V, θ, Y)
    
    if success_JuMP == true
        s_JuMP = s_JuMP +1
    end
    if success_newt == true
        s_newt = s_newt +1
    end
end
    
push!(angle_amplitude, i)    
push!(ss_JuMP, s_JuMP)
push!(ss_newt, s_newt)
       
end

plot(angle_amplitude,ss_JuMP,label = "JuMP", xlims = (0, 65), title = "Power Flow Convergence")
plot!(angle_amplitude, ss_newt, label = "Newton", xlabel = "Phase Angle Amplitude (degrees)", ylabel = "Convergence (%)")

res = 0.015*ones(6)
react = 0.01*ones(6)
baseMVA = 100.
Zij = zeros(6,6)
Zij = complex(Zij)
Zij[1,4] = (1/baseMVA)*complex(res[1],react[1]); Zij[4,1] = (1/baseMVA)*complex(res[1],react[1]);
Zij[4,5] = (1/baseMVA)*complex(res[2],react[2]); Zij[5,4] = (1/baseMVA)*complex(res[2],react[2]);
Zij[3,5] = (1/baseMVA)*complex(res[3],react[3]); Zij[5,3] = (1/baseMVA)*complex(res[3],react[3]);
Zij[5,6] = (1/baseMVA)*complex(res[4],react[4]); Zij[6,5] = (1/baseMVA)*complex(res[4],react[4]);
Zij[6,2] = (1/baseMVA)*complex(res[5],react[5]); Zij[2,6] = (1/baseMVA)*complex(res[5],react[5]);
Zij[6,4] = (1/baseMVA)*complex(res[6],react[6]); Zij[4,6] = (1/baseMVA)*complex(res[6],react[6]);

Y = zeros(6,6)
Y = complex(Y)

for i=1:6, j=1:6
    if i != j
        if Zij[i,j] != 0
            Y[i,j] = -Zij[i,j]^-1
        else
            Y[i,j] = 0.0
        end
    end
end
for i = 1:6
    Y[i,i] = -sum(Y[i,:])
end
Y = (Y+conj(Y'))/2;
Y

# generator buses starting at bus 1
vconstraints = [1., 1., 1.]
vlow = 0.9
vhigh = 1.1

# start at bus 2
pconstraints = [150, 75, -90, -100, -125]

# start at bus 4
qconstraints = [-30, -35, -50]

# JuMP does not recognize real and imaginary functions so I have to define the Real and Imaginary Parts of Y
YR = real(Y); YI = imag(Y);

# Initialize Ipopt Solver
m = Model(solver=IpoptSolver(print_level=2))

# Define our x variable
@variable(m, x[1:12])

# Set initial guess
for i in 1:6
    setvalue(x[i], 0.0)
end

for i in 7:12
    setvalue(x[i], 1.0)
end

@NLobjective(m, Max, 10) #doesn't matter we can write as constant

# Voltage magnitude and voltage angle linear constraints

# slackbus reference angle
@constraint(m, x[1] == 0)

# generator voltage constraints
@constraint(m, vconstr[i=1:3], x[i+6] == vconstraints[i])
@constraint(m, vmin[i=1:6], x[i+6] >= vlow)
@constraint(m, vhig[i=1:6], x[i+6] <= vhigh)

# Active Power & Reactive Power non-linear constraints

#active power constraints
@NLconstraint(m, pconstr[i=1:5], sum(x[i+7]*x[j+6]*(YR[i+1,j]*cos(x[i+1]-x[j])+
            YI[i+1,j]*sin(x[i+1]-x[j])) for j = 1:6) == pconstraints[i])

#reactive power constraints
@NLconstraint(m, qconstr[i=1:3], sum(x[i+9]*x[j+6]*(YR[i+3,j]*sin(x[i+3]-x[j])-
            YI[i+3,j]*cos(x[i+3]-x[j])) for j = 1:6) == qconstraints[i])

solve(m)
x_PF = getvalue(x)

# All buses starting at bus 1
vlow = 0.9
vhigh = 1.1

#Generators (buses 1-3 - includes slack)

# generators cost coefficients
a = [0.11, 0.085, 0.1225]
b = [5, 1.2, 1.]
c = [150., 600., 335.];

# active power
pgenmin = [10,10,10]
pgenmax = [200,100,300]

# reactive power
qgenmin = [-300,-300,-300]
qgenmax = [300,300,300]

# Loads (buses 4-6)
plconstraints = [-90, -100, -125]
qlconstraints = [-30, -35, -50]

# Lines

l14constraint = 95

# Initialize Ipopt Solver
m = Model(solver=IpoptSolver(print_level=0))

# Define our x variable
@variable(m, x[1:12])

# Set initial guess
for i in 1:6
    setvalue(x[i], 0.0)
end

for i in 7:12
    setvalue(x[i], 1.0)
end

# Define power expression P = f(x)
@NLexpression(m, pi[i=1:3], sum(x[i+6]*x[j+6]*(YR[i,j]*cos(x[i]-x[j])+
            YI[i,j]*sin(x[i]-x[j])) for j = 1:6))

# minimize generation costs
@NLobjective(m, Min, sum(a[i]*pi[i]^2+b[i]*pi[i]+c[i] for i = 1:3))

# slackbus reference angle
@constraint(m, x[1] == 0)

# All buses starting at bus 1
                        
@constraint(m, vmin[i=1:6], x[i+6] >= vlow)
@constraint(m, vhig[i=1:6], x[i+6] <= vhigh)

#Generators (buses 1-3 - includes slack)     
                        
#active power
@NLconstraint(m, pgenminJP[i=1:3], sum(x[i+6]*x[j+6]*(YR[i,j]*cos(x[i]-x[j])+
            YI[i,j]*sin(x[i]-x[j])) for j = 1:6) >= pgenmin[i])
@NLconstraint(m, pgenmaxJP[i=1:3], sum(x[i+6]*x[j+6]*(YR[i,j]*cos(x[i]-x[j])+
            YI[i,j]*sin(x[i]-x[j])) for j = 1:6) <= pgenmax[i])
                                                
#reactive power
@NLconstraint(m, qgenminJP[i=1:3], sum(x[i+6]*x[j+6]*(YR[i,j]*sin(x[i]-x[j])-
            YI[i,j]*cos(x[i]-x[j])) for j = 1:6) >= qgenmin[i])
@NLconstraint(m, qgenmaxJP[i=1:3], sum(x[i+6]*x[j+6]*(YR[i,j]*sin(x[i]-x[j])-
            YI[i,j]*cos(x[i]-x[j])) for j = 1:6) <= qgenmax[i])
                                                                        

#Loads (buses 4-6)  

#active power constraints
@NLconstraint(m, plconstr[i=1:3], sum(x[i+9]*x[j+6]*(YR[i+3,j]*cos(x[i+3]-x[j])+
            YI[i+3,j]*sin(x[i+3]-x[j])) for j = 1:6) == plconstraints[i])

#reactive power constraints
@NLconstraint(m, qlconstr[i=1:3], sum(x[i+9]*x[j+6]*(YR[i+3,j]*sin(x[i+3]-x[j])-
            YI[i+3,j]*cos(x[i+3]-x[j])) for j = 1:6) == qlconstraints[i])

#line constraint

@NLexpression(m, p14, x[1+6]*x[1+6]*YR[1,1] + x[1+6]*x[4+6]*(YR[1,4]*cos(x[1]-x[4])+YI[1,4]*sin(x[1]-x[4])))
@NLexpression(m, q14, -x[1+6]*x[1+6]*YI[1,1] + x[1+6]*x[4+6]*(YR[1,4]*sin(x[1]-x[4])-YI[1,4]*cos(x[1]-x[4])))

# line constraint must be in terms of complex power
@NLexpression(m, s14sq, p14^2+q14^2)                                                                                                
                                                                                       
@NLconstraint(m, l14constrJP, s14sq <= l14constraint^2)
                                                                                                
solve(m)
x_OPF = getvalue(x)

getobjectivevalue(m)

v_θ_jump = x_OPF[1:6]; v_jump = x_OPF[7:12]

P_jp = zeros(6,6); Q_jp = zeros(6,6); 
P_jump = zeros(6); Q_jump = zeros(6); 

for i = 1:6, j = 1:6
    P_jp[i,j] = v_jump[i]*v_jump[j]*(YR[i,j]*cos(v_θ_jump[i]-v_θ_jump[j]) + 
        YI[i,j]*sin(v_θ_jump[i]-v_θ_jump[j]))
    Q_jp[i,j] = v_jump[i]*v_jump[j]*(YR[i,j]*sin(v_θ_jump[i]-v_θ_jump[j]) - 
        YI[i,j]*cos(v_θ_jump[i]-v_θ_jump[j]))
end

for i = 1:6
    P_jump[i] = sum(P_jp[i,:])
    Q_jump[i] = sum(Q_jp[i,:])
end

[P_jump Q_jump]

LMP_jump = zeros(6)

# For loads
LMP_jump[4:6] = abs(getdual(plconstr))

# For generators

μ_pmin = getdual(pgenminJP)
μ_pmax = abs(getdual(pgenmaxJP))

for i = 1:3
    LMP_jump[i] = 2*a[i]*P_jump[i]+b[i]-μ_pmin[i]+μ_pmax[i]
end
    
LMP_jump

JuMP_r = DataFrame(Bus = 1:6, Voltage_pu = v_jump, θ_deg = v_θ_jump*180/π, P_MW = P_jump, Q_MW = Q_jump, 
    LMP_JuMP = LMP_jump)

v_mat = [1.090,1.095,1.100,1.076, 1.078,1.076]; v_θ = [0.000,-0.472,-0.243,-0.407, -0.470,-0.532];
p_mat = [94.69,100.00,125.51,-90.00,-100.00,-125.00]; q_mat = [7.67,58.47,52.33,-30.00, -35.00,-50.00];
LMP_mat = [25.832,32.141,31.751,32.872,32.805,32.960]
Matlab = DataFrame(Bus = 1:6, Voltage_pu = v_mat, θ_deg = v_θ, P_MW = p_mat, Q_MW = q_mat, LMP_Matlab = LMP_mat)

# generator cost functions
f1(x) = 0.25*x^2 + 5*x + 1000
f2(x) = 0.20*x^2 + 7*x + 1100
f3(x) = 0.05*x^2 + 1*x + 100

plot(f1,0,300, label = "Generator 1", xlabel = "Power (MW)", ylabel = "Generation Costs (dollars/hr)",
 xlims = (0, 340),ylims = (-1000,33000))
plot!(f2, label = "Generator 2")
plot!(f3, label = "Generator with Market Power")

function mrkt_power(P_cheap)

# All buses starting at bus 1
vlow = 0.9
vhigh = 1.1

#Generators (buses 1-3 - includes slack)

# generators cost coefficients
a = [0.25, 0.20, 0.05]
b = [5., 7., 1.]
c = [1000., 1100., 100.]

# active power
pgenmin = [10,10,10]
pgenmax = [600,600,P_cheap]

# reactive power
qgenmin = [-300,-300,-300]
qgenmax = [300,300,300]

# Loads (buses 4-6)
plconstraints = [-150, -150, -150]
qlconstraints = [-30, -35, -50]

l14constraint = 95    
    
# Initialize Ipopt Solver
m = Model(solver=IpoptSolver(print_level=0))

# Define our x variable
@variable(m, x[1:12])

# Set initial guess
for i in 1:6
    setvalue(x[i], 0.0)
end

for i in 7:12
    setvalue(x[i], 1.0)
end

# Define power expression P = f(x)
@NLexpression(m, pi[i=1:3], sum(x[i+6]*x[j+6]*(YR[i,j]*cos(x[i]-x[j])+
            YI[i,j]*sin(x[i]-x[j])) for j = 1:6))

# minimize generation costs
@NLobjective(m, Min, sum(a[i]*pi[i]^2+b[i]*pi[i]+c[i] for i = 1:3))

# slackbus reference angle
@constraint(m, x[1] == 0)

# All buses starting at bus 1
                        
@constraint(m, vmin[i=1:6], x[i+6] >= vlow)
@constraint(m, vhig[i=1:6], x[i+6] <= vhigh)

#Generators (buses 1-3 - includes slack)     
                        
#active power
@NLconstraint(m, pgenminJP[i=1:3], sum(x[i+6]*x[j+6]*(YR[i,j]*cos(x[i]-x[j])+
            YI[i,j]*sin(x[i]-x[j])) for j = 1:6) >= pgenmin[i])
@NLconstraint(m, pgenmaxJP[i=1:3], sum(x[i+6]*x[j+6]*(YR[i,j]*cos(x[i]-x[j])+
            YI[i,j]*sin(x[i]-x[j])) for j = 1:6) <= pgenmax[i])
                                                
#reactive power
@NLconstraint(m, qgenminJP[i=1:3], sum(x[i+6]*x[j+6]*(YR[i,j]*sin(x[i]-x[j])-
            YI[i,j]*cos(x[i]-x[j])) for j = 1:6) >= qgenmin[i])
@NLconstraint(m, qgenmaxJP[i=1:3], sum(x[i+6]*x[j+6]*(YR[i,j]*sin(x[i]-x[j])-
            YI[i,j]*cos(x[i]-x[j])) for j = 1:6) <= qgenmax[i])
                                                                        

#Loads (buses 4-6)  

#active power constraints
@NLconstraint(m, plconstr[i=1:3], sum(x[i+9]*x[j+6]*(YR[i+3,j]*cos(x[i+3]-x[j])+
            YI[i+3,j]*sin(x[i+3]-x[j])) for j = 1:6) == plconstraints[i])

#reactive power constraints
@NLconstraint(m, qlconstr[i=1:3], sum(x[i+9]*x[j+6]*(YR[i+3,j]*sin(x[i+3]-x[j])-
            YI[i+3,j]*cos(x[i+3]-x[j])) for j = 1:6) == qlconstraints[i])

#line constraint

@NLexpression(m, p14, x[1+6]*x[1+6]*YR[1,1] + x[1+6]*x[4+6]*(YR[1,4]*cos(x[1]-x[4])+YI[1,4]*sin(x[1]-x[4])))
@NLexpression(m, q14, -x[1+6]*x[1+6]*YI[1,1] + x[1+6]*x[4+6]*(YR[1,4]*sin(x[1]-x[4])-YI[1,4]*cos(x[1]-x[4])))

# line constraint must be in terms of complex power
@NLconstraint(m, s14sq, p14^2+q14^2 <= l14constraint^2)                                                                                                
                                                                                                
solve(m)
x_OPF = getvalue(x)

v_θ_jump = x_OPF[1:6]; v_jump = x_OPF[7:12]

P_jp = zeros(6,6); Q_jp = zeros(6,6); 
P_jump = zeros(6); Q_jump = zeros(6); 

for i = 1:6, j = 1:6
    P_jp[i,j] = v_jump[i]*v_jump[j]*(YR[i,j]*cos(v_θ_jump[i]-v_θ_jump[j]) + 
        YI[i,j]*sin(v_θ_jump[i]-v_θ_jump[j]))
    Q_jp[i,j] = v_jump[i]*v_jump[j]*(YR[i,j]*sin(v_θ_jump[i]-v_θ_jump[j]) - 
        YI[i,j]*cos(v_θ_jump[i]-v_θ_jump[j]))
end

for i = 1:6
    P_jump[i] = sum(P_jp[i,:])
    Q_jump[i] = sum(Q_jp[i,:])
end

LMP_jump = zeros(6)

# For loads
LMP_jump[4:6] = abs(getdual(plconstr))

# For generators

μ_pmin = getdual(pgenminJP)
μ_pmax = abs(getdual(pgenmaxJP))

for i = 1:3
    LMP_jump[i] = 2*a[i]*P_jump[i]+b[i]-μ_pmin[i]+μ_pmax[i]
end
                                                                                                    
obj = getobjectivevalue(m)
                                                                                                    
return P_jump, Q_jump, LMP_jump, obj, x_OPF                                                                                                    
                                                                                                    
end                                                                                               

(P_opt, Q_opt, LMP_opt, optval)  = mrkt_power(600)

gen_p = round(P_opt[3],2)
optval = round(optval,2)

println("The optimal production of the generator with market power is $gen_p")
println("The total cost of generation the $optval dollars/hr")
println("The Deadweight loss is 0")

P_iter = 10:2:328

length = size(P_iter)[1]

P_run = zeros(length,6)
Q_run = zeros(length, 6)
LMP_run = zeros(length, 6)
obj_run = zeros(length)

c = 0

for i in P_iter
    c = c + 1
    (P, Q, LMP, obj)  = mrkt_power(i)
    P_run[c,:] = P
    Q_run[c,:] = Q
    LMP_run[c,:] = LMP
    obj_run[c] = obj
end

# Revenue = PQ = locational marginal prices x production level of generator
rev_gen = LMP_run[:,3].*P_run[:,3]

# Cost = generation cost function evaluated at production level
costs_gen = f3.(P_run[:,3]);

# Profit = Revenue - Costs
profit_gen = rev_gen - costs_gen

# Optimal generation point = profit maximizing generation level
gen_opt_ind = indmax(profit_gen)
gen_opt_P = round(P_run[gen_opt_ind,3], 2)

# DWL = associated system-wide costs - system-wide costs at society's optimal point
gen_opt_DWL = round(obj_run[gen_opt_ind]-optval,2);

plot(P_run[:,3],rev_gen, label = "Revenue generator",xlims = (10, 340),ylims = (-1000,22000),
    xticks = 30:30:330, yticks = 0:5000:25000, xlabel = "Generator Power (MW)", ylabel = "Value (dollars/hour)", 
    w = 3, title = "Market Power of a Generator")
plot!(P_run[:,3],costs_gen, label = "Costs generator", w = 3)
plot!(P_run[:,3],profit_gen, label = "Profit generator", w = 3)
scatter!([gen_opt_P],[profit_gen[gen_opt_ind]],color=[:red],marker=([:d],6,0.8,stroke(3,:gray)), 
    label = "Profit-Maximizing point")
scatter!([gen_opt_P+3], [10_500], series_annotations = ["P = $gen_opt_P MW"], markersize = 0, color = [:white]
    , label = "Annotation")

default(legend=true)
plot(P_run[:,3],obj_run-optval, label = "DWL to Society",xlims = (10, 360),ylims = (-1000,26000),
    xticks = 30:30:330, yticks = 0:5000:25000, xlabel = "Generator Power (MW)", ylabel = "Value (dollars/hour)", 
    w = 3, title = "Market Power of a Generator",fill=(0,:skyblue))
plot!(P_run[:,3],rev_gen-costs_gen, label = "Profit generator", w = 3)
scatter!([P_opt[3]],[0],color=[:green],marker=([:d],6,0.8,stroke(3,:gray)), label = "Society Optimal")
scatter!([gen_opt_P],[gen_opt_DWL],color=[:red],marker=([:d],6,0.8,stroke(3,:gray)), label = "Generator Optimal")
vline!([gen_opt_P, P_opt[3]],color = [:purple], label = "Support")
scatter!([gen_opt_P+3], [14_500], series_annotations = ["DWL = $gen_opt_DWL dollars/hr"], markersize = 0, color = [:white]
, label = "Annotation")


