
using JuMP, KNITRO, DataFrames

# Line data for Garver 6 bus example
Line_ID = 1:15; 
Corridor = ["1-2","1-3","1-4","1-5","1-6","2-3","2-4","2-5","2-6","3-4","3-5","3-6","4-5","4-6","5-6"]
Resistance = [0.10,0.09,0.15,0.05,0.17,0.05,0.10,0.08,0.08,0.15,0.05,0.12,0.16,0.08,0.15]
Reactance = [0.40,0.38,0.60,0.20,0.68,0.20,0.40,0.31,0.30,0.59,0.20,0.48,0.63,0.30,0.61]
Line_cost = [40.,38,60,20,68,20,40,31,30,59,20,48,63,30,61]
Capacity_MW = [100,100,80,100,70,100,100,100,100,82,100,100,75,100,78]

Garver_line = DataFrame(Line_ID = Line_ID, Corridor = Corridor, Resistance = Resistance, 
    Reactance = Reactance, Investment_Cost = Line_cost, Capacity_MW = Capacity_MW)

# Existing Lines Garver
Ex_ID = [1,3,4,6,7,11]
Ex_lines = DataFrame(Line_ID = Line_ID[Ex_ID], Corridor = Corridor[Ex_ID], 
    Resistance = Resistance[Ex_ID], Reactance = Reactance[Ex_ID], Capacity_MW = Capacity_MW[Ex_ID])

# Bus data for Garver 6 bus example with reactive power
Bus_ID = 1:6
PG_min = [0.,0,0,0,0,0]
PG_max = [150.,0,360,0,0,600]
QG_min = [-10.,0,-10,0,0,-10]
QG_max = [65.,0,150,0,0,200]
LMP_G = [10.,0,20,0,0,30]
PD = [80,240,40,160,240,0]
QD = [16, 48, 8, 32, 48,0]

Garver_bus_ac = DataFrame(Bus_ID = Bus_ID, PD_MW = PD, QD_MVAr = QD, PG_min = PG_min, PG_max = PG_max,
    QG_min = QG_min, QG_max = QG_max, LMP_G = LMP_G)

# line data
C1 = [1,1,1,1,1,2,2,2,2,3,3,3,4,4,5]
C2 = [2,3,4,5,6,3,4,5,6,4,5,6,5,6,6]

# resized line data vectors (45-long vectors)
Line_ID_jump = 1:45
Corridor_jump = [Corridor;Corridor;Corridor]
lcost_jump = [Line_cost;Line_cost;Line_cost]
Capacity_MW_jump = [Capacity_MW;Capacity_MW;Capacity_MW]
C1_jump = [C1;C1;C1]
C2_jump = [C2;C2;C2]

# NPV factor for 20 years with 10% discount rate to make fixed and operating costs compareable
df = 0.0
for i = 1:20
    df += 1/(1.1)^i
end
fy = 1/df

# bus data
generators = [1,3,6]
loads = [2,4,5];

# sets of bus and neighboring line IDs
bus_1_c1 = Int64[]; bus_2_c1 = Int64[]; bus_3_c1 = Int64[]; bus_4_c1 = Int64[]; bus_5_c1 = Int64[]; bus_6_c1 = Int64[]
bus_1_c2 = Int64[]; bus_2_c2 = Int64[]; bus_3_c2 = Int64[]; bus_4_c2 = Int64[]; bus_5_c2 = Int64[]; bus_6_c2 = Int64[]

for i = 1:45
    if C1_jump[i] == 1
        bus_1_c1 = push!(bus_1_c1,i)
    elseif C1_jump[i] == 2
        bus_2_c1 = push!(bus_2_c1,i)
    elseif C1_jump[i] == 3
        bus_3_c1 = push!(bus_3_c1,i)
    elseif C1_jump[i] == 4
        bus_4_c1 = push!(bus_4_c1,i)       
    elseif C1_jump[i] == 5
        bus_5_c1 = push!(bus_5_c1,i)
     elseif C1_jump[i] == 6
        bus_6_c1 = push!(bus_6_c1,i)
    end
end   
   
for i = 1:45
    if C2_jump[i] == 1
        bus_1_c2 = push!(bus_1_c2,i)
    elseif C2_jump[i] == 2
        bus_2_c2 = push!(bus_2_c2,i)
    elseif C2_jump[i] == 3
        bus_3_c2 = push!(bus_3_c2,i)
    elseif C2_jump[i] == 4
        bus_4_c2 = push!(bus_4_c2,i)     
    elseif C2_jump[i] == 5
        bus_5_c2 = push!(bus_5_c2,i)
     elseif C2_jump[i] == 6
        bus_6_c2 = push!(bus_6_c2,i)
    end
end  

# Set of lines going from bus s
bus_mat_c1 = [bus_1_c1,bus_2_c1,bus_3_c1,bus_4_c1,bus_5_c1,bus_6_c1]

# Set of lines going into bus s
bus_mat_c2 = [bus_1_c2,bus_2_c2,bus_3_c2,bus_4_c2,bus_5_c2,bus_6_c2]

# line data susceptance and conductance
b = -Reactance./(Reactance.^2+Resistance.^2)
g = Resistance./(Reactance.^2+Resistance.^2)

# resized line data vectors (45-long vectors)
b_jump = 100*[b;b;b]
g_jump = 100*[g;g;g]

# buses
delta_lim = 20*π/180;

# All buses starting at bus 1
vlow = 0.95
vhigh = 1.05

n = Model(solver=KnitroSolver(ms_enable = 1,ms_maxsolves = 500))

@variable(n, 0 <= w[1:45] <= 1) # 45 variables because each line can be built up to 3 times
@variable(n, -delta_lim <= δ[1:6] <= delta_lim)
@variable(n, vlow <= v[1:6] <= vhigh)

@NLexpression(n, pst[i=1:45], w[i]*(v[C1_jump[i]]^2*g_jump[i]-v[C1_jump[i]]*v[C2_jump[i]]*(g_jump[i]*
            cos(δ[C1_jump[i]]-δ[C2_jump[i]])+b_jump[i]*sin(δ[C1_jump[i]]-δ[C2_jump[i]]))))

@NLexpression(n, qst[i=1:45], w[i]*(-v[C1_jump[i]]^2*b_jump[i]+v[C1_jump[i]]*v[C2_jump[i]]*(b_jump[i]*
            cos(δ[C1_jump[i]]-δ[C2_jump[i]])-g_jump[i]*sin(δ[C1_jump[i]]-δ[C2_jump[i]]))))

# Set a constraint for exisiting lines
@constraint(n, exconstr[i=1:6], w[Ex_ID[i]] == 1)

# Set constraint for reference bus
@constraint(n, refbus, δ[1] == 0)

# Calculate an expression for the power flows through each bus
@NLexpression(n, pk[i=1:6],sum(pst[j] for j in bus_mat_c1[i])-sum(pst[j] for j in bus_mat_c2[i]))
@NLexpression(n, qk[i=1:6],sum(qst[j] for j in bus_mat_c1[i])-sum(qst[j] for j in bus_mat_c2[i]))
                                                                               
@NLconstraint(n, load_balance_ac[i in loads] , pk[i]  == -PD[i])
@NLconstraint(n, load_balance_re[i in loads] , qk[i]  == -QD[i])
                                        
@NLconstraint(n, gen_min_ac[i in generators] , pk[i]  >= -PD[i] + PG_min[i])
@NLconstraint(n, gen_min_re[i in generators] , qk[i]  >= -QD[i] + QG_min[i])

@NLconstraint(n, gen_max_ac[i in generators] , pk[i]  <= -PD[i] + PG_max[i])
@NLconstraint(n, gen_max_re[i in generators] , qk[i]  <= -QD[i] + QG_max[i])                                        
                                        
@NLconstraint(n, max_flow_cstr_ac[i=1:45], pst[i]^2 + qst[i]^2 <= w[i]*Capacity_MW_jump[i]^2)
 
@NLconstraint(n, relax[i=1:45],w[i]*(1-w[i]) <= 10^-9)                                        
                                        
# Define the minimization objective
@NLexpression(n, gen_costs_ac, sum(pk[i]*LMP_G[i] for i in generators))
@NLexpression(n, line_costs_ac, sum(w[i]*lcost_jump[i] for i = 1:45))
@NLexpression(n, obj_ac, line_costs_ac + gen_costs_ac*(fy*8760/10^6))
                                                    
@NLobjective(n, Min, obj_ac)

@NLexpression(n,pg[i=1:6],pk[i]+PD[i])
                                                                
solve(n)

W_jump_ac = getvalue(w)
pg_jump_ac = getvalue(pg)
qst_jump_ac = getvalue(qst)
pst_jump_ac = getvalue(pst)
line_costs_jump_ac = getvalue(line_costs_ac)

getvalue(gen_costs_ac)*(fy*8760/10^6)

getvalue(obj_ac)

ac_lines = Int64[]
for i = 1:45
    if W_jump_ac[i] > 0.110
        push!(ac_lines,i)
    end
end

ac_lines

5+5

ac_new_lines = setdiff(ac_lines,Ex_ID)

line_flows = DataFrame(Line_ID = Line_ID_jump[ac_lines], Corridor = Corridor_jump[ac_lines], 
    Active_Power_MW = pst_jump_ac[ac_lines], Reactive_Power_MVAr = qst_jump_ac[ac_lines])

inv_ac = sum(lcost_jump[ac_new_lines])
println("The total investment cost of new lines is $inv_ac M dollars")